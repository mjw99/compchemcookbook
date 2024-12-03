from openmm.app import *
from openmm import *
from openmm.unit import *
from openmmtools.testsystems import LennardJonesFluid


# Simulation settings
pressure = 80*atmospheres
temperature = 120*kelvin
collision_rate = 5/picoseconds
timestep = 2.5*femtoseconds


# Create a Lennard Jones test fluid
sigma = 3.4*angstrom
epsilon = 0.238 * kilocalories_per_mole
fluid = LennardJonesFluid(sigma=sigma, epsilon=epsilon)
[topology, system, positions] = [fluid.topology, fluid.system, fluid.positions]

# Add a barostat
barostat = MonteCarloBarostat(pressure, temperature)
system.addForce(barostat)

# Retrieve the NonbondedForce
forces = { force.__class__.__name__ : force for force in system.getForces() }
nbforce = forces['NonbondedForce']

# Add a CustomNonbondedForce to handle only alchemically-modified interactions

# Make two sets of particles, one that contains just the particle we will alchemically annihilate
# and the other which contains all the other particles.
alchemical_particles = set([0])
chemical_particles = set(range(system.getNumParticles())) - alchemical_particles


# Define the energy function for the CustomNonbondedForce
# when lambda is 1.0 it is a normal LJ potential, when lambda is 0.0 the interaction vanishes 
energy_function = 'lambda*4*epsilon*x*(x-1.0); x = (sigma/reff_sterics)^6;'
energy_function += 'reff_sterics = sigma*(0.5*(1.0-lambda) + (r/sigma)^6)^(1/6);'
energy_function += 'sigma = 0.5*(sigma1+sigma2); epsilon = sqrt(epsilon1*epsilon2);'
custom_force = CustomNonbondedForce(energy_function)

# Set the custom force to use a cutoff
custom_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
# and set the cutoff distance to be the same as the existing nbforce
custom_force.setCutoffDistance(nbforce.getCutoffDistance())

# Add lambda as a parameter we can change during the simulation
custom_force.addGlobalParameter('lambda', 1.0)

# set the values of sigma and epsilon by copying them from the existing NonBondedForce
custom_force.addPerParticleParameter('sigma')
custom_force.addPerParticleParameter('epsilon')
for index in range(system.getNumParticles()):
    [charge, sigma, epsilon] = nbforce.getParticleParameters(index)
    custom_force.addParticle([sigma, epsilon])
    if index in alchemical_particles:
        # remove the alchemical particle from the existing NonBondedForce
        nbforce.setParticleParameters(index, charge*0, sigma, epsilon*0)

# Set the custom force to occur between just the alchemical particle and the other particles
custom_force.addInteractionGroup(alchemical_particles, chemical_particles)
system.addForce(custom_force)

# Create an integrator
integrator = LangevinIntegrator(temperature, collision_rate, timestep)

# Create a simulation
simulation = Simulation(topology, system, integrator)
simulation.context.setPositions(positions)

# Minimize energy
print('Minimizing energy...')
#LocalEnergyMinimizer.minimize(context)
simulation.minimizeEnergy()


# Collect data

# number of steps per sample
nsteps = 2500

# number of samples to collect per alchemical state
niterations = 500

import numpy as np
lambdas = np.linspace(1.0, 0.0, 10) # alchemical lambda schedule
nstates = len(lambdas)
u_kln = np.zeros([nstates,nstates,niterations], np.float64)
kT = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB * integrator.getTemperature()
for k in range(nstates):
    for iteration in range(niterations):
        print('state %5d iteration %5d / %5d' % (k, iteration, niterations))
        # Set alchemical state
        simulation.context.setParameter('lambda', lambdas[k])
        # Run some dynamics
        simulation.step(nsteps)
        # Compute energies at all alchemical states
        for l in range(nstates):
            simulation.context.setParameter('lambda', lambdas[l])
            u_kln[k,l,iteration] = simulation.context.getState(getEnergy=True).getPotentialEnergy() / kT


# Estimate free energy of Lennard-Jones particle insertion
from pymbar import MBAR, timeseries

## Subsample data to extract uncorrelated equilibrium timeseries
N_k = np.zeros([nstates], np.int32) # number of uncorrelated samples
for k in range(nstates):
    [nequil, g, Neff_max] = timeseries.detect_equilibration(u_kln[k,k,:])
    indices = timeseries.subsample_correlated_data(u_kln[k,k,:], g=g)
    N_k[k] = len(indices)
    u_kln[k,:,0:N_k[k]] = u_kln[k,:,indices].T

# Compute free energy differences
mbar = MBAR(u_kln, N_k)

# If this fails try setting compute_uncertainty to false
# See this issue: https://github.com/choderalab/pymbar/issues/419
results = mbar.compute_free_energy_differences(compute_uncertainty=True)

print("Free energy change to insert a particle = ", results['Delta_f'][nstates-1,0])
print("Statistical uncertainty = ", results['dDelta_f'][nstates-1,0])

with open("data.txt", "a") as file:
    file.write(str(results['Delta_f'][nstates-1,0]) + " " + str(results['dDelta_f'][nstates-1,0]) + "\n")
