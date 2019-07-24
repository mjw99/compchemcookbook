# AMBER optimises its electrostatic calculation between pairs (EE) by 
# incorporating constants into the charges values within a prmtop file
# See: http://ambermd.org/Questions/units.html
# The charge values are multipled by 18.2223 and this 
# is defined as AMBER_ELECTROSTATIC in SANDER

# However, repeating the EE in other engines, one sees a slight 
# difference against SANDER and PMEMD.


# For example, with the ACE system:
# EE OpenMM 4.4.1 Reference -27.8352952389 kcal/mol
# EE JAMBER                 -27.8342137002 kcal/mol
#                                 ^^^^^^^^^^^^^^^^


# The origin of this lies with the fact that EPSILON_ZERO used by AMBER
# to calculate this is slightly wrong; CHARMM suffers from something similar...

# AMBER 	332.0522173
# CHARMM 	332.0716
# "Exact"	332.063408432

# This script carries out the back 
# calculation to work out this difference:

#  ATKINS's molecular Quantum Mechanics  8.854187817e-12
#  AMBER                                 8.85449607758e-12
#                                             ^^^^^^^^^^^

# Hence AMBER is wrong in the 4th decimal place.

# Refs
# http://www.charmm.org/ubbthreads/ubbthreads.php?ubb=showflat&Number=12674





# AMBER_ELECTROSTATIC
AMBER_CHARGE_PREFACTOR =  18.2223 
# This is  AMBER_CHARGE_PREFACTOR**2
AMBER_PREFACTOR = 332.0522173

KCAL_TO_JOULE = 4184 
JOULES_TO_KILOJOULE = 0.001


## OpenMM
PI = 3.14159265358979323846
E_CHARGE = 1.60217733e-19
AVOGADRO = 6.0221367e23


## ATKINS's molecular Quantum Mechanics
#E_CHARGE     = 1.602176e-19
EPSILON_ZERO = 8.854187817e-12
#AVOGADRO     = 6.02214e23

## SANDER Constants
#PI = 3.14159265358979323846
#E_CHARGE = 1.60217733e-19
#AVOGADRO = 6.0221367e23

AMBER_EPSILON_ZERO = 8.85449607758e-12




print("AMBER's prefactor ", AMBER_PREFACTOR )
# Exact
# 332.063408432
print("'Exact' prefactor ", str((( ( 1 * E_CHARGE  * E_CHARGE) / ( 4 * PI * EPSILON_ZERO * 1E-10 ) ) * AVOGADRO) / KCAL_TO_JOULE))




# OpenMM's 
# ./OpenMM4.1.1-Source/platforms/reference/src/SimTKUtilities/SimTKOpenMMRealType.h
# ONE_4PI_EPS0 138.935456
print("OpenMM's ONE_4PI_EPS0 is ", str(138.935456))

print("OpenMM's ONE_4PI_EPS0 should be " ,   str((( ( 1 * E_CHARGE  * E_CHARGE) / ( 4 * PI * EPSILON_ZERO * 1E-9 ) ) * AVOGADRO) * JOULES_TO_KILOJOULE))


# 332.0522173
# What is AMBER Using for EPSILON_ZERO?
print("Back calculating AMBER_EPSILON_ZERO ", str((( ( 1 * E_CHARGE  * E_CHARGE) / ( 4 * PI * AMBER_PREFACTOR * 1E-10 ) ) * AVOGADRO) / KCAL_TO_JOULE))


# Hence, what would OpenMM's ONE_4PI_EPS0's be, using AMBER_EPSILON_ZERO ?
print("Broken ONE_4PI_EPS0 would be " , str( (( ( 1 * E_CHARGE  * E_CHARGE) / ( 4 * PI * AMBER_EPSILON_ZERO * 1E-9 ) ) * AVOGADRO) * JOULES_TO_KILOJOULE))
