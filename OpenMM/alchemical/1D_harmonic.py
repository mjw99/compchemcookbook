from pymbar import testsystems
from pymbar import MBAR

[x_n, u_kn, N_k, s_n] = testsystems.HarmonicOscillatorsTestCase().sample()

mbar = MBAR(u_kn, N_k)

results = mbar.compute_free_energy_differences()

print(results)
