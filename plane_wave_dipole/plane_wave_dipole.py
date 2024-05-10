import bempp.api
import numpy as np
import matplotlib
from matplotlib import pylab as plt
matplotlib.use("TkAgg")
import time
import math
from scipy import integrate
import sys, os

sys.path.append(os.path.dirname(__file__) + "/..")

import plotting_functions

current = os.path.dirname(__file__)
time1 = time.time()
GEOM_FACTOR = 1e-3
EPS_0 = 8.854187817E-12
MU_0 = 4 * np.pi * 1E-7
Z_0 = np.sqrt( MU_0 / EPS_0) # impedance of free space

omega = 2 * np.pi * 2.3e9
wavenumber = omega * np.sqrt(EPS_0 * MU_0)

# physical parameters
voltage0=complex(1,0)
impedance0 = 50

gap_size = 5 * GEOM_FACTOR
dipole_length = 55 * GEOM_FACTOR
dipole_diameter = 4 * GEOM_FACTOR
mesh_size = 2 * GEOM_FACTOR

grid = bempp.api.import_grid(current+"/flat_dipole_mesh.msh")
div_space = bempp.api.function_space(grid, "RWG", 0)
curl_space = bempp.api.function_space(grid, "SNC", 0)

direction = np.array([1, 0, 0])
polarization = np.array([0, 0, 1])

@bempp.api.complex_callable(jit=False)
def tangential_trace(x, n, domain_index, result):
    value = np.array(polarization * np.exp(-1j * wavenumber * np.dot(x, direction)))
    result[:] = np.cross(value, n, axis=0)

trace_fun = bempp.api.GridFunction(div_space, fun=tangential_trace, dual_space=curl_space)


impedances = []
currents = []
currents_raw = []
frequencies = np.arange(2.0, 3.0, 0.1)
for frequency in frequencies:
    omega = 2 * np.pi * frequency * 1e9
    wavenumber = omega * np.sqrt(EPS_0 * MU_0)
    efie = bempp.api.operators.boundary.maxwell.electric_field(
        div_space, div_space, curl_space, wavenumber)

    #Solve System
    solved_system = bempp.api.linalg.lu(efie, trace_fun)


    print(max(solved_system.coefficients))

plotting_functions.plot_squared_field_density(solved_system, div_space, wavenumber, v_max=5, extent=0.05)
# plotting_functions.plot_impedance(frequencies, impedances)
plotting_functions.plot_directivity(solved_system, div_space, wavenumber, plane="YZ")
plotting_functions.plot_directivity(solved_system, div_space, wavenumber, plane="XZ")
plotting_functions.plot_directivity(solved_system, div_space, wavenumber, plane="XY")
plt.show()