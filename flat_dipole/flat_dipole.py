import bempp.api
import numpy as np
import matplotlib
from matplotlib import pylab as plt
matplotlib.use("TkAgg")
import time
import math
from scipy import integrate
import sys, os

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

# get edge
source_edges = []
for edges, dofs in zip(grid.element_edges.T, div_space.local2global):
    for idx, edge_id in enumerate(edges):
        vertex1 = grid.vertices[:, grid.edges[0,edge_id]]
        vertex2 = grid.vertices[:, grid.edges[1,edge_id]]
        if math.isclose(vertex1[2], 0, abs_tol=1e-4) and math.isclose(vertex2[2], 0, abs_tol=1e-4):
            if dofs[idx] not in source_edges:
                source_edges.append(dofs[idx])

assert source_edges
coefficients=np.zeros(div_space.global_dof_count, dtype = 'complex')
for source_edge in source_edges:
    coefficients[source_edge] = voltage0

trace_fun = bempp.api.GridFunction(div_space, coefficients= coefficients, dual_space=curl_space)


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


    current_raw = solved_system.coefficients[source_edges[0]]
    impedance = voltage0 / current_raw
    impedances.append(impedance)
    currents_raw.append(current_raw)


plotting_functions.plot_impedance(frequencies, impedances)
plotting_functions.plot_directivity(solved_system, div_space, wavenumber, plane="YZ")
plotting_functions.plot_directivity(solved_system, div_space, wavenumber, plane="XZ")
plotting_functions.plot_directivity(solved_system, div_space, wavenumber, plane="XY")
plt.show()