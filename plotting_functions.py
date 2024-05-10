"""
Functions for plotting data with bempp library.
"""

import bempp.api
import numpy as np
import matplotlib.pyplot as plt
import math

Z_0 = 376.730313668

def plot_directivity(solved_system, div_space, wavenumber, radiated_power=1, number_of_angles=400, plane='YZ'):
    """
    Plot the directivity of a solved efie
    Planes: 'YZ', 'XZ', 'XY'
    """
    number_of_angles = 400
    angles = np.linspace(0, 2*np.pi, number_of_angles)
    if plane in ['YZ', 'ZY']:
        unit_points = np.array([np.zeros(len(angles)),np.sin(angles),np.cos(angles)])
    elif plane in ['XZ', 'ZX']:
        unit_points = np.array([np.sin(angles),np.zeros(len(angles)), np.cos(angles)])
    elif plane in ['XY', 'YX']:
        unit_points = np.array([np.sin(angles),np.cos(angles), np.zeros(len(angles))])
    else:
        raise ValueError(f"Invalid plane: {plane}")

    electric_far = bempp.api.operators.far_field.maxwell.electric_field(div_space, unit_points, wavenumber)
    far_field = electric_far * solved_system

    far_field_squared = far_field*far_field.conj()
    far_field_amplitude = np.real(np.sum(far_field_squared, axis=0))
    radiation_intensity = far_field_amplitude**2 / (2 * Z_0)
    directivity = 10 * np.log10(radiation_intensity * (4 * math.pi / radiated_power))

    fig1, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    plt.rcParams['figure.figsize'] = (10, 8) # Resize the figure
    plt.plot(angles, directivity)
    plt.title(f"Directivity [dB] - {plane}")
    _ = plt.xlabel('Angle (Degrees)')

def plot_impedance(frequencies, impedances):
    """
    Plot the impedance vs frequency
    """
    z_real = [z.real for z in impedances]
    z_imag = [z.imag for z in impedances]
    plt.plot(frequencies, z_real, label='Resistance')
    plt.plot(frequencies, z_imag, label='Reactance')

    # Adding labels and legend
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Impedance (Ohm)')
    plt.title('Impedance vs Frequency')
    plt.legend()

def plot_s_param(frequencies, impedances):
    impedances = np.array(impedances)
    z_0 = 50 # ohms
    gamma = (impedances-z_0)/(impedances+z_0)
    s_param = -20*np.log10(gamma)

    plt.plot(frequencies, s_param)

    plt.xlabel('Frequency (GHz)')
    plt.ylabel('S-Parameter, dB')

def plot_squared_field_density(solved_system, div_space, wavenumber, nx=200, nz=200, v_max=5):
    extent = 0.1
    x, y, z = np.mgrid[-extent:extent:nx * 1j, 0:0:1j, -extent:extent:nz * 1j]
    points = np.vstack((x.ravel(), y.ravel(), z.ravel()))
    el_field = bempp.api.operators.potential.maxwell.electric_field(div_space, points, wavenumber)
    # mag_field = bempp.api.operators.potential.maxwell.magnetic_field(div_space, points, wavenumber)
    e_field_data = -el_field * solved_system
    e_field_density = np.real(np.sqrt(np.sum(e_field_data * e_field_data.conj(), axis=0)))
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    im = ax.imshow(e_field_density.reshape((nx, nz)).T,
                   cmap='coolwarm', origin='lower',
                   extent=[-extent, extent, -extent, extent], vmax=v_max)

    cbar = fig.colorbar(im, ax=ax, format='%.2e')
    cbar.set_label("Electric Field Density (V/m)")
    _ = ax.set_title("Electric Field Density")
    ax.set_xlabel('x, mm')
    ax.set_ylabel('z, mm')