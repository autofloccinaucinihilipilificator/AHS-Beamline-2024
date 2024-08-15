# Calculates, graphs angular power distribution of emitted SPR


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter # smooth data


# Get velocity as a fraction of c from electron energy in GeV
def calc_beta(E):
    # c = 1  # Speed of light
    m = 5.11e-4  # Electron mass in GeV
    # return np.sqrt(1 - (m**2 * c**4) / E**2)
    return np.sqrt(1 - (m ** 2) / E ** 2)


def calc_lambda(theta, phi, D, beta, n=-1):
    return (D / np.abs(n)) * (1 / beta - np.cos(theta) * np.sin(phi))


def calc_lambda_e(beta, gamma, theta, phi):
    term1 = calc_lambda(theta, phi, D, beta, n) / (2 * np.pi)
    term2 = (beta * gamma) / (np.sqrt(1 + beta**2 * gamma**2 * np.sin(theta)**2 * np.sin(phi)**2))
    return term1 * term2


# Calculate R_n when phi = 0
# theta - Angle 'up' wrt. beam
# D - echelle grating period
# beta - velocity as fraction of c
# n - diffraction order
# h - height of diffraction grating
# alpha - blaze angle
# lambda_e - evanescent wavelength
# Equations from https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.8.091301
def calc_R_n2(theta, beta, n, h, alpha):
    l = h / np.tan(alpha)
    k = 2 * np.pi * n / (l * (1 / beta - np.cos(theta)))

    # Equation (7)
    term1 = (1 / (np.sin(theta) - np.tan(alpha) * (1 + np.cos(theta))))**2
    term2 = np.tan(alpha)**2 * np.sin(theta)**2
    term3 = np.sinc(k * h * np.sin(theta) / (2 * np.pi))**2
    return term1 * term2 * term3


# Calculate angular power distribution of SPR radiation along normal plane (phi = 0)
# theta - Angle 'up' wrt. beam direction (rad)
# phi - Angle 'left/right' wrt. beam direction (rad)
# I - Beam current (A)
# n - Diffraction order
# L - Total length of grating (m)
# D - Grating period (m)
# alpha - Blaze angle of echelle grating (rad)
# h - Height of grating (peak to trough)
# E - Beam energy (GeV)
# d - Height of beam above grating (m)
def calc_distribution(theta, n, L, D, alpha, h, E, d):
    beta = calc_beta(E)
    gamma = 1 / np.sqrt(1 - beta**2)

    lambda_this = calc_lambda(theta, 0, D, beta, n)

    h_int = lambda_this * beta * gamma / (4 * np.pi)

    term1 = 1 / 137 * np.abs(n) * L / D
    term2 = np.sin(theta)**2 / (1 / beta - np.cos(theta))**2
    term3 = calc_R_n2(theta, beta, n, h, alpha)
    term4 = np.exp(d / h_int)

    return term1 * term2 * term3 * term4


if __name__ == '__main__':

    # n = 1           # Diffraction order
    L = 5.74e-1     # Total grating length (m)
    D = 6.58e-4     # Grating period (m)
    E = 0.855       # Beam energy (GeV)
    d = 1.27e-4      # Height of beam above grating (m)
    alpha = np.pi / 36  # Blaze angle of echelle grating
    h = 0.0575675406  # Height of echelle grating

    beta = calc_beta(E)

    thetas = np.linspace(0, np.pi, 10000)

    power_dist_1 = np.array([calc_distribution(theta, 1, L, D, alpha, h, E, d) for theta in thetas])
    power_dist_2 = np.array([calc_distribution(theta, 2, L, D, alpha, h, E, d) for theta in thetas])
    power_dist_3 = np.array([calc_distribution(theta, 3, L, D, alpha, h, E, d) for theta in thetas])

    thetas = thetas * 180 / np.pi

    fig, ax = plt.subplots()
    ax.plot(thetas, power_dist_1)
    ax.plot(thetas, power_dist_2)
    ax.plot(thetas, power_dist_3)

    ax.set_title("Expected angular distribution of \nSPR per electron along normal plane")

    ax.set_xlabel('$\\theta$ (°)')
    ax.set_ylabel('$\\frac{dN}{d\\Omega}$')

    plt.yscale("log")

    plt.show()