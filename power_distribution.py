# Calculates, graphs angular power distribution of emitted SPR


import numpy as np
import matplotlib.pyplot as plt
# from scipy.signal import savgol_filter # smooth data


# Get velocity as a fraction of c from electron energy in GeV
def calc_beta(E):
    # c = 1  # Speed of light
    m = 5.11e-4  # Electron mass in GeV
    # return np.sqrt(1 - (m**2 * c**4) / E**2)
    return np.sqrt(1 - (m ** 2) / E ** 2)


def calc_lambda(theta, phi, D, beta, n):
    return (D / n) * (1 / beta - np.cos(theta) * np.sin(np.abs(np.pi / 2 - phi)))


def calc_lambda_e(theta, phi, beta, gamma, lambda_this):
    term1 = lambda_this / (2 * np.pi)
    term2 = (beta * gamma) / np.sqrt(1 + beta**2 * gamma**2 * np.sin(theta)**2 * np.sin(phi)**2)
    return term1 * term2


# Calculate R_n when phi = 0
# theta - Angle 'up' wrt. beam
# phi - Angle 'left
# Equations from https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.8.091301
def calc_R2(theta, phi, beta, gamma, N, L, n, alpha):

    D = L / N # grating period
    h = D * np.tan(alpha) # grating height

    lambda_this = calc_lambda(theta, phi, D, beta, n)
    lambda_e = calc_lambda_e(theta, phi, beta, gamma, lambda_this)

    k = 2 * np.pi * D / lambda_this # wavenumber
    k_x = k * np.cos(phi) * np.sin(theta)
    k_y = -1 * k * np.sin(phi)
    k_z = k * np.cos(phi) * np.cos(theta)

    D_j = k / beta - k_z - k_x * np.tan(alpha) - 1j * np.tan(alpha) / lambda_e

    G_bar_term_1 = N * np.array([np.tan(alpha), 2j * k_y * lambda_e * np.tan(alpha), 1])
    G_bar_term_2 = np.exp((1 / lambda_e - 1j * k_x) * h + 1j * (k / beta - k_z) * D)
    G_bar_term_3 = (np.exp(-1j * D_j * D) - 1) / (1j * D_j * D)

    G_bar = G_bar_term_1 * G_bar_term_2 * G_bar_term_3

    eps_bar_par = np.array([np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -1 * np.sin(theta)])
    eps_bar_perp = np.array([-1 * np.sin(phi), np.cos(phi), 0])

    R2_par = np.abs(np.dot(eps_bar_par, G_bar))**2
    R2_perp = np.abs(np.dot(eps_bar_perp, G_bar))**2

    return R2_par + R2_perp

    # Old
    # l = h / np.tan(alpha)
    # k = 2 * np.pi * n / (l * (1 / beta - np.cos(theta)))
    #
    # # Equation (7)
    # term1 = (1 / (np.sin(theta) - np.tan(alpha) * (1 + np.cos(theta)))) ** 2
    # term2 = np.tan(alpha) ** 2 * np.sin(theta) ** 2
    # term3 = np.sinc(k * h * np.sin(theta) / (2 * np.pi)) ** 2
    # return term1 * term2 * term3


# Calculate angular power distribution of SPR radiation along normal plane (phi = 0)
# theta - Angle 'up' wrt. beam direction (rad)
# phi - Angle 'left/right' wrt. beam direction - pi/2 is parallel to beam (rad)
# I - Beam current (A)
# n - Diffraction order
# L - Total length of grating (m)
# N - Number of grating periods
# alpha - Blaze angle of echelle grating (rad)
# h - Height of grating (peak to trough)
# E - Beam energy (GeV)
# d - Height of beam above grating (m)
def calc_distribution(theta, phi, n, L, N, alpha, h, E, d):
    D = L / N # grating period

    beta = calc_beta(E)
    gamma = 1 / np.sqrt(1 - beta**2)

    lambda_this = calc_lambda(theta, np.pi / 2, D, beta, n)

    h_int = lambda_this * beta * gamma / (4 * np.pi)

    term1 = 1 / 137 * np.abs(n) * L / D
    term2 = np.sin(theta)**2 * np.sin(np.pi / 2 - phi)**2 / (1 / beta - np.cos(theta) * np.sin(np.abs(np.pi / 2 - phi)))**2
    term3 = calc_R2(theta, phi, beta, gamma, N, L, n, alpha)
    term4 = np.exp(-1 * d / h_int * np.sqrt(1 + (beta * gamma * np.cos(np.pi / 2 - phi))**2))

    return term1 * term2 * term3 * term4


if __name__ == '__main__':

    # n = 1             # Diffraction order
    L = 2.5e-2          # Total grating length (m)
    N = 30000           # Total number of grating periods
    E = 0.855           # Beam energy (GeV)
    d = 1.6e-3          # Height of beam above grating (m)
    alpha = np.pi / 36  # Blaze angle of echelle grating (rad)
    h = 6.58e-4         # Height of echelle grating (m)

    beta = calc_beta(E)

    thetas = np.linspace(0, np.pi, 1000)

    phi = 0
    gamma = 1 / np.sqrt(1 - beta**2)

    # power_dist_1 = np.array([calc_distribution(theta, phi, 1, L, N, alpha, h, E, d) for theta in thetas])
    # power_dist_2 = np.array([calc_distribution(theta, phi,2, L, N, alpha, h, E, d) for theta in thetas])
    # power_dist_3 = np.array([calc_distribution(theta, phi, 3, L, N, alpha, h, E, d) for theta in thetas])


    R_ns_1 = np.array([calc_R2(theta, phi, beta, gamma, N, L, 1, alpha) for theta in thetas])
    R_ns_2 = np.array([calc_R2(theta, phi, beta, gamma, N, L, 2, alpha) for theta in thetas])
    R_ns_3 = np.array([calc_R2(theta, phi, beta, gamma, N, L, 3, alpha) for theta in thetas])

    thetas = thetas * 180 / np.pi

    fig, ax = plt.subplots()
    # ax.plot(thetas, power_dist_3)
    # ax.plot(thetas, power_dist_2)
    # ax.plot(thetas, power_dist_1)

    ax.plot(thetas, R_ns_3)
    ax.plot(thetas, R_ns_2)
    ax.plot(thetas, R_ns_1)

    ax.set_title("Expected angular distribution of \nSPR per electron along normal plane")

    ax.set_xlabel('$\\theta$ (°)')
    ax.set_ylabel('$\\frac{dN}{d\\Omega}$')

    plt.yscale("log")

    plt.show()