# Calculates, graphs angular power distribution of emitted SPR


import numpy as np
import matplotlib.pyplot as plt


# Get velocity as a fraction of c from electron energy in GeV
def calc_beta(E):
    # c = 1  # Speed of light
    m = 5.11e-4  # Electron mass in GeV
    # return np.sqrt(1 - (m**2 * c**4) / E**2)
    return np.sqrt(1 - (m ** 2) / E ** 2)


def calc_lambda(theta, phi, D, beta, n=-1):
    return (D / np.abs(n)) * (1 / beta - np.cos(theta) * np.sin(phi))


# Calculate R_n when phi = 0
# theta - Angle 'up' wrt. beam
# D - echelette grating period
# beta - velocity as fraction of c
# n - diffraction order
# h - height of diffraction grating
# alpha - blaze angle
# Equations from https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.8.091301
def calc_R_n2(theta, D, beta, n, h, alpha):
    k = 2 * np.pi * n / (1 / beta - np.cos(theta))
    l = h / np.tan(alpha)
    gamma = 1 / np.sqrt(1 - beta**2)
    lambda_this = calc_lambda(theta, 0, D, beta, n)
    lambda_e = (lambda_this * beta * gamma) / (2 * np.pi * np.sqrt(1 + beta**2 * gamma**2 * np.sin(theta)**2 * np.sin(0)**2))

    # Equation (7)
    term1 = (1 / (np.sin(theta) - np.tan(alpha) * (1 + np.cos(theta))))**2
    term2 = np.tan(alpha)**2 * np.sin(theta)**2
    term3 = np.sinc(k * h * np.sin(theta) / (2 * np.pi))**2
    return term1 * term2 * term3

    # Equation (6) (equivalent to (7)
    # term1 = (1 / (np.sin(theta) - np.tan(alpha) * (1 + np.cos(theta))))**2
    # term2 = (2 * np.exp(-1 * k * h / (beta * gamma))
    #          * (np.cosh(k * h / (beta * gamma)) - np.cos(k * h * np.sin(theta)))
    #          / (k**2 * l**2))
    # return term1 * term2

    # Equation (5) (different
    # term1 = 2 * beta**2 / (k**2 * D**2)
    # term2 = (np.sin(theta - alpha)**2) / ((1 - beta * np.cos(theta)) * (1 - beta * np.cos(theta - 2 * alpha)))
    # term3 = np.exp(2 * h / lambda_e)
    # term4 = np.cosh(D * np.tan(alpha) / lambda_e) - np.cos(D * k * (1 / beta - np.cos(theta - alpha) / np.cos(alpha)))
    # return term1 * term2 * term3 * term4


# Calculate angular power distribution of SPR radiation along normal plane (phi = 0)
# theta - Angle 'up' wrt. beam direction (rad)
# phi - Angle 'left/right' wrt. beam direction (rad)
# I - Beam current (A)
# n - Diffraction order
# L - Total length of grating (m)
# D - Grating period (m)
# alpha - Blaze angle of echelette grating
# h - Height of grating (peak to trough)
# E - Beam energy (GeV)
# d - Height of beam above grating (m)
def calc_distribution(theta, I, n, L, D, alpha, h, E, d, S_inc):
    pass


if __name__ == '__main__':


    I = 7e-2        # Beam current (A)
    n = 1           # Diffraction order
    L = 1e-2        # Total grating length (m)
    D = 6.58e-4     # Grating period (m)
    # E = 1         # Beam energy (GeV)
    d = 1.6e-3      # Height of beam above grating (m)

    # ATTN: Get exact exact values from manufacturer
    alpha = np.pi / 36  # Blaze angle of echelle grating
    h = 0.0575675406  # Height of echelle grating

    # ATTN: I plugged it into Wolfram Alpha and got this
    # I have no idea if this is correct
    # Please send help :)
    S_inc = 1.27e-4

    # temp
    beta = np.sqrt(1 - 1 / 1000**2)

    thetas = np.linspace(0, np.pi / 2, 501)

    R_ns = np.array([calc_R_n2(theta, D, beta, n, h, alpha) for theta in thetas])

    print(R_ns)

    fig, ax = plt.subplots()
    ax.plot(thetas, R_ns)

    ax.set_title("Estimated angular power distribution of \nSPR of diffraction order 1 along normal plane (DO NOT USE)")

    ax.set_xlabel('Theta (rad)')
    ax.set_ylabel('R^2')

    plt.yscale("log")

    plt.show()