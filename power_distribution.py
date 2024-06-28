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
def calc_R_n2(theta, D, beta, n, h, alpha):
    gamma = 1 / np.sqrt(1 - beta ** 2)

    lambda_this = calc_lambda(theta, 0, D, beta, n)

    lambda_e = (lambda_this / (2 * np.pi)) * (
                (beta * gamma) / (1 + beta ** 2 * gamma ** 2 * np.sin(theta) ** 2 * np.sin(0) ** 2))

    l = h / np.tan(alpha)

    k = 2 * np.pi * n / (l * (1 / beta - np.tan(theta)))

    return (2
            * beta ** 2 / (k ** 2 * D ** 2)
            * np.sin(theta - alpha) ** 2
            / ((1 - beta * np.cos(theta)) * (1 - beta * np.cos(theta - 2 * alpha)))
            * np.exp((2 * l - D) / lambda_e)
            * (np.cosh(0 * np.tan(alpha) / lambda_e) - np.cos(
                (1 / beta - (np.cos(theta - alpha) / np.cos(alpha))) * k * 0))
            )


# Calculate angular power distribution of SPR radiation along normal plane (phi = 0)
# theta - Angle 'up' wrt. beam direction (rad)
# phi - Angle 'left/right' wrt. beam direction (rad)
# I - Beam current (A)
# n - Diffraction order
# L - Total length of grating (m)
# D - Grating period (m)
# alpha - Blaze angle of echelle grating
# h - Height of grating (peak to trough)
# E - Beam energy (GeV)
# d - Height of beam above grating (m)
def calc_distribution(theta, I, n, L, D, alpha, h, E, d, S_inc):
    phi = 0

    beta = calc_beta(E)
    gamma = 1 / np.sqrt(1 - beta ** 2)
    lambda_this = calc_lambda(theta, phi, D, beta, n)

    # Elementary charge
    q = 1.602e-19

    lambda_e = lambda_this / (2 * np.pi) * beta * gamma

    R_n2 = calc_R_n2(theta, D, beta, n, h, alpha)

    term2 = (n ** 2 * beta ** 3) / (1 - beta * np.cos(theta)) ** 3

    return 2 * np.pi * q ** 2 * I * L / (D ** 2) * term2 * R_n2 * np.exp(-2 * d / lambda_e) * S_inc


if __name__ == '__main__':
    print('Don\'t use this graph! Make sure the variable values are correct.')

    I = 7e-2  # Beam current (A)
    n = 1  # Diffraction order
    L = 2.5e-2  # Total grating length (m)
    D = 8.33e-7  # Grating period (m)
    E = 1  # Beam energy (GeV)
    d = 1.6e-4  # Height of beam above grating (m)

    # ATTN: Get exact exact values from manufacturer
    alpha = np.pi / 36  # Blaze angle of echelle grating
    h = 0.0575675406  # Height of echelle grating

    # ATTN: I plugged it into Wolfram Alpha and got this
    # I have no idea if this is correct
    # Please send help :)
    S_inc = 1.27e-4

    thetas = np.linspace(0, np.pi / 2, 501)

    power_dist = np.array([calc_distribution(theta, I, n, L, D, alpha, h, E, d, S_inc) for theta in thetas])

    R_ns = np.array([calc_R_n2(theta, D, calc_beta(E), n, h, alpha) for theta in thetas])

    print(R_ns)

    fig, ax = plt.subplots()
    ax.plot(thetas, R_ns)

    ax.set_title("Estimated angular power distribution of SPR along normal plane (DO NOT USE)")

    ax.set_xlabel('Theta (rad)')
    ax.set_ylabel('R_n^2')

    plt.show()
