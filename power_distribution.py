# Calculates, graphs angular power distribution of emitted SPR


import numpy as np
import matplotlib.pyplot as plt


# Get velocity as a fraction of c from electron energy in GeV
def calc_beta(E):
    # c = 1  # Speed of light
    m = 5.11e-4  # Electron mass in GeV
    # return np.sqrt(1 - (m**2 * c**4) / E**2)
    return np.sqrt(1 - (m**2) / E**2)


def calc_lambda(theta, phi, D, beam_energy, n=-1):
    beta = calc_beta(beam_energy)
    return (D / np.abs(n)) * (1 / beta - np.cos(theta) * np.sin(phi))


def calc_distribution(I, n, L, D, eps_0, theta, phi, beam_energy, d):
    beta = calc_beta(beam_energy)
    gamma = 1 / np.sqrt(1 - beta**2)
    lambda_this = calc_lambda(theta, phi, D, beam_energy, n)
    h_int = lambda_this * beta * gamma / (4 * np.pi)

    R_n = 'some function of theta and phi idk'

    term1 = (I * n**2 * L) / (2 * D**2 * eps_0)
    term2 = (np.sin(theta)**2 + np.sin(phi)**2) / (((1 / beta) - np.cos(theta) * np.sin*phi)**3)
    term3 = R_n**2
    term4 = np.exp((-1 * d / h_int) * np.sqrt(1 + (beta * gamma * np.cos(phi))**2))
    return term1 * term2 * term3 * term4


if __name__ == '__main__':
    print('Hello world!')
