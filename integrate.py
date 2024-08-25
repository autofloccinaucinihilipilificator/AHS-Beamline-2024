# Using SPR distribution data, integrate numerically over surface of detector

import json
import numpy as np


# Convert coordinates in theta and phi to Cartesian coordinates
def to_cartesian(theta, phi):
    return np.array([np.cos(phi) * np.sin(theta), np.sin(phi), np.cos(phi) * np.cos(theta)])


# Get distance between two points in 3D space
def dist_3d(A, B):
    return np.sqrt(np.sum((A - B)**2))


def get_solid_angle(theta):
    return np.abs(np.sin(theta + np.pi / 1000) - np.sin(theta - np.pi / 1000)) * np.pi / 30 / 500


if __name__ == '__main__':

    filepath = input('Enter data file path: ')

    with open(filepath, 'r') as f:
        data = json.loads(f.read())

    detector_center_theta = float(input('Enter theta of center of detector: '))
    detector_center_phi = float(input('Enter phi of center of detector: '))
    detector_distance = float(input('Enter distance between detector and grating: '))
    detector_radius = float(input('Enter radius of detector: '))

    phis, thetas = np.meshgrid(np.linspace(-np.pi / 60, np.pi / 60, 501), np.linspace(0, np.pi, 501))

    detector_center_coords = np.sqrt(detector_distance**2 - detector_radius**2) / detector_distance * to_cartesian(detector_center_theta, detector_center_phi)

    result = 0

    for i in range(501):
        for j in range(501):
            theta = i * np.pi / 500
            phi = (j - 250) * np.pi / 60 / 250

            coord_current = to_cartesian(theta, phi)

            if dist_3d(coord_current, detector_center_coords) <= detector_radius / detector_distance:
                result += data[i][j] * get_solid_angle(theta)

    print(result)
