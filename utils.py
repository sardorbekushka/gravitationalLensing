import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist

def readData(filename='src/data_15kHz.txt'):
    ## в new 12 точек (0,0), в old 16
    with open(filename) as f:
        lines = f.readlines()
    data = np.transpose([list(map(float, line.split())) for line in lines])

    return data.T[data[3] != 0].T
    # return data

def dist(p1, p2):
    return np.min(cdist(p2, p1) ** 2, axis=1)

def dev(p1, p2, dp2):
    weights = 1 / (np.linalg.norm(dp2, axis=1) ** 2)
    return np.sum(weights * dist(p1, p2)) / len(p2)

def flux_dev(p1, p2, m1, m2, dm2):
    weights = 1 / dm2 ** 2
    nearest_indices = np.argmin(cdist(p2, p1) ** 2, axis=1)

    return np.sum(weights * (m2 - m1[nearest_indices]) ** 2) / len(m2)

