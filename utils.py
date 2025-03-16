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

# для горизонтального джета dev около 0.75, для вертикального  0.35
def dev(p1, p2, dp2):
    weights = 1 / (np.linalg.norm(dp2, axis=1) ** 2)
    return np.sum(weights * dist(p1, p2)) / len(p2)

# data = readData('src/data_old15kHz.txt')
# print(min(data[4] ** 2 + data[2] ** 2))

