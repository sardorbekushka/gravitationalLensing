import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist

def readData(filename='src/data_15kHz.txt', cut=True, year=2016):
    ## в new 12 точек (0,0), в old 16
    with open(filename) as f:
        lines = f.readlines()
    data = np.transpose([list(map(float, line.split())) for line in lines])
    data = data.T[data[0] >= year].T

    return data.T[data[3] != 0].T if cut else data
    # return data

def dist(p1, p2):
    return np.min(cdist(p2, p1) ** 2, axis=1)

def dev(p1, p2, dp2):
    weights = 1 / (np.linalg.norm(dp2, axis=1) ** 2)
    return np.sqrt((np.sum(weights * dist(p1, p2)) / len(p2)))

def flux_dev(p1, p2, m1, m2, dm2):
    weights = 1 / dm2 ** 2
    nearest_indices = np.argmin(cdist(p2, p1) ** 2, axis=1)

    return np.sum(weights * (m2 - m1[nearest_indices]) ** 2) / len(m2)

def huber_loss(e, c=1.345):
    return np.where(e <= c, e**2/2, c *e - c**2/2)

def robust_fit(p1, p2, dp2):
    weights = 1 / (np.linalg.norm(dp2, axis=1))
    e = dist(p1, p2) ** 0.5 * weights   # Нормированные отклонения
    return np.sum(huber_loss(e)) / len(e)

def log_cosh(p1, p2, dp2):
    weights = 1 / (np.linalg.norm(dp2, axis=1))
    e = dist(p1, p2) ** 0.5 * weights  # Нормированные отклонения
    return np.sum(np.log(np.cosh(e))) / len(e)

#
# from scipy.interpolate import interp1d
#
# x = [1, 2, 5, 6]
# y = [1, 4, 3, 5]
# #
# # f = interp1d(x, y, kind='cubic')
# # x_ = np.linspace(1, 6, 20)
# # plt.plot(x_, f(x_))
# # plt.scatter(x, y)
# # plt.show()



