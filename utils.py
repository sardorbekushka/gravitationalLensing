import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
from settings import model
from scipy.interpolate import interp1d


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

def D2z_convertor(z1=0, z2=0.3365):
    z = np.linspace(z1, z2, 10000)
    Dl = model.angular_diameter_distance(z).to('kpc').value

    return interp1d(Dl, z, kind='cubic')


def Dls2Dl_convertor(zs=0.3365, z1=0, z2=0.3365):
    z = np.linspace(z1, z2, 10000)
    Dls = model.angular_diameter_distance_z1z2(z, zs).to('kpc').value
    Dl = model.angular_diameter_distance(z).to('kpc').value

    return interp1d(Dls, Dl, kind='cubic')

def Ds2Dls_convertor(zl, z2=0.3365):
    z = np.linspace(zl, z2, 10000)
    Dls = model.angular_diameter_distance_z1z2(zl, z).to('kpc').value
    Ds = model.angular_diameter_distance(z).to('kpc').value

    return interp1d(Ds, Dls, kind='cubic')

# a = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]], [[9, 10], [11, 12]]])
# print(np.concatenate([a[:, 0], a[:, 1]]))
