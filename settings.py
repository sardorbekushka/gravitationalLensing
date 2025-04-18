import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as constants
import astropy.cosmology.funcs as cosmology
from astropy.cosmology import LambdaCDM
from matplotlib import rc
import matplotlib.colors as colors
import os

# model = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
model = LambdaCDM(H0=67.8, Om0=0.308, Ode0=0.692)

g = np.array([0.8, 0.8, 0.9])
# rc('font', **{'family': 'Times new roman'})
# rc('font', **{'family': 'DejaVu Sans'})
 # [по горизонтали, по вертикали]. или None чтобы ограничение шло по углу Эйнштейна
# lim = [[-1, 1], [-2.5, 1]] # [по горизонтали, по вертикали]. или None чтобы ограничение шло по углу Эйнштейна

# if lim
# if lim:
#     k = (lim[0][1] - lim[0][0]) / (lim[1][1] - lim[1][0])
#     figsize = (7 * k, 7)
# else:
#     figsize = (8.7, 7)
#
sec2deg = 1 / 3600
deg2rad = np.pi / 180
asec = 0.000012120342027738

arcsec = 180 * 3600 / np.pi
rad = np.pi / 180
L = 2e-2 # kpc
# plt.rcParams.update({'font.size': 14})


kpc2m = 3.086e19
G = constants.G.value
c2 = constants.c.value ** 2
CONVERSION_FACTOR = 3.086e19 * (deg2rad * sec2deg / 1e3) ** 2
flux = 0.01448780487804878

# data_old15kHz.txt = readData()
# mask = data_old15kHz.txt[3] < -3
# # print(data_old15kHz.txt[3] < -3)
# print(data_old15kHz.txt.T[mask])

## в среднем в районе y = -4 потол flux=0.044