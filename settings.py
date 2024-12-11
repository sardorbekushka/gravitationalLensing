import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as constants
import astropy.cosmology.funcs as cosmology
from astropy.cosmology import LambdaCDM
from matplotlib import rc
import os


# model = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
model = LambdaCDM(H0=67.8, Om0=0.308, Ode0=0.692)

g = np.array([0.8, 0.8, 0.9])
# rc('font', **{'family': 'Times new roman'})
# rc('font', **{'family': 'DejaVu Sans'})
rc('text', usetex=True)
# rc('text.latex',unicode=True)
rc('text.latex', preamble=r'\usepackage[utf8]{inputenc}')
rc('text.latex', preamble=r'\usepackage[russian]{babel}')
rc('axes', edgecolor=g)

lim = [[-2, 2], [-5, 2]] # [по горизонтали, по вертикали]. или None чтобы ограничение шло по углу Эйнштейна

# if lim
# if lim:
#     k = (lim[0][1] - lim[0][0]) / (lim[1][1] - lim[1][0])
#     figsize = (7 * k, 7)
# else:
#     figsize = (8.7, 7)
#
figsize = ((lim[0][1] - lim[0][0]) / (lim[1][1] - lim[1][0]) * 8 * 1.24, 8) if lim else (8.7, 7)
plt.figure(figsize=figsize, facecolor='black')
arcsec = 180 * 3600 / np.pi
rad = np.pi / 180
L = 2e-2 # kpc
plt.style.use('dark_background')
