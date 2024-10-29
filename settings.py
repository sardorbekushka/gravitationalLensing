import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as constants
import astropy.cosmology.funcs as cosmology
from astropy.cosmology import LambdaCDM
from matplotlib import rc

model = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

g = np.array([0.8, 0.8, 0.9])
# rc('font', **{'family': 'Times new roman'})
# rc('font', **{'family': 'DejaVu Sans'})
rc('text', usetex=True)
# rc('text.latex',unicode=True)
rc('text.latex', preamble=r'\usepackage[utf8]{inputenc}')
rc('text.latex', preamble=r'\usepackage[russian]{babel}')
rc('axes', edgecolor=g)

plt.figure(figsize=(8.7, 7), facecolor='black')
arcsec = 180 * 3600 / np.pi