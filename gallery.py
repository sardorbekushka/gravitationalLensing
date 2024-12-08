import os

import matplotlib.pyplot as plt
import numpy as np

from settings import *
from renderer import *
from solver import *

N = 3
m_ = np.linspace(-1, 1, N) ** 3
dls_ = np.linspace(-1, 1, N)
x_ = np.linspace(-0.9, 0.9, N) ** 3
# x_ = np.array([-0.5, -0.3, 0, 0.5])                             ###################################
# y_ = np.linspace(-4 ** (1/3), 2 ** (1/3), N) ** 3
y_ = np.linspace(-3, 1, N)
# y_ = np.array([0])                                              ###################################
a_ = 10 * np.linspace(0.4, 1, N) ** 4

M_range = constants.M_sun.value * 10 ** (2 * m_ + 10)        # kg
# M_range = [constants.M_sun.value * 1e10]
x_range = 2 * x_ * 1e-3                                        # as
y_range = y_ * 1e-3                                            # as
# X, Y = np.meshgrid(x_range, y_range)
Dls_range = 5 ** dls_                                      # в долях от среднего значения ###############################
# Dls_range = np.array([0.1, 0.5, 0.8, 1, 2, 5, 10])
# Dls_range = [15]
angle_range = 5 ** dls_                                            # угол в долях угла эйнштейна. перерасчет будет в коде плюс дополнить еще углами 45 и 90
# angle_range = [10]
D_s = (model.angular_diameter_distance(0.3365).to('kpc')).value
L = 2e-2  # kpc

def createGallery(self, directory, X, Y):
    self.showData = False
    cbar = plt.colorbar(self.show(), ax=self.ax)
    cbar.set_label('ln(magnification)', color=g)
    cbar.ax.tick_params(labelcolor=g)

    for y in Y:
        for x in X:
            self.solver.setLens(np.array([x, y]))
            self.show()
            plt.savefig(f'{directory}Y{round(y * 1e3, 3)}X{round(x * 1e3, 3)}.png', dpi=150)


Renderer.createGallery = createGallery

def createGallery_(directory, M, D_ls, A, X, Y):
    for m in M:
        dir = f'{"M{:.4e}".format(m)}kg/'
        directorym = directory + dir
        # os.mkdir(directorym)
        ea = 0.75
        D_ls_ = (ea / 1000 / arcsec * D_s * constants.c.value) ** 2 / 4 / constants.G.value / m * 3.086e19 * D_ls

        for d in D_ls_:
            einst = np.sqrt(4 * constants.G.value * m / constants.c.value ** 2 * d / (D_s - d) / D_s / 3.086e19) * arcsec
            dir = f'Dls{round(d, 2)}kpc_ea{round(einst * 1e3, 3)}mas/'
            directoryd = directorym + dir
            # if not (os.path.isdir(directoryd)):
            #     os.mkdir(directoryd)
            X = 2 * einst * np.linspace(-1, 1, N) ** 3
            for a in A:
                angle = einst * (D_s - d) / d / arcsec / rad * a
                # angle = 0.1
                dir = f'dir{round(angle, 5)}_{a}ea/'
                directorya = directoryd + dir
                os.makedirs(directorya)

                length = L / (np.sin(np.radians(angle)) + L / D_s * np.cos(np.radians(angle)))

                Y = np.linspace(-a * einst, einst / 2, N)

                solver = Solver(Lens(mass=m, D_ls=d, center=np.array([X[0], Y[0]])),
                                Source(z=0.3365, center=[0, 0], direction=-angle, angle=0,
                                       num=5000, source_type='cylinder', length=length, radius=3e-5))
                ax = plt.axes()
                renderer = Renderer(solver, ax)

                renderer.createGallery(directorya, X, Y)


dir = 'gallery6/'
# print(Dls_range)
# print(angle_range)
os.mkdir(dir)
createGallery_(dir, M_range, Dls_range, angle_range, x_range, y_range)
# print(x_)