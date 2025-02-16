import numpy as np

from settings import *
from solver import *
from renderer import *

D_s = (model.angular_diameter_distance(0.3365).to('kpc')).value
# M = constants.M_sun.value * 1.5e10
M = 1e40
D_ls = 7
einst2 = 4 * constants.G.value * M / constants.c.value ** 2 * D_ls / (D_s - D_ls) / D_s / 3.086e19 * (arcsec ** 2)
einst = 0.75e-3
print(einst2 ** 0.5)

# angle = einst2 ** 0.5 * (D_s - D_ls) / D_ls / arcsec * 10
angle = 1.57
L = 2e-2 # kpc
length = L / (np.sin(np.radians(angle)) + L / D_s * np.cos(np.radians(angle)))

if einst2 < einst ** 2:
    print('WARNING. PARAMETERS SHOULD BE CORRECTED')

# solver = Solver(Lens(mass=M, D_ls=D_ls, center=np.array([0, -2e-3])),
#                 Source(z=0.3365, center=[0, 0], direction=-angle, angle=0,
#                        num=1000, source_type='line', length=length, radius=3e-3))
#
ea = 1
k = 4 * constants.G.value * M / constants.c.value ** 2 / D_s / 3.086e19 * (arcsec ** 2)
# D_ls = D_s / (1 + k / ea * 1e3 ** 2)
# D_ls = 7
solver = Solver(Lens(mass=M, D_ls=D_ls, center=np.array([0.0, -0.2])),
                Source(z=0.3365, source_type='line', x0=0, x1=0.1, y0=0, y1=-4, d0=0, d1=20))
#


# solver = Solver(Lens(mass=M, D_ls=D_ls, center=np.array([0, 0])),
#                 Source(z=0.3365, source_type='circle', x0=0, x1=0.5, y0=0.05, y1=-0.05, d0=0, d1=2 * D_ls))
#

# y0 = -3.2
# l = 0
# y1 = y0 - 0.1
# d0 = y0 / -4 * 7
# d1 = y1 / -4 * 7
# solver = Solver(Lens(mass=M, D_ls=D_ls, center=np.array([-1e-2, -3.95])),
#                 Source(z=0.3365, source_type='line', x0=0, x1=0, y0=y0, y1=y1, d0=d0, d1=d1))

ax = plt.axes()
renderer = Renderer(solver, ax)
renderer.start()
