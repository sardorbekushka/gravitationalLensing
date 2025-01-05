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

# angle = einst2 ** 0.5 * (D_s - D_ls) / D_ls / arcsec * 10
angle = 1.57
L = 2e-2 # kpc
length = L / (np.sin(np.radians(angle)) + L / D_s * np.cos(np.radians(angle)))

if einst2 < einst ** 2:
    print('WARNING. PARAMETERS SHOULD BE CORRECTED')

# solver = Solver(Lens(mass=M, D_ls=D_ls, center=np.array([0, -2e-3])),
#                 Source(z=0.3365, center=[0, 0], direction=-angle, angle=0,
#                        num=1000, source_type='line', length=length, radius=3e-3))

solver = Solver(Lens(mass=M, D_ls=D_ls, center=np.array([-1e-2, -2.4])),
                Source(z=0.3365, source_type='line', x0=0, x1=0, y0=0, y1=-4, d0=0, d1=7))
ax = plt.axes()
renderer = Renderer(solver, ax)
renderer.start()
