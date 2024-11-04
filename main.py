from settings import *
from solver import *
from renderer import *

#
# solver = Solver(Lens(mass=constants.M_sun.value * 4e9, D_ls=6.8, center=np.array([0, 0])),
#                 Source(z=0.3365, center=[0, 0], direction=1e-2, angle=0, num=2000, source_type='line'))
# solver = Solver(Lens(mass=constants.M_sun.value * 4e10, D_ls=7, center=np.array([0, 0])),
#                 Source(z=0.3365, center=[0, 0], direction=1, angle=0, num=1000, source_type='line', length=10))
#

D_s = (model.angular_diameter_distance(0.3365).to('kpc')).value
M = constants.M_sun.value * 1.5e10
D_ls = 10

einst2 = 4 * constants.G.value * M / constants.c.value ** 2 * D_ls / (D_s - D_ls) / D_s / 3.086e19 * (arcsec ** 2)
einst = 1e-3
angle = einst * einst2 / (einst2 - einst ** 2) * (D_s - D_ls) / D_ls / 3600

if einst2 < einst ** 2:
    print('WARNING. PARAMETERS SHOULD BE CORRECTED')

solver = Solver(Lens(mass=M, D_ls=D_ls, center=np.array([0, 0])),
                Source(z=0.3365, center=[0, 0], direction=-angle, angle=0,
                       num=10000, source_type='cylinder', length=D_ls * 1.2, radius=5e-5))

ax = plt.axes()
renderer = Renderer(solver, ax)
renderer.start()

# у величиин astropy есть размерности!!!
