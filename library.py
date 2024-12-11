from solver import *
from renderer import *

D_s = (model.angular_diameter_distance(0.3365).to('kpc')).value
# M = constants.M_sun.value * 1e10
M = 1.591e40
D_ls = 15

einst2 = 4 * constants.G.value * M / constants.c.value ** 2 * D_ls / (D_s - D_ls) / D_s / 3.086e19 * (arcsec ** 2)
# einst = 0.75e-3
# angle = einst * einst2 / (einst2 - einst ** 2) * (D_s - D_ls) / D_ls / 3600
angle = -0.1
# angle = 45
# angle = -0.2
length1 = np.abs(einst2 ** 0.5 / arcsec * (D_s - D_ls) / (angle / 180 * np.pi) * 4)
length2 = 2e-2 / np.abs(np.sin(angle * rad))
length = min(length1, length2)
# angle = -90
# length = 2e-2
# if einst2 < einst ** 2:
#     print('WARNING. PARAMETERS SHOULD BE CORRECTED')

solver = Solver(Lens(mass=M, D_ls=D_ls, center=np.array([-2e-3, -8e-4])),
                Source(z=0.3365, center=[0, 0], direction=angle, angle=0,
                       num=5000, source_type='line', length=length, radius=3e-5))

ax = plt.axes()
renderer = Renderer(solver, ax)
renderer.gif()