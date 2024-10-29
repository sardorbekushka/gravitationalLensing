from settings import *
from solver import *
from renderer import *

#
# solver = Solver(Lens(mass=constants.M_sun.value * 4e9, D_ls=6.8, center=np.array([0, 0])),
#                 Source(z=0.3365, center=[0, 0], direction=1e-2, angle=0, num=2000, source_type='line'))
solver = Solver(Lens(mass=constants.M_sun.value * 4e10, D_ls=6.8, center=np.array([0, 0])),
                Source(z=0.3365, center=[0, 0], direction=1e-1, angle=0, num=100, source_type='line'))

ax = plt.axes()
renderer = Renderer(solver, ax)
renderer.start()

# у величиин astropy есть размерности!!!

# можно источник задать сразу как (угол угол расстояние) (т.е. по углу angle задать направление источника а потом просто задать длину источника исходя из них и направления direction)