from settings import *
from solver import *
from renderer import *


solver = Solver(Lens(mass=constants.M_sun.value * 1e12, D_ls=5, center=np.array([0, 0])),
                Source(z=0.3365, center=[1e-3, 0], direction=0.5, angle=0, num=10000, source_type='line'))

ax = plt.axes()
renderer = Renderer(solver, ax)
renderer.show()

# у величиин astropy есть размерности!!!

# можно источник задать сразу как (угол угол расстояние) (т.е. по углу angle задать направление источника а потом просто задать длину источника исходя из них и направления direction)