import numpy as np

from settings import *
from solver import *
from renderer import *

M = 1e40
D_ls = 7

solver = Solver(Lens(mass=M, D_ls=D_ls, center=np.array([0.0, 0.0])),
                Source(z=0.3365, source_type='line', x0=0, x1=0.1, y0=0, y1=-4, d0=0, d1=20))

ax = plt.axes()
renderer = Renderer(solver, ax)
renderer.start()
