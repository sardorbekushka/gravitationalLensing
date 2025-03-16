import numpy as np

from settings import *
from solver import *
from renderer import *

M = 2e40
D_ls = 15
D_s = model.angular_diameter_distance(0.3365).to('kpc').value
D_l = D_s - D_ls
ea2 = k = 4 * constants.G.value * M / constants.c.value ** 2 / D_l / 3.086e19 / (
            deg2rad * sec2deg / 1e3) ** 2 * D_ls / D_s
x0 = 0.1
x1 = (x0 + (x0 ** 2 + ea2) ** 0.5) / 2

M = 3.456e40
m = 5e9 * 2e30
y = -1.2
x = -7.5e-2
# x = -5.7456e-02
# y = -1.0086
d =  4.96e-3 * 4 / np.deg2rad(0.0284)
d1 = 41.9464
print(np.deg2rad(0.0284))
D_ls = 15
d_ls = D_ls / m * M
# d = d * d_ls / D_ls
y = -0.94898
x = 1.40816e-04
x1 = 0.05
M = 3.11040e+40
# 0.3469_0.5909.png
#
# M = 3.11040000e+40
# D_ls = 1.67283542e+01
# x = 2.77434895e-02
# y = -1.09520349e+00
# x1 = 1.00837621e-01
# d1 = 4.19240139e+01
# M = 3.11040000e+40
# D_ls = 2.95912305e+01
# x = 8.97892366e-02
# y = -7.74792448e-01
# x1 = 7.00000000e-01
# d1 = 1.12490905e+02
M = 1e40
# D_ls = 20
# D_ls = 27.44190116
# x = 0.11003923
# y = -0.92784841
# x1 = 0.7
# d1 = 35.80732602
#
# D_ls = 14.51821294
# x = 0.03396592
# y = -1.12152954
# x1 = 0.08797423
# d1 = -4.30708758

# dev = 80.5
M = 1e40
D_ls = 2.71441354e+01
x = 2.74456158e-02
y = -1.10237773e+00
x1 = 9.30838467e-02
d1 = 3.95408677e+01

# dev =
M = 1e39

solver = Solver(Lens(mass=M, D_ls=D_ls, center=np.array([x, y])),
                Source(z=0.3365, source_type='line', x0=0.0, x1=x1, y0=0.0, y1=-4, d0=0, d1=d1),
                filename='src/data_15kHz.txt')

ax = plt.axes()
renderer = Renderer(solver, ax)


renderer.start()