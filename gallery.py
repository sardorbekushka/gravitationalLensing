from renderer import *
from solver import *

# d1_range in parts of D_ls
def createGallery(self, directory, m, ea_range, d1_range, num):
    self.showData = False
    self.showSource = False

    k = 4 * constants.G.value * m / constants.c.value ** 2 / self.solver.source.D_s / 3.086e19 / (
                deg2rad * sec2deg) ** 2

    for ea in ea_range:
        D_ls = self.solver.source.D_s  / (1 + k / (ea * 1e-3) ** 2)

        self.solver.setMass(m)
        self.solver.setDls(D_ls)

        dir = f'{"ea_{:.2e}".format(self.solver.getEinsteinRadius())}mas_{"Dls_{:.2e}".format(self.solver.lens.D_ls)}kpc/'
        directoryea = directory + dir

        for d in d1_range:
            self.solver.setLength(d * D_ls)

            dir = f'{"angle_{:.2e}".format(self.solver.getSourceDirection())}/'
            directoryd = directoryea + dir
            os.makedirs(directoryd)

            x_range = np.linspace(-1, 0, num) ** 3 * ea
            y_range = np.linspace(min(D_ls / self.solver.source.d1, 1) * self.solver.source.y1, 0.5, num)
            for y in y_range:
                for x in x_range:
                    self.solver.setLens(np.array([x+1e-3, y]))
                    self.show()
                    plt.savefig(f'{directoryd}y{round(y, 2)}_x{round(x, 2)}.png', dpi=150)

Renderer.createGallery = createGallery


M = 1e40
# ea_range = [0.75]
ea_range = [0.5, 0.75, 1]
d1_range = [10, 2, 1, 0.5, 0]
num = 6

solver = Solver(Lens(), Source())
ax = plt.axes()
renderer = Renderer(solver, ax)

dir = f'{"M_{:.2e}".format(M)}kg/'

# os.mkdir(dir)
renderer.createGallery(dir, M, ea_range, d1_range, num)



