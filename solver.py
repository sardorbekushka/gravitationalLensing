import numpy as np

from settings import *
from utils import *
from scipy.interpolate import interp1d

rs_coeff = 2 * constants.G.value / constants.c.value ** 2
mas2m = deg2rad * sec2deg / 1e3
kpc2m = 3.086e19
CONVERSION_FACTOR = 3.086e19 * (deg2rad * sec2deg / 1e3) ** 2

class Lens:
    def __init__(self, mass=1e40, D_ls=7,
                 center=np.array([0.0, 0.0])) -> None:
        """
        :param mass: the mass of the Lens (in kg)
        :param D_ls: distance between Lens and Source (in kpc)
        :param center: center of the Lens in seconds (in angle plane)
        """
        self.D_ls = D_ls
        self.m = mass
        self.center = np.array(center, dtype=float)
        self.D_l = None
        self.z = None

    def set_Dl(self, Dl):
        self.D_l = Dl

    def set_z(self, z):
        self.z = z

    @property
    def get_schwarzschild_radius(self):
        return rs_coeff * self.m

class Source:
    def __init__(self, z=0.3365, source_type='line',
                 x0=0, x1=0, y0=0, y1=-4, d0=0, d1=7) -> None:
        """
        :param z: redshift of Source
        """
        self.z = z
        self.D_s = model.angular_diameter_distance(self.z).to('kpc').value
        self._setup_geometry(source_type, x0, x1, y0, y1, d0, d1)

    def _setup_geometry(self, source_type, x0, x1, y0, y1, d0, d1):
        self.x0, self.x1 = x0, x1
        self.y0, self.y1 = y0, y1
        self.d0, self.d1 = d0, d1
        self.source_type = source_type
        geometry_creators = {
            'line': self._create_line_geometry
        }

        self.x, self.d, self.direction = geometry_creators[source_type](x0, x1, y0, y1, d0, d1)

    def _create_line_geometry(self, x0, x1, y0, y1, d0, d1):
        x = lambda y: (y - y0) / (y1 - y0) * (x1 - x0) + x0
        d = lambda y: ((self.D_s * (y - y0) * (d1 - d0) + d0 * (y1 * (self.D_s - d1) - y0 * self.D_s)) /
                       (y * (d1 - d0) + y1 * (self.D_s - d1) - y0 * self.D_s))
        angle = np.arctan(np.sqrt((y1 - y0) ** 2 + (x1 - x0) ** 2) * mas2m * self.D_s / (d1 - d0))

        return x, d, np.rad2deg(angle)

    def update_geometry(self, **params):
        self._setup_geometry(self.source_type, **params)


class Solver:
    def __init__(self, lens: Lens, source: Source) -> None:
        """
        :param lens: an object of the Lens class. the lens of the system
        :param source: an object of the Source class. the source of the system
        """

        self.lens = lens
        self.source = source
        self._setup_constants()
        self._setup_convertors()
        self._update_k()
        self.points = np.array([])
        self.images = np.array([])
        self.magnifications = np.array([])

    def _setup_constants(self):
        self.step_mass = 1.2
        self.step_pos = 1e-1
        self.step_length = 1

    def _update_k(self):
        self.k = 2 * rs_coeff * self.lens.m / self.lens.D_l / CONVERSION_FACTOR

    def _setup_convertors(self):
        self.Dls2Dl = Dls2Dl_convertor(self.source.z)
        self.Dl2z = D2z_convertor()
        self._setup_lens(self.lens.D_ls)

    def _setup_lens(self, Dls):
        # print(Dls)
        self.lens.D_ls = Dls
        self.lens.set_Dl(self.Dls2Dl(Dls))
        self.lens.set_z(self.Dl2z(self.lens.D_l))
        self.Ds2Dls = Ds2Dls_convertor(self.lens.z)

    def einstein_radius(self, D_s):
        """
        the Einstein radius for point with angular diameter distance D_s.
        :param D_s: an angular diameter distance of point
        :return: the Einstein radius in arcseconds
        """
        # print(D_s, self.Ds2Dls(D_s), self.lens.D_ls)

        return 0 if D_s < self.lens.D_l else np.sqrt(self.k * self.Ds2Dls(D_s) / D_s)

    def einstein_radius_z(self, z):
        """
        the Einstein radius for point with redshift equal to z.
        :param z: redshift of point
        :return: the Einstein radius in arcseconds
        """
        D_s = model.angular_diameter_distance(z)

        return self.einstein_radius(D_s)

    def process_point(self, y):
        x = self.source.x(y)
        d = self.source.d(y)
        p = np.array([x, y])
        dp = p - self.lens.center

        beta = np.linalg.norm(dp)
        beta2 = beta ** 2
        ea2 = self.einstein_radius(self.Dls2Dl(d)) ** 2

        theta1 = (beta + np.sqrt(beta2 + 4 * ea2)) / 2
        theta2 = (beta - np.sqrt(beta2 + 4 * ea2)) / 2

        im1_ = theta1 / beta * dp + self.lens.center
        im2_ = theta2 / beta * dp + self.lens.center

        m1_ = 1 / (1 - (ea2 / theta1 ** 2) ** 2)
        m2_ = -1 / (1 - (ea2 / theta2 ** 2) ** 2)

        return p, np.array([im1_, im2_]), np.array([m1_, m2_])

    def process_image(self, d_min=5e-3, d_max=1e-2):
        y_values, dy = [self.source.y0], 1e-2
        images, magnifications, points = [], [], []
        # magnifications = []
        # points = []

        y_new = y_values[-1] - dy
        y_values.append(y_new)
        p, im, m = self.process_point(y_new)
        points.append(p)
        images.append(im)
        magnifications.append(m)

        while y_values[-1] > self.source.y1:
            y_new = y_values[-1] - dy
            y_values.append(y_new)
            p, im, m = self.process_point(y_new)
            points.append(p)
            images.append(im)
            magnifications.append(m)

            if np.any(np.linalg.norm(images[-1] - images[-2]) > d_max):
                dy /= 2
            elif np.all(np.linalg.norm(images[-1] - images[-2]) < d_min):
                dy *= 2

        self.points = np.array(points)
        images = np.array(images)
        self.images = np.concatenate([images[:, 0], images[:, 1]])
        self.magnifications = np.array(magnifications).T.flatten()

        return self.images, self.magnifications

    def update_source(self, x1, d1):
        self.source.update_geometry(x0=self.source.x0, x1=x1,
                                    y0=self.source.y0, y1=self.source.y1,
                                    d0=self.source.d0, d1=d1)

    def move_lens(self, shift):
        self.set_lens_center(self.lens.center + np.array(shift) * self.step_pos)

    def set_lens_center(self, pos):
        self.lens.center = np.array(pos)

    def set_Dls(self, Dls):
        self._setup_lens(Dls)

    def move_Dls(self, k):
        self.set_Dls(self.lens.D_ls + k * 10 * self.step_length)

    def set_mass(self, m):
        self.lens.m = m
        self._update_k()

    def change_mass(self, k):
        self.set_mass(self.lens.m * self.step_mass ** k)

    def set_length(self, d1):
        self.update_source(self.source.x1, d1)

    def decline_source(self, k):
        self.set_length(self.source.d1 * 1.1 ** k)

    @property
    def get_einstein_radius(self):
        return self.einstein_radius(self.source.D_s)

    @property
    def get_lens_mass(self):
        return self.lens.m

    @property
    def get_source_redshift(self):
        return self.source.z

    @property
    def get_source_distance(self):
        return self.source.D_s

    @property
    def get_lens_distance(self):
        return self.lens.D_l

    @property
    def get_source_direction(self):
        return self.source.direction

    @property
    def get_lens_center(self):
        return self.lens.center

    @property
    def get_source_points(self):
        return self.points


