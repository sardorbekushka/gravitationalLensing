import numpy as np

from settings import *

class Lens:
    def __init__(self, mass=1e40, D_ls=7, center=np.array([0.0, 0.0])) -> None:
        '''
        :param mass: the mass of the Lens (in kg)
        :param D_ls: distance between Lens and Source (in kpc)
        :param center: center of the Lens in seconds (in angle plane)
        '''
        self.D_ls = D_ls
        self.m = mass
        self.center = np.array(center, dtype=float)

    def getSchwarzschildRadius(self):
        return 2 * constants.G.value * self.m / constants.c.value ** 2

class Source:
    def __init__(self, z=0.3365, source_type='line', x0=0, x1=0, y0=0, y1=-4, d0=0, d1=7, width=0) -> None:
        '''
        :param z: redshift of Source
        :param center: center of the Source in angle plane (in arcseconds)
        :param direction: angle between the jet direction and the observer view direction (in degrees)
        :param angle: rotation angle measured from the vertical in clockwise arrow (in degrees)
        '''

        self.z = z
        self.D_s = model.angular_diameter_distance(self.z).to('kpc').value

        self.points = np.array([[0, 0]])

        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.d0 = d0
        self.d1 = d1
        self.width = width

        self.x, self.d, self.direction = self.createLineSource(x0, x1, y0, y1, d0, d1, width) if source_type == 'line' else (
                                         self.createCircleSource(x0, x1, y0, y1, d0, d1) if source_type == 'circle' else
                                         self.createParabolicSource(x0, x1, y0, y1, d0, d1) if source_type == 'parabolic' else None)

    def createLineSource(self, x0, x1, y0, y1, d0, d1, width):
        '''
        creates the Source in shape of line .

        output format [[x1, y1, d1], [x2, y2, d2], ... ], \
        where xi, yi are scaled into angles in arcseconds, di - angular diameter distance in kpc.

        :param length: length of the source i.e. cylinder height (in kpc)
        :param direction: angle between the jet direction and the observer view direction (in degrees)
        :param angle: rotation angle measured from the vertical in clockwise arrow (in degrees)
        :param num: amount of points in the Source
        :return: an array of Source points
        '''

        x = lambda y: (y - y0) / (y1 - y0) * (x1 - x0) + x0 + width * ((int(10000 * y) % 15) - 7) / 500
        # d = lambda y: (y - y0) / (y1 - y0) * (d1 - d0) + d0
        d = lambda y: ((self.D_s * (y - y0) * (d1 - d0) + d0 * (y1 * (self.D_s - d1) - y0 * self.D_s)) /
                       (y * (d1 - d0) + y1 * (self.D_s - d1) - y0 * self.D_s))
        angle = np.arctan(np.sqrt((y1 - y0) ** 2 + (x1 - x0) ** 2) * sec2deg * deg2rad / 1000 * self.D_s / (d1 - d0))

        return x, d, np.rad2deg(angle)

    def createParabolicSource(self, x0, x1, y0, y1, d0, d1):
        x = lambda y: (y - y0) * (y - y1) * x1 * (2 / (y1 - y0)) ** 2 + x0
        d = lambda y: ((self.D_s * (y - y0) * (d1 - d0) + d0 * (y1 * (self.D_s - d1) - y0 * self.D_s)) /
                       (y * (d1 - d0) + y1 * (self.D_s - d1) - y0 * self.D_s))

        angle = np.arctan(np.sqrt((y1 - y0) ** 2 + (x1 - x0) ** 2) * sec2deg * deg2rad / 1000 * self.D_s / (d1 - d0))

        return x, d, np.rad2deg(angle)

    def createCircleSource(self, x0, x1, y0, y1, d0, d1):
        x = lambda y: (x1 - x0) * np.sqrt(1 - (y - (y1 + y0) / 2) ** 2 / ((y1 - y0) / 2) ** 2) + x0
        d = lambda y: ((self.D_s * (y - y0) * (d1 - d0) + d0 * (y1 * (self.D_s - d1) - y0 * self.D_s)) /
                       (y * (d1 - d0) + y1 * (self.D_s - d1) - y0 * self.D_s))
        angle = np.arctan((y1 - y0) * sec2deg * deg2rad / 1000 * self.D_s / (d1 - d0))

        return x, d, np.rad2deg(angle)

    def updateLineSource(self, x0, x1, y0, y1, d0, d1, width):
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.d0 = d0
        self.d1 = d1
        self.width = width

        self.x, self.d, self.direction = self.createLineSource(x0, x1, y0, y1, d0, d1, width)

        return self.direction

    def setPoints(self, p):
        self.points = np.array(p)

class Solver:
    def __init__(self, lens: Lens, source: Source, filename='src/data_15kHz.txt') -> None:
        """
        :param lens: an object of the Lens class. the lens of the system
        :param source: an object of the Source class. the source of the system
        """
        self.lens = lens
        self.source = source
        self.stepMass = 1.2
        self.stepPos = 5e-3
        self.stepLength = 1
        self.D_l = self.source.D_s - self.lens.D_ls
        self.k = 4 * constants.G.value * self.lens.m / constants.c.value ** 2 / self.D_l / 3.086e19 / (deg2rad * sec2deg / 1e3) ** 2
        self.points = []
        self.image_points = np.array([])

        self.real_data = readData(filename)
        self.countMatch = None
        self.createMatch()

    def createMatch(self):
        t, x, dx, y, dy, l, dl = self.real_data
        x1 = (x - dx)[:, None]
        x2 = (x + dx)[:, None]
        y1 = (y - dy)[:, None]
        y2 = (y + dy)[:, None]
        xy = np.array([x, y]).transpose()
        ### поиск по точкам модели
        # self.countMatch = lambda p: np.any(((x1 <= p[0]) & (p[0] <= x2) & (y1 <= p[1]) & (p[1] <= y2)), axis=0).sum()
        # self.countMatch = lambda p: self.image_points.transpose()[np.any(((x1 <= p[0]) & (p[0] <= x2) & (y1 <= p[1]) & (p[1] <= y2)), axis=0)]

        ### поиск по реальным точкам
        self.countMatch = lambda p: np.any(((x1 <= p[0]) & (p[0] <= x2) & (y1 <= p[1]) & (p[1] <= y2)), axis=1).sum()
        # self.countMatch = lambda p: xy[np.any(((x1 <= p[0]) & (p[0] <= x2) & (y1 <= p[1]) & (p[1] <= y2)), axis=1)]

    def einsteinRadius(self, D_s):
        """
        the Einstein radius for point with angular diameter distance D_s.
        :param D_s: an angular diameter distance of point
        :return: the Einstein radius in arcseconds
        """
        D_ls = D_s - self.D_l

        return 0 if D_ls < 0 else np.sqrt(self.k * D_ls / D_s)

    def einsteinRadiusZ(self, z):
        """
        the Einstein radius for point with redshift equal to z.
        :param z: redshift of point
        :return: the Einstein radius in arcseconds
        """
        D_s = model.angular_diameter_distance(z)

        return self.einsteinRadius(D_s)

    def processPoint(self, y):
        x = self.source.x(y)
        d = self.source.d(y)
        self.points.append([x, y])
        p = np.array([x, y])
        dp = p - self.lens.center

        beta = np.linalg.norm(dp)
        beta2 = beta ** 2
        ea2 = self.einsteinRadius(self.source.D_s - d) ** 2

        theta1 = (beta + np.sqrt(beta2 + 4 * ea2)) / 2
        theta2 = (beta - np.sqrt(beta2 + 4 * ea2)) / 2

        im1_ = theta1 / beta * dp + self.lens.center
        im2_ = theta2 / beta * dp + self.lens.center

        m1_ = 1 / (1 - (ea2 / theta1 ** 2) ** 2)
        m2_ = -1 / (1 - (ea2 / theta2 ** 2) ** 2)

        return im1_, im2_, m1_, m2_

    def processImage(self, dy_min=1e-6, dy_max=1e-2):
        self.points = []

        y = self.source.y0
        dy = dy_max

        image1 = []
        magn1 = []
        image2 = []
        magn2 = []

        im1, im2, m1, m2 = self.processPoint(y)
        image1.append(im1)
        image2.append(im2)
        magn1.append(m1)
        magn2.append(m2)
        i = 0

        dmin = 5e-3
        dmax = 1e-2
        # dmin = 1e-3
        # dmax = 1e-2

        while y > self.source.y1:
            i += 1
            y -= dy
            im1, im2, m1, m2 = self.processPoint(y)
            image1.append(im1)
            image2.append(im2)
            magn1.append(m1)
            magn2.append(m2)

            if np.linalg.norm(im1 - image1[i - 1]) > dmax or np.linalg.norm(im2 - image2[i - 1]) > dmax:
                if dy > dy_min:
                    dy /= 2
            elif np.linalg.norm(im1 - image1[i - 1]) < dmin or np.linalg.norm(im2 - image2[i - 1]) < dmin:
                if dy < dy_max:
                    dy *= 2

        im = np.concatenate([image1, image2])
        im = np.transpose(im)
        magn = np.concatenate([magn1, magn2])
        self.source.setPoints(self.points)
        self.image_points = np.array(im)

        # print(len(image1))
        return im, magn

    def getEfficiency(self):
        return self.countMatch(self.image_points) / len(self.real_data[0])

    def changeWidth(self, k):
        self.source.updateLineSource(self.source.x0, self.source.x1, self.source.y0, self.source.y1, self.source.d0, self.source.d1, self.source.width + 0.5 * k)

    def moveLens(self, shift):
        self.lens.center += np.array(shift) * self.stepPos

    def setLens(self, pos):
        self.lens.center = pos

    def setDls(self, Dls):
        self.lens.D_ls = Dls
        self.D_l = self.source.D_s - self.lens.D_ls

    def moveDls(self, k):
        self.setDls(self.lens.D_ls + k * self.stepLength)

    def changeMass(self, k):
        self.lens.m *= self.stepMass ** k
        self.k = 4 * constants.G.value * self.lens.m / constants.c.value ** 2 / self.source.D_s / 3.086e19 / (deg2rad * sec2deg / 1e3) ** 2

    def decreaseMass(self):
        self.lens.m /= self.stepMass
        self.k = 4 * constants.G.value * self.lens.m / constants.c.value ** 2 / self.source.D_s / 3.086e19 / (deg2rad * sec2deg / 1e3) ** 2

    def setMass(self, M):
        self.lens.m = M
        self.k = 4 * constants.G.value * self.lens.m / constants.c.value ** 2 / self.source.D_s / 3.086e19 / (deg2rad * sec2deg / 1e3) ** 2

    def setLength(self, d1):
        self.source.updateLineSource(self.source.x0, self.source.x1, self.source.y0, self.source.y1, self.source.d0, d1, self.source.width)

    def declineSource(self, k):
        self.setLength(self.source.d1 + k * self.stepLength)

    def getEinsteinRadius(self):
        return self.einsteinRadius(self.source.D_s)

    def getLensMass(self):
        return self.lens.m

    def getSourceRedshift(self):
        return self.source.z

    def getSourceDistance(self):
        return self.source.D_s

    def getLensDistance(self):
        return self.D_l

    def getSourceDirection(self):
        return self.source.direction

    def getLensCenter(self):
        return self.lens.center

    def getSourcePoints(self):
        return np.array(self.points)


