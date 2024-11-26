import numpy as np

from settings import *


class Lens:
    def __init__(self, mass, D_ls, center) -> None:
        '''
        :param mass: the mass of the Lens (in kg)
        :param D_ls: distance between Lens and Source (in kpc)
        :param center: center of the Lens in seconds (in angle plane)
        '''
        self.D_ls = D_ls
        self.m = mass
        self.center = np.array(center, dtype=float)


class Source:
    def __init__(self, z=0.3365, center=[0, 0], direction=0, angle=0, num=1000, source_type='line', length=10, radius=3e-5) -> None:
        '''
        :param z: redshift of Source
        :param center: center of the Source in angle plane (in arcseconds)
        :param direction: angle between the jet direction and the observer view direction (in degrees)
        :param angle: rotation angle measured from the vertical in clockwise arrow (in degrees)
        '''

        self.z = z
        self.D_s = model.angular_diameter_distance(self.z).to('kpc').value
        self.center = np.array(center)
        self.direction = direction
        self.angle = angle
        self.length = length

        self.num = num
        self.points = self.createCylinderSource(length, radius, direction, angle, num) if source_type == 'cylinder' \
                 else self.createCircleSource(5e-2, 100) if source_type == 'circle' \
                 else self.createLineSource(length, direction, angle, num)

    def createLineSource(self, length, direction, angle, num):
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
        v = np.linspace([0, 0, 0], [0, 0, -length], num)
        v = self.rotateX(v, direction)
        v = self.rotateZ(v, angle)

        return self.scale(v)

    def createCircleSource(self, radius, num):
        '''
        creates the Source in shape of circle.

        output format [[x1, y1, d1], [x2, y2, d2], ... ], \
        where xi, yi are scaled into angles in arcseconds, di - angular diameter distance in kpc.

        :param radius: radius of the circle (in kpc)
        :param num: amount of points in the Source
        :return: an array of Source points
        '''
        theta = np.linspace(0, 2 * np.pi, num)
        radius = self.D_s * radius / arcsec
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)

        d = np.zeros(num)

        return self.scale(np.array([x, y, d]).T)

    def createCylinderSource(self, length, radius, direction, angle, num):
        '''
        creates the Source in shape of cylinder.

        output format [[x1, y1, d1], [x2, y2, d2], ... ], \
        where xi, yi are scaled into angles in arcseconds, di - angular diameter distance in kpc.

        :param length: length of the source i.e. cylinder height (in kpc)
        :param radius: radius of the source (in kpc)
        :param direction: angle between the jet direction and the observer view direction (in degrees)
        :param angle: rotation angle measured from the vertical in clockwise arrow (in degrees)
        :param num:  amount of points in the Source
        :return: an array of Source points
        '''

        # z = np.linspace(0, -length, round(num ** (1/3) * 10))
        # theta = np.linspace(0, 2 * np.pi, round(num ** (1/3) / 20))
        # radii = np.linspace(0, radius, round(num ** (1/3) * 2))
        z = np.linspace(-length, 0, num // 40)
        theta = np.linspace(0, 2 * np.pi, 4)
        radii = np.linspace(0, radius, 10)

        R, A = np.meshgrid(radii, theta)
        x = R * np.cos(A)
        y = R * np.sin(A)

        X, Z = np.meshgrid(x, z)
        Y, _ = np.meshgrid(y, z)

        points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T

        return self.scale(self.rotateZ(self.rotateX(points, direction), angle))

    def rotateX(self, points, angle):
        cos = np.cos(np.radians(angle))
        sin = np.sin(np.radians(angle))
        A = np.array([[1, 0, 0], [0, cos, -sin], [0, sin, cos]])

        return list(map(lambda p: A @ p, points))

    def rotateZ(self, points, angle):
        cos = np.cos(np.radians(angle))
        sin = np.sin(np.radians(angle))
        A = np.array([[cos, sin, 0], [-sin, cos, 0], [0, 0, 1]])

        return list(map(lambda p: A @ p, points))

    def rotateY(self, points, angle):
        cos = np.cos(np.radians(angle))
        sin = np.sin(np.radians(angle))
        A = np.array([[cos, 0, sin], [0, 1, 0], [-sin, 0, cos]])

        return list(map(lambda p: A @ p, points))

    def scale(self, points):
        '''
        scales the coordinates from kpc into angles and moves the center of the Source to the given point
        :param points: an array of points
        :return:
        '''
        s = np.array(points).T
        s[0:2] *= arcsec / self.D_s
        s = s.T
        s += np.array([self.center[0], self.center[1], self.D_s])

        return s

    def unscale(self, points):
        s = np.array(points)
        s -= np.array([self.center[0], self.center[1], self.D_s])
        s = s.T
        s[0:2] /= arcsec / self.D_s

        return s.T

    def update(self, direction, angle):
        s = self.unscale(self.points)
        self.rotateX(self.rotateZ(s, -self.angle), -self.direction)
        self.direction = direction
        self.angle = angle
        self.scale(self.rotateZ(self.rotateX(s, self.direction), self.angle))
        self.points = s

    def setDirection(self, d):
        self.direction = d
        self.setPoints(self.createCylinderSource(self.length, 3e-5, d, self.angle, self.num))

    def setPoints(self, points):
        self.points = points


class Solver:
    def __init__(self, lens: Lens, source: Source) -> None:
        '''
        :param lens: an object of the Lens class. the lens of the system
        :param source: an object of the Source class. the source of the system
        '''
        self.lens = lens
        self.source = source
        self.stepMass = 2
        self.stepPos = 1e-5
        self.stepRot = 0.1
        self.D_l = self.source.D_s - self.lens.D_ls

    def einsteinRadius(self, D_s):
        """
        the Einstein radius for point with angular diameter distance D_s.
        :param D_s: an angular diameter distance of point
        :return: the Einstein radius in arcseconds
        """
        D_ls = D_s - self.D_l

        return 0 if D_ls < 0 else np.sqrt(4 * constants.G.value * self.lens.m / constants.c.value ** 2 * D_ls / self.D_l / D_s / 3.086e19) * arcsec

    def einsteinRadiusZ(self, z):
        """
        the Einstein radius for point with redshift equal to z.
        :param z: redshift of point
        :return: the Einstein radius in arcseconds
        """
        D_s = model.angular_diameter_distance(z)

        return self.einsteinRadius(D_s)

    def processPoint(self, p_):
        '''
        processes the images of a single point.
        :param p_: the point to be lensed. format [x1, y1, d1]
        :return: the coordinates of two images in format [im1, im2]
        '''
        p = np.matrix(p_[0:2])
        D_s = p_[2]

        if D_s <= self.D_l:
            pp = np.array(p)[0]

            return np.array([pp, pp])

        dp = p - self.lens.center
        beta = np.linalg.norm(dp)
        beta2 = beta ** 2
        einstAngle = self.einsteinRadius(D_s)
        v = np.matrix([1, -1])

        angle = (beta + v * np.sqrt(beta2 + 4 * einstAngle ** 2)) / 2
        points = angle.T * dp / beta + self.lens.center

        return np.array(points)


    def magnification(self, p, ea):
        m = 1 / (1 - (ea / np.linalg.norm(p - self.lens.center)) ** 4)
        return 1 if ea == 0 else np.abs(m)

    def processImage(self):
        '''
        processes a whole image of the source
        :return: image coordinates of each point in source. format [[im1_x, im1_y], [im2_x, im2_y], ... ]
        '''
        points_ = list(map(lambda p: self.processPoint(p), self.source.points)) # здесь формат [[[p1_1.x, p1_1.y], [p1_2.x, p1_2.y]],   [[p2_1.x, ..]]]
        points = [pp for p in points_ for pp in p]     # здесь формат [[p1_1.x, p1_1.y], [p1_2.x, p1_2.y], [p2_1.x, ..] ..]

        eas = list(map(lambda Ds: self.einsteinRadius(Ds), self.source.points.T[2]))
        eas2 = [val for pair in zip(eas, eas) for val in pair]

        m = list(map(lambda p, ea: self.magnification(p, ea), points, eas2))

        return np.transpose(points), m

    def moveLens(self, shift):
        self.lens.center += np.array(shift) * self.stepPos

    def setLens(self, pos):
        self.lens.center = pos

    def increaseMass(self):
        self.lens.m *= self.stepMass

    def decreaseMass(self):
        self.lens.m /= self.stepMass

    def setMass(self, M):
        self.lens.m = M

    def turnSource(self, ang, direct):
        self.source.update(self.source.direction + direct * self.stepRot,
                           self.source.angle + ang * self.stepRot)

    def setDirection(self, d):
        self.source.setDirection(d)

    def declineSource(self, k):
        self.setDirection(self.source.direction + k * self.stepRot)

    def getEinsteinRadius(self):
        return self.einsteinRadius(self.source.D_s)

    def getLensMass(self):
        return self.lens.m

    def getSourceRedshift(self):
        return self.source.z

    def getSourceDistance(self):
        return self.source.D_s

    def getLensDistance(self):
        return self.source.D_s - self.lens.D_ls

    def getSourceDirection(self):
        return self.source.direction

    def getLensCenter(self):
        return self.lens.center

    def getSourceCenter(self):
        return self.source.center

