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
        self.center = center


class Source:
    def __init__(self, z=0.3365, center=[0, 0], direction=90, angle=0, num=1000, source_type='line',) -> None:
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
        length = self.D_s / np.radians(max(direction, 0.5)) * 5e-8

        self.points = self.createCylinderSource(length, 3e-3, direction, angle, num) if source_type == 'cylinder' \
                 else self.createCircleSource(1e-2, 100) if source_type == 'circle' \
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
        # z = np.linspace(0, length, round(num ** (1 / 3) / 9))
        # theta = np.linspace(0, 2 * np.pi, round(num ** (1 / 3)) * 9)
        # radii = np.linspace(0, radius, round(num ** (1 / 3)))

        z = np.linspace(0, -length, round(num ** (1/3) * max(direction, 1) + 1))
        theta = np.linspace(0, 2 * np.pi, 3 * round(num ** (1/3) + 1))
        radii = np.linspace(0, radius, round(num ** (1/3) / max(direction, 1) / 3 + 1))

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


class Solver:
    def __init__(self, lens: Lens, source: Source) -> None:
        '''
        :param lens: an object of the Lens class. the lens of the system
        :param source: an object of the Source class. the source of the system
        '''
        self.lens = lens
        self.source = source

    def einsteinRadius(self, D_s):
        '''
        the Einstein radius for point with angular diameter distance equal to z.
        :param D_s: an angular diameter distance of point
        :return: the Einstein radius in degrees
        '''
        D_ls = self.lens.D_ls
        D_l = D_s - D_ls

        return np.sqrt(4 * constants.G.value * self.lens.m / constants.c.value ** 2 * D_ls / D_l / D_s / 3e19) * arcsec

    def einsteinRadiusZ(self, z):
        '''
        the Einstein radius for point with redshift equal to z.
        :param z: redshift of point
        :return: the Einstein radius in degrees
        '''
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
        if self.source.D_s - D_s > self.lens.D_ls:
            pp = np.array(p)[0]
            # print('FFF', pp)
            return np.array([pp, pp])

        dp = p - self.lens.center
        beta = np.linalg.norm(dp)
        beta2 = beta ** 2
        einstAngle = self.einsteinRadius(D_s)
        v = np.matrix([1, -1])

        angle = beta + v * np.sqrt(beta2 + 4 * einstAngle ** 2) / 2
        points = angle.T * dp / beta + self.lens.center

        return np.array(points)

    def processImage(self):
        '''
        processes a whole image of the source
        :return: image coordinates of each point in source. format [[im1_1, im1_2], [im2_1, im2_2], ... ]
        '''
        points = list(map(lambda p: self.processPoint(p), self.source.points))

        return np.transpose([pp for p in points for pp in p])

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

