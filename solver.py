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
    def __init__(self, z=0.3365, center=[0, 0], direction=90, angle=0, length=0.1, num=1000) -> None:
        '''
        :param z: redshift of Source
        :param center: center of the Source in angle plane (in arcseconds)
        :param direction: the angle of the Source radiation (in degrees)
        :param angle: the rotation angle of the Source (in degrees)
        '''

        self.z = z
        self.D_s = model.angular_diameter_distance(self.z).to('kpc').value
        self.center = np.array(center)
        self.direction = direction
        length = self.D_s / np.radians(max(direction, 0.5)) * 5e-8
        # self.points = self.createSource(length, direction, angle, num)
        self.points = self.createFulledSource(length, 3e-3, direction, angle, num)
        # self.points = self.createCircleSource(1e-2, 100)

    def createSource(self, length, direction, angle, num):
        v = np.linspace([0, 0, 0], [0, 0, -length], num)
        v = self.rotateX(v, direction)
        v = self.rotateZ(v, angle)

        return self.scale(v)

    def createCircleSource(self, radius, num):
        theta = np.linspace(0, 2 * np.pi, num)
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)

        d = np.zeros(num)

        return self.scale(np.array([x, y, d]).T)

    def createFulledSource(self, length, radius, direction, angle, num):
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
        s = np.array(points).T
        s[0:2] *= arcsec / self.D_s
        s = s.T
        s += np.array([self.center[0], self.center[1], self.D_s])

        return s


class Solver:
    def __init__(self, lens, source) -> None:
        self.lens = lens
        self.source = source
        # self.image = np.empty(len(self.source.points))

    def einsteinRadius(self, D_s):
        D_ls = self.lens.D_ls
        D_l = D_s - D_ls

        return np.sqrt(4 * constants.G.value * self.lens.m / constants.c.value ** 2 * D_ls / D_l / D_s / 3e19) * arcsec

    def einsteinRadiusZ(self, z):
        D_s = model.angular_diameter_distance(z)

        return self.einsteinRadius(D_s)

    def processPoint(self, p_):
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

