import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
import numpy as np

from settings import *


class Renderer:
    def __init__(self, solver, ax) -> None:
        self.solver = solver
        self.ax = ax
        self.showSource = True
        self.showMagnification = True
        self.showData = True

    def handleKeyEvent(self, event):
        # print(event.key)
        shift = [0.0, 0]
        shift_size = 0.5
        if event.key == 'right':
            shift = [shift_size, 0]
        elif event.key == 'left':
            shift = [-shift_size, 0]
        if event.key == 'up':
            shift = [0, shift_size]
        elif event.key == 'down':
            shift = [0, -shift_size]
        self.solver.moveLens(shift)
        if event.key == 'm' or event.key == 'ь':
            self.showMagnification = not self.showMagnification
        if event.key == 'enter':
            s = self.solver
            m = s.getLensMass()
            pos = s.getLensCenter()
            a = s.getSourceDirection()
            plt.savefig(f'library/M{"{:.4e}".format(m)}_X{"{:.1e}".format(pos[1])}Y{"{:.1e}".format(pos[0])}_A{"{:.3g}".format(a)}_S{int(self.showSource)}.png', dpi=150)

        if event.key == '=':
            self.solver.increaseMass()
        elif event.key == '-':
            self.solver.decreaseMass()

        if event.key == 's' or event.key == 'ы':
            self.showSource = not self.showSource
        # angle = 0
        # if event.key == 'd':
        #     angle = 1
        # elif event.key == 'a':
        #     angle = -1
        # self.solver.rotateSource(angle)

        if event.key == 'w' or event.key == 'z':
            direction = 1 if event.key == 'w' else -1
            self.solver.declineSource(direction)

        self.show()

    def handleMouseEvent(self, event):
        if event.inaxes:
            if event.button is MouseButton.LEFT:
                mouse_pos = (np.array([event.xdata, event.ydata]))
                self.solver.setLens(mouse_pos)

        self.show()

    def show(self):
        self.ax.clear()
        p, m = self.solver.processImage()

        order = np.arange(0, len(self.solver.source.points))

        m = np.log(m)
        order2 = [val for pair in zip(order, order) for val in pair]

        scatter = self.ax.scatter(p[0], p[1], c=m, cmap='inferno', s=5, vmin=-3, vmax=3) if self.showMagnification\
            else self.ax.scatter(p[0], p[1], c=order2, cmap='viridis', s=5)

        pp = self.solver.source.points.T
        if self.showSource:
            if self.showMagnification:
                self.ax.scatter(pp[0], pp[1], color=g, s=5)
            else:
                self.ax.scatter(pp[0], pp[1], c=order, cmap='viridis', alpha=0.2, s=5)

        ll = self.solver.getLensCenter()
        e_an = self.solver.getEinsteinRadius()
        theta = np.linspace(0, 2 * np.pi, 100)

        self.ax.plot(e_an * np.cos(theta) + ll[0], e_an * np.sin(theta) + ll[1], color=g, linestyle=':')

        data = (rf'mass: {"{:.3e}".format(self.solver.getLensMass())} kg' + '\n' +
                rf'$\theta_E$: {"{:.3e}".format(self.solver.getEinsteinRadius())} arcsec' + '\n' +
                rf'$z_s$: {self.solver.source.z}' + '\n' +
                r'$D_{ls}:$' + f'{self.solver.lens.D_ls} kpc \n' +
                rf'jet direction: {round(self.solver.getSourceDirection(), 2)}$^\circ$' + 'добавить центр линзы')# + '\n' +
                # f'Re = {"{:.3e}".format(self.solver.getLensDistance() * self.solver.getEinsteinRadius() / arcsec)} kpc \n' +
                # f'Dl = {"{:.3e}".format(self.solver.getLensDistance())} kpc')
        title = rf'$H_0$: {model.H0}, $\Omega_M$: {model.Om0}, $\Omega_0$: {model.Ode0} '

        bbox = dict(boxstyle='round', fc=g, ec=g/2, alpha=0.3)
        textpos = [0.73, 0.95]
        if self.showData:
            self.ax.text(0.73, 0.95, data, fontsize=9, bbox=bbox, color=g,
                         horizontalalignment='left', verticalalignment='top',
                         transform=self.ax.transAxes)

        self.ax.set_title(title, color=g)

        self.ax.set_facecolor([0.05, 0.05, 0.1])
        limx = np.array(lim[0]) * 1e-3 if lim else np.array([-e_an, e_an]) * 3
        limy = np.array(lim[1]) * 1e-3 if lim else np.array([4 * e_an, -2 * e_an]) * np.sign(self.solver.getSourceDirection())

        self.ax.set_xlim(self.solver.getSourceCenter()[0] + limx)
        self.ax.set_ylim(self.solver.getSourceCenter()[1] + limy)
        plt.grid(linestyle='--', color=g/2)
        plt.tick_params(axis='x', colors=g)
        plt.tick_params(axis='y', colors=g)
        plt.draw()

        return scatter

    def start(self):
        key_id = plt.connect('key_press_event', self.handleKeyEvent)
        mouseMove_id = plt.connect('motion_notify_event', self.handleMouseEvent)
        mouseClick_id = plt.connect('button_press_event', self.handleMouseEvent)

        cbar = plt.colorbar(self.show(), ax=self.ax)
        if self.showMagnification:
            cbar.set_label('ln(magnification)', color=g)
        else:
            cbar.remove()
        cbar.ax.tick_params(labelcolor=g)
        plt.show()

    def gif(self):
        cbar = plt.colorbar(self.show(), ax=self.ax)
        cbar.set_label('ln(magnification)', color=g)
        cbar.ax.tick_params(labelcolor=g)
        # while self.solver.getLensCenter[0] < 2e-3:
        for i in range(100):
            plt.savefig(f'images/image{i}.png', dpi=150)
            self.solver.moveLens([5, 0])
            self.show()

        # для галлереи аналогично можно в цикле посоздавать изображения. условно for angle in range(A1, A2):
        #                                                                            for m in range(M1, M2):
        #                                                                                 for D_ls in range(D1, D2):
        #                                                                                      ... и координаты