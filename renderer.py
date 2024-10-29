import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
import numpy as np

from settings import *


class Renderer:
    def __init__(self, solver, ax) -> None:
        self.solver = solver
        self.ax = ax
        self.showSource = True

    def handleKeyEvent(self, event):
        # print(event.key)
        shift = [0.0, 0]
        shift_size = 0.5
        # print(event.key)
        if event.key == 'right':
            shift = [shift_size, 0]
            # print(self.solver.getLensCenter())
        elif event.key == 'left':
            shift = [-shift_size, 0]
        if event.key == 'up':
            shift = [0, shift_size]
        elif event.key == 'down':
            shift = [0, -shift_size]
        self.solver.moveLens(shift)
        if event.key == 'm':
            s = self.solver
            m = s.getLensMass()
            pos = s.getLensCenter()
            a = s.getSourceDirection()
            plt.savefig(f'library/M{"{:.4e}".format(m)}_X{"{:.1e}".format(pos[1])}Y{"{:.1e}".format(pos[0])}_A{"{:.3g}".format(a)}_S{int(self.showSource)}.png', dpi=150)
        angle = 0
        if event.key == 'd':
            angle = 1
            # print(self.solver.source.angle)
        elif event.key == 'a':
            angle = -1
        # self.solver.rotateSource(angle)

        direction = 0
        if event.key == 'w':
            direction = 1
        elif event.key == 'z':
            direction = -1

        # self.solver.turnSource(angle, direction)
        # self.solver.declineSource(direction)

        if event.key == '=':
            self.solver.increaseMass()
            # print(self.solver.getLensMass())
        elif event.key == '-':
            self.solver.decreaseMass()

        if event.key == 'enter':
            self.showSource = not self.showSource

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

        pp = self.solver.source.points.T
        if self.showSource:
            # self.ax.scatter(pp[0], pp[1], c=pp[2], cmap='viridis', alpha=5e-2, s=5)
            self.ax.scatter(pp[0], pp[1], color=g, alpha=1, s=5)

        # print(p, m)
        # print(max(m))
        # m = np.array([min(mm / 10, 1) for mm in m])
        # m = np.log(m)
        # print(m)
        # m = 1
        # color2 = [val for pair in zip(pp[2], pp[2]) for val in pair]
        # color2 = color

        # self.ax.scatter(p[0], p[1], color=[0.9, 0.9, 0.9], s=5, alpha=m)
        scatter = self.ax.scatter(p[0], p[1], c=m, cmap='inferno', s=5, vmin=0, vmax=2)
        # plt.colorbar(label='D_s')

        ll = self.solver.getLensCenter()
        e_an = self.solver.getEinsteinRadius()
        theta = np.linspace(0, 2 * np.pi, 100)

        # self.ax.scatter(ll[0], ll[1], c='w')
        # color = ['r', 'g', 'c']

        self.ax.plot(e_an * np.cos(theta) + ll[0], e_an * np.sin(theta) + ll[1], color=g, linestyle=':')

        # print(pp[0], pp[1], pp[2])
        # print(self.solver.getLensDistance())

        # for i in range(len(pp[2])):
        #     e_an = self.solver.einsteinRadius(pp[2][i])
        #     # print(e_an)
        #     self.ax.plot(e_an * np.cos(theta) + ll[0], e_an * np.sin(theta) + ll[1], color=color[i], linestyle=':')

        data = (rf'mass: {"{:.3e}".format(self.solver.getLensMass())} kg' + '\n' +
                rf'$\theta_E$: {"{:.3e}".format(self.solver.getEinsteinRadius())} arcsec' + '\n' +
                # r'$D_{l}$: ' + f'{"{:.3e}".format(self.solver.getLensDistance())} kpc \n' +
                # rf'$D_s$: {"{:.3e}".format(self.solver.getSourceDistance())} kpc' + '\n' +
                rf'$z_s$: {self.solver.source.z}' + '\n' +
                r'$D_{ls}:$' + f'{self.solver.lens.D_ls} kpc \n' +
                rf'jet direction: {self.solver.getSourceDirection()}$^\circ$')# + '\n' +
                # f'Re = {"{:.3e}".format(self.solver.getLensDistance() * self.solver.getEinsteinRadius() / arcsec)} kpc \n' +
                # f'Dl = {"{:.3e}".format(self.solver.getLensDistance())} kpc')
        title = rf'$H_0$: {model.H0}, $\Omega_M$: {model.Om0}, $\Omega_0$: {model.Ode0} '

        # data = (f'mass: {"{:.3e}".format(self.solver.getLensMass())} kg \n'
        #         f'e angle: {"{:.3e}".format(self.solver.getEinsteinRadius())} \'\' \n' +
        #         f'Dl: {"{:.3e}".format(self.solver.getLensDistance())} kpc \n' +
        #         f'Ds: {"{:.3e}".format(self.solver.getSourceDistance())} kpc \n' +
        #         f'jet direction: {self.solver.getSourceDirection()}')
        # title = f'H0: {model.H0.value}, omegaM: {model.Om0}, omegaL: {model.Ode0}'

        bbox = dict(boxstyle='round', fc=g, ec=g/2, alpha=0.3)
        self.ax.text(0.73, 0.95, data, fontsize=9, bbox=bbox, color=g,
                     horizontalalignment='left', verticalalignment='top',
                     transform=self.ax.transAxes)
        self.ax.set_title(title, color=g)

        self.ax.set_facecolor([0.05, 0.05, 0.1])
        lim = np.array([-e_an, e_an]) * 3
        # self.ax.set_xlim(self.solver.getLensCenter()[0] + lim)
        # self.ax.set_ylim(self.solver.getLensCenter()[1] + lim)

        self.ax.set_xlim(self.solver.getSourceCenter()[0] + lim)
        self.ax.set_ylim(self.solver.getSourceCenter()[1] + lim)
        plt.grid(linestyle='--', color=g/2)
        plt.tick_params(axis='x', colors=g)
        plt.tick_params(axis='y', colors=g)
        # plt.legend(loc='best')
        plt.draw()

        return scatter

    def start(self):
        key_id = plt.connect('key_press_event', self.handleKeyEvent)
        mouseMove_id = plt.connect('motion_notify_event', self.handleMouseEvent)
        mouseClick_id = plt.connect('button_press_event', self.handleMouseEvent)

        cbar = plt.colorbar(self.show(), ax=self.ax)
        cbar.set_label('ln(magnification)', color=g)
        cbar.ax.tick_params(labelcolor=g)  # Цвет меток
        # cbar.ax.yaxis.label.set_color(g)  # Цвет подписи
        # plt.colorbar(self.show(), ax=self.ax)
        plt.show()
