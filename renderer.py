import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
import numpy as np
from settings import *


class Renderer:
    def __init__(self, solver, ax) -> None:
        self.solver = solver
        self.ax = ax
        self.ax.set_xlim([-2, 2])
        self.ax.set_ylim([-5, 2])
        self.showSource = True
        self.showMagnification = True
        self.showData = True

        self.scatter_points = None
        self.scatter_source = None
        self.lens_circle = None
        order = np.arange(1)
        self.order = order
        self.order2 = [val for pair in zip(order, order) for val in pair]

        self.initialize_plot()

        self.data = (rf'mass: {"{:.3e}".format(self.solver.getLensMass())} kg' + '\n' +
                     rf'$\theta_E$: {"{:.3e}".format(self.solver.getEinsteinRadius())} mas' + '\n' +
                     rf'$z_s$: {self.solver.source.z}' + '\n' +
                     r'$D_{ls}:$' + f'{self.solver.lens.D_ls} kpc \n' +
                     rf'jet direction: {round(self.solver.getSourceDirection(), 2)}$^\circ$')

        self.title = rf'$H_0$: {model.H0}, $\Omega_M$: {model.Om0}, $\Omega_0$: {model.Ode0} '
        self.ax.set_title(self.title, fontsize=12)

        bbox = dict(boxstyle='round', fc=g, ec=g/2, alpha=0.3)
        self.text_block = self.ax.text(0.625, 0.99, self.data, fontsize=9, bbox=bbox,
                                       color=g,
                                       horizontalalignment='left', verticalalignment='top',
                                       transform=self.ax.transAxes)

    def updateOrder(self):
        order = np.arange(len(self.solver.points))
        self.order = order

    def initialize_plot(self):
        self.scatter_points = self.ax.scatter([], [], c=[], cmap='inferno', s=0.5, norm=colors.LogNorm(vmin=0.1, vmax=10))
        pp = self.solver.source.points.T

        self.scatter_source = self.ax.scatter(pp[0], pp[1] , c=self.order, cmap='viridis', s=0.5, alpha=0.5)

        self.lens_circle, = self.ax.plot([], [], color=g, linestyle=':')
        self.scatter_source.set_zorder(1)
        self.scatter_points.set_zorder(2)

    def generate_filename(self):
        s = self.solver
        m = s.getLensMass()
        pos = s.getLensCenter()
        a = s.getSourceDirection()
        return f'library/M{"{:.4e}".format(m)}_X{"{:.1e}".format(pos[1])}Y{"{:.1e}".format(pos[0])}_A{"{:.3g}".format(a)}_S{int(self.showSource)}.png'

    def handleKeyEvent(self, event):
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
            plt.savefig(self.generate_filename(), dpi=150)
        if event.key == '=':
            self.solver.increaseMass()
        elif event.key == '-':
            self.solver.decreaseMass()

        if event.key == 'h' or event.key == 'р':
            self.showSource = not self.showSource
        if event.key == 'd' or event.key == 'в':
            self.showData = not self.showData

        if event.key == 'w' or event.key == 'z':
            direction = 1 if event.key == 'w' else -1
            self.solver.declineSource(direction)

        if event.key == 'W' or event.key == 'Z':
            k = 1 if event.key == 'Z' else -1
            self.solver.moveDls(k)

        self.updateData()
        self.show()

    def handleMouseEvent(self, event):
        if event.inaxes:
            if event.button is MouseButton.LEFT:
                mouse_pos = (np.array([event.xdata, event.ydata]))
                self.solver.setLens(mouse_pos)

        self.show()

    def checkFlags(self, m):
        if self.showMagnification:
            self.scatter_points.set_array(m)
            self.scatter_points.set_cmap('inferno')
            self.scatter_source.set_facecolor(g)
            self.scatter_points.set_norm(colors.LogNorm(vmin=0.1, vmax=10))

        else:
            self.scatter_points.set_array(self.order)
            self.scatter_points.set_cmap('viridis')
            self.scatter_source.set_array(self.order)
            self.scatter_source.set_cmap('viridis')
            self.scatter_points.set_norm(colors.Normalize(vmin=self.order[0], vmax=self.order[-1]))

        self.scatter_source.set_visible(self.showSource)
        self.text_block.set_visible(self.showData)

    def updateData(self):
        self.data = (rf'mass: {"{:.3e}".format(self.solver.getLensMass())} kg' + '\n' +
                     rf'$\theta_E$: {"{:.3e}".format(self.solver.getEinsteinRadius())} mas' + '\n' +
                     rf'$z_s$: {self.solver.source.z}' + '\n' +
                     r'$D_{ls}:$' + f'{self.solver.lens.D_ls} kpc \n' +
                     rf'jet direction: {round(self.solver.getSourceDirection(), 2)}$^\circ$')
        pp = self.solver.source.points
        self.scatter_source.set_offsets(pp.T)

    def show(self):
        p, m = self.solver.processImage()
        ll = self.solver.getLensCenter()
        e_an = self.solver.getEinsteinRadius()
        theta = np.linspace(0, 2 * np.pi, 100)
        self.updateOrder()
        self.checkFlags(m)

        self.scatter_points.set_offsets(p.T)
        pp = self.solver.source.points
        self.scatter_source.set_offsets(pp.T)
        # print(pp)
        x = e_an * np.cos(theta) + ll[0]
        y = e_an * np.sin(theta) + ll[1]
        self.lens_circle.set_data(x, y)

        self.text_block.set_text(self.data)

        self.ax.set_facecolor(g / 10)
        plt.grid(linestyle='--', color=g / 2)
        plt.tick_params(axis='x', colors=g)
        plt.tick_params(axis='y', colors=g)

        plt.draw()

        return self.scatter_points

    def start(self):
        key_id = plt.connect('key_press_event', self.handleKeyEvent)
        mouseMove_id = plt.connect('motion_notify_event', self.handleMouseEvent)
        mouseClick_id = plt.connect('button_press_event', self.handleMouseEvent)

        cbar = plt.colorbar(self.show(), ax=self.ax)

        cbar.set_label('magnification', color=g)

        cbar.ax.tick_params(labelcolor=g)
        plt.show()

    def gif(self):
        cbar = plt.colorbar(self.show(), ax=self.ax)
        cbar.set_label('ln(magnification)', color=g)
        cbar.ax.tick_params(labelcolor=g)
        for i in range(100):
            plt.savefig(f'images/image{i}.png', dpi=150)
            self.solver.moveLens([5, 0])
            self.show()
