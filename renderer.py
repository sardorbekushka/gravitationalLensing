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
        self.ax.set_xlabel('mas')
        self.ax.set_ylabel('mas')
        # self.ax.set_xlim([-0.25, 0.25])
        # self.ax.set_ylim([-4.2, -3.2])
        self.showSource = True
        self.showImage = False
        self.showMagnification = True
        self.showData = True

        self.scatter_image = None
        self.scatter_source = None
        self.lens_circle = None
        self.order = None
        self.real_image = None
        self.real_data = None

        # self.order2 = [val for pair in zip(order, order) for val in pair]
        self.cbar = None

        self.initialize_plot()

        self.data = (rf'mass: {"{:.3e}".format(self.solver.getLensMass())} kg' + '\n' +
                     rf'$\theta_E$: {"{:.3e}".format(self.solver.getEinsteinRadius())} mas' + '\n' +
                     rf'$z_s$: {self.solver.source.z}' + '\n' +
                     r'$D_{ls}:$' + f'{round(self.solver.lens.D_ls, 2)} kpc \n' +
                     rf'jet direction: {round(self.solver.getSourceDirection(), 2)}$^\circ$')

        self.title = rf'$H_0$: {model.H0}, $\Omega_M$: {model.Om0}, $\Omega_0$: {model.Ode0} '
        self.ax.set_title(self.title, fontsize=12)

        bbox = dict(boxstyle='round', fc=g, ec=g/2, alpha=0.3)
        self.text_block = self.ax.text(0.625, 0.99, self.data, fontsize=9, bbox=bbox,
                                       color=g,
                                       horizontalalignment='left', verticalalignment='top',
                                       transform=self.ax.transAxes)

    def updateOrder(self):
        order = np.arange(len(self.solver.getSourcePoints()))
        self.order = order

    def initialize_plot(self):
        self.scatter_image = self.ax.scatter([], [], c=[], cmap='inferno', s=0.5, norm=colors.LogNorm(vmin=0.1, vmax=100))
        self.cbar = plt.colorbar(self.scatter_image, ax=self.ax)
        self.cbar.set_label('magnification')
        self.scatter_source = self.ax.scatter([], [] , c=[], cmap='viridis', s=0.5, alpha=0.01)

        self.lens_circle, = self.ax.plot([], [], color=g, linewidth=1, linestyle=':')
        self.scatter_source.set_zorder(1)
        self.scatter_image.set_zorder(2)
        # self.real_image = self.ax.imshow(plt.imread('src/image.png'),
        #                                  extent=(-2, 2, -7, 2), alpha=0.3)

        # self.real_image = self.ax.imshow(plt.imread('src/image2.png'),
        #                                  extent=(-1, 1, -4.5, 1), alpha=0.3)
        data = readData()
        # self.real_data = self.ax.scatter(data[1], data[3], c='c', s=0.5)
        self.real_data = self.ax.errorbar(data[1], data[3], data[4], data[2], c=[0, 1, 1], fmt='o', elinewidth=0.1, markersize=0.5)
        # self.real_data = self.ax.scatter(data[1], data[3], c=data[5], cmap='inferno', s=0.5, norm=colors.LogNorm(vmin=0.001, vmax=1))


    def generate_filename(self):
        s = self.solver
        m = s.getLensMass()
        pos = s.getLensCenter()
        a = s.getSourceDirection()
        return f'interesting/test/M{"{:.4e}".format(m)}_X{"{:.1e}".format(pos[1])}Y{"{:.1e}".format(pos[0])}_A{"{:.3g}".format(a)}_y1{"{:.2e}".format(s.source.y1)}_S{int(self.showSource)}.png'

    def handleKeyEvent(self, event):
        shift = [0.0, 0]
        shift_size = 0.5
        if event.key == 'right':
            shift = [shift_size, 0]
        elif event.key == 'left':
            shift = [-shift_size, 0]
        elif event.key == 'up':
            shift = [0, shift_size]
        elif event.key == 'down':
            shift = [0, -shift_size]

        self.solver.moveLens(shift)

        if event.key in ['m', 'ь']:
            self.showMagnification = not self.showMagnification
        elif event.key == 'enter':
            plt.savefig(self.generate_filename(), dpi=150)
        elif event.key == '=':
            self.solver.increaseMass()
        elif event.key == '-':
            self.solver.decreaseMass()

        elif event.key in ['h', 'р']:
            self.showSource = not self.showSource
        elif event.key in ['d', 'в']:
            self.showData = not self.showData

        elif event.key in ['w', 'z', 'ц', 'я']:
            direction = 1 if event.key in ['w', 'ц'] else -1
            self.solver.declineSource(direction)

        elif event.key in ['W', 'Z', 'Ц', 'Я']:
            k = 1 if event.key in ['Z', 'Я'] else -1
            self.solver.moveDls(k)

        elif event.key in ['i', 'ш']:
            self.showImage = not self.showImage

        elif event.key in ['u', 'г', 'U', 'Г']:
            k = 1 if event.key in ['u', 'г'] else -1
            self.solver.increaseWidth(k)

        elif event.key == 'r':
            print(self.solver.getEfficiency())

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
            self.scatter_image.set_array(m)
            self.scatter_image.set_cmap('inferno')
            self.scatter_image.set_norm(colors.LogNorm(vmin=0.1, vmax=100))
            self.scatter_source.set_array([10])
            self.scatter_source.set_norm(colors.Normalize(vmin=1, vmax=2))
            self.cbar.set_label('magnification', color=g)

        else:
            self.scatter_image.set_array(self.order)
            self.scatter_image.set_cmap('viridis')
            self.scatter_image.set_norm(colors.Normalize(vmin=self.order[0], vmax=self.order[-1]))
            self.scatter_source.set_array(self.order)
            self.scatter_source.set_cmap('viridis')
            self.scatter_source.set_norm(colors.Normalize(vmin=self.order[0], vmax=self.order[-1]))
            self.cbar.set_label('point order', color=g)

        # self.real_image.set_visible(self.showImage)
        # for line in self.real_data.lines:
        #     line.set_visible(self.showImage)
        for item in self.real_data.lines:
            if hasattr(item, 'set_visible'):
                item.set_visible(self.showImage)
            elif isinstance(item, tuple):  # Если внутри есть вложенные элементы
                for sub_item in item:
                    if hasattr(sub_item, 'set_visible'):
                        sub_item.set_visible(self.showImage)

        self.scatter_source.set_visible(self.showSource)
        self.text_block.set_visible(self.showData)

    def updateData(self):
        self.data = (rf'mass: {"{:.3e}".format(self.solver.getLensMass())} kg' + '\n' +
                     rf'$\theta_E$: {"{:.3e}".format(self.solver.getEinsteinRadius())} mas' + '\n' +
                     rf'$z_s$: {self.solver.source.z}' + '\n' +
                     r'$D_{ls}:$' + f'{round(self.solver.lens.D_ls, 2)} kpc \n' +
                     rf'jet direction: {round(self.solver.getSourceDirection(), 2)}$^\circ$')

    def show(self):
        p, m = self.solver.processImage()
        ll = self.solver.getLensCenter()
        e_an = self.solver.getEinsteinRadius()
        theta = np.linspace(0, 2 * np.pi, 100)
        self.updateOrder()
        self.checkFlags(m)
        self.scatter_image.set_offsets(p.T)
        pp = self.solver.getSourcePoints()
        self.scatter_source.set_offsets(pp)

        x = e_an * np.cos(theta) + ll[0]
        y = e_an * np.sin(theta) + ll[1]
        self.lens_circle.set_data(x, y)

        self.text_block.set_text(self.data)

        self.ax.set_facecolor(g / 10)
        plt.grid(linestyle=':', color=g / 2)
        plt.tick_params(axis='x', colors=g)
        plt.tick_params(axis='y', colors=g)

        plt.draw()

        return self.scatter_image

    def start(self):
        key_id = plt.connect('key_press_event', self.handleKeyEvent)
        mouseMove_id = plt.connect('motion_notify_event', self.handleMouseEvent)
        mouseClick_id = plt.connect('button_press_event', self.handleMouseEvent)

        self.cbar.update_normal(self.show())

        plt.show()

    def gif(self):
        cbar = plt.colorbar(self.show(), ax=self.ax)
        cbar.set_label('ln(magnification)', color=g)
        cbar.ax.tick_params(labelcolor=g)
        for i in range(100):
            plt.savefig(f'images/image{i}.png', dpi=150)
            self.solver.moveLens([5, 0])
            self.show()
