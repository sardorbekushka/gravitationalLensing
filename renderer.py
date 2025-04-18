import matplotlib.pyplot as plt
from PIL.ImageOps import scale
from matplotlib.backend_bases import MouseButton
import numpy as np
from solver import*
# from reverseSolver import source
from settings import *
from scipy.ndimage import gaussian_filter

class Renderer:
    def __init__(self, solver) -> None:
        figsize = ((lim[0][1] - lim[0][0]) / (lim[1][1] - lim[1][0]) * 7 * 1.24, 7) if lim else (5, 7)
        plt.figure(figsize=figsize, facecolor='black')
        plt.style.use('dark_background')
        ax = plt.axes()
        self.solver = solver
        self.ax = ax
        # self.ax.set_xlim([-2, 2])
        # self.ax.set_ylim([-5, 2])
        self.ax.set_xlim(lim[0])
        self.ax.set_ylim(lim[1])
        self.ax.set_xlabel('mas')
        self.ax.set_ylabel('mas')
        # self.ax.set_xlim([-0.25, 0.25])
        # self.ax.set_ylim([-4.2, -3.2])
        self.showSource = False
        self.showImage = True
        self.showMagnification = True
        self.showData = True

        self.scatter_image = None
        self.scatter_source = None
        self.lens_circle = None
        self.lens_pos = None
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
                     rf'jet direction: {round(self.solver.getSourceDirection(), 5)}$^\circ$')

        self.title = rf'$H_0$: {model.H0}, $\Omega_M$: {model.Om0}, $\Omega_0$: {model.Ode0} '
        self.ax.set_title(self.title, fontsize=12)

        bbox = dict(boxstyle='round', fc=g, ec=g/2, alpha=0.3)
        self.text_block = self.ax.text(0.64, 0.99, self.data, fontsize=8, bbox=bbox,
                                       color=g,
                                       horizontalalignment='left', verticalalignment='top',
                                       transform=self.ax.transAxes)

    def updateOrder(self):
        order = np.arange(len(self.solver.getSourcePoints()))
        self.order = order

    def initialize_plot(self):
        self.lens_pos = self.ax.scatter([], [], marker='x', c='r')
        #
        # self.scatter_image = self.ax.scatter([], [], c=[], cmap='inferno', s=5, norm=colors.LogNorm(vmin=0.1, vmax=100))
        # self.cbar = plt.colorbar(self.scatter_image, ax=self.ax)
        # self.cbar.set_label('magnification')

        self.scatter_image = self.ax.scatter([], [], c=[], cmap='inferno', s=5, norm=colors.LogNorm(vmin=0.1 * flux, vmax=100 * flux))
        self.cbar = plt.colorbar(self.scatter_image, ax=self.ax)
        self.cbar.set_label('flux, Jy')

        self.scatter_source = self.ax.scatter([], [] , c=[], cmap='viridis', s=0.5, alpha=0.2)

        self.lens_circle, = self.ax.plot([], [], color=g, linewidth=1, linestyle=':')


        self.scatter_source.set_zorder(1)
        self.scatter_image.set_zorder(2)

        # data_old15kHz.txt = readData()
        data = self.solver.real_data
        self.real_data = self.ax.errorbar(data[1], data[3], data[4], data[2], c=[0, 1, 1], fmt='o', elinewidth=0.15, markersize=0.5)
        # self.real_data = self.ax.scatter(data_old15kHz.txt[1], data_old15kHz.txt[3], c=data_old15kHz.txt[5], cmap='inferno', s=0.5, norm=colors.LogNorm(vmin=0.001, vmax=1.6))

    def generate_filename(self):
        s = self.solver
        m = s.getLensMass()
        pos = s.getLensCenter()
        a = s.getSourceDirection()
        return f'top10/M{"{:.4e}".format(m)}_X{"{:.4e}".format(pos[0])}Y{"{:.4e}".format(pos[1])}_A{"{:.4g}".format(a)}_D1_{"{:.2g}".format(s.source.d1)}_X1_{"{:.4g}".format(s.source.x1)}.png'

    def generate_filename_(self):
        e = self.solver.getEfficiency()
        return self.generate_filename()[:-4] + f'_{"{:.4g}".format(e[0])}_{"{:.4g}".format(e[1])}' + '.png'

    def handleKeyEvent(self, event):
        shift_size = 0.5
        shift_map = {
            'right': [shift_size, 0],
            'left': [-shift_size, 0],
            'up': [0, shift_size],
            'down': [0, -shift_size]
        }

        shift = shift_map.get(event.key, [0.0, 0])

        if event.key in shift_map:
            self.solver.moveLens(shift)

        toggle_map = {
            'm': 'showMagnification', 'ь': 'showMagnification',
            'h': 'showSource', 'р': 'showSource',
            'd': 'showData', 'в': 'showData',
            'i': 'showImage', 'ш': 'showImage'
        }

        if event.key in toggle_map:
            setattr(self, toggle_map[event.key], not getattr(self, toggle_map[event.key]))

        action_map = {
            'enter': lambda: plt.savefig(self.generate_filename_(), dpi=150),
            '+': lambda: self.solver.changeMass(1),
            '_': lambda: self.solver.changeMass(-1),
            'e': lambda: print(self.solver.getEfficiency()),
            'E': lambda: print(self.solver.getEfficiency())
        }

        if event.key in action_map:
            action_map[event.key]()

        decline_map = {
            'w': 1, 'z': -1, 'ц': 1, 'я': -1
        }

        if event.key in decline_map:
            direction = decline_map[event.key]
            self.solver.declineSource(direction)

        move_lens_map = {
            'W': -1, 'Z': 1, 'Ц': -1, 'Я': 1
        }

        if event.key in move_lens_map:
            shift = move_lens_map[event.key]
            self.solver.moveDls(shift)

        width_change_map = {
            'u': 1, 'г': 1, 'U': -1, 'Г': -1
        }

        if event.key in width_change_map:
            k = width_change_map[event.key]
            self.solver.changeWidth(k)

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
            self.scatter_image.set_norm(colors.LogNorm(vmin=0.1 * flux, vmax=100 * flux))
            self.scatter_source.set_array([10])
            self.scatter_source.set_norm(colors.Normalize(vmin=1, vmax=2))
            self.cbar.set_label('flux, Jy', color=g)

        else:
            self.scatter_image.set_array(self.order)
            self.scatter_image.set_cmap('viridis')
            self.scatter_image.set_norm(colors.Normalize(vmin=self.order[0], vmax=self.order[-1]))
            self.scatter_source.set_array(self.order)
            self.scatter_source.set_cmap('viridis')
            self.scatter_source.set_norm(colors.Normalize(vmin=self.order[0], vmax=self.order[-1]))
            self.cbar.set_label('point order', color=g)

        for item in self.real_data.lines:
            if hasattr(item, 'set_visible'):
                item.set_visible(self.showImage)
            elif isinstance(item, tuple):
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
                     rf'jet direction: {round(self.solver.getSourceDirection(), 5)}$^\circ$')

    def show(self):
        p, m = self.solver.processImage()
        m *= flux
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
        self.lens_pos.set_offsets(ll)

        self.text_block.set_text(self.data)


        plt.draw()

        return self.scatter_image

    def show_blur(self):
        self.ax.clear()
        p, m = self.solver.processImage()
        m *= flux

        N = 400
        M = 700
        sx = 0.48
        sy = 1.06

        image = np.zeros((M, N), dtype=np.float64)
        xx = np.round((p[0] + 2) / 4 * N).astype(int)
        yy = np.round((p[1] + 5) / 7 * M).astype(int)

        image[yy, xx] = m
        image = np.nan_to_num(image)

        sigma_x = sx / 4 * N  # Параметр размытия (чем больше, тем сильнее размытие)
        sigma_y = sy / 7 * M  # Параметр размытия (чем больше, тем сильнее размытие)
        sigma = np.array([sigma_y, sigma_x]) / 2.35482

        blurred_image = gaussian_filter(image, sigma=sigma)
        blurred_image[blurred_image <= 0] = 1e-10
        self.scatter_image = plt.imshow(blurred_image, cmap='hot', extent=[-2, 2, -5, 2], origin='lower', vmin=0, vmax=1e-3)
        # self.scatter_image = plt.imshow(blurred_image, cmap='hot', extent=[-2, 2, -5, 2], origin='lower', norm=colors.LogNorm(vmin=1e-5, vmax=1e-3))

        plt.ylabel('mas')
        plt.xlabel('mas')

        lens_center = self.solver.lens.center
        circle = plt.Circle(lens_center, self.solver.getEinsteinRadius(), color='gray', linestyle='dotted',
                            fill=False)  # Центр (0.5, 0.5), радиус 0.2
        self.ax.add_patch(circle)
        self.ax.scatter(lens_center[0], lens_center[1], marker='x', c='r')

        # plt.show()
        plt.draw()
        return self.scatter_image

    def start(self):

        self.ax.set_facecolor(g / 10)
        plt.grid(linestyle=':', color=g / 2)
        plt.tick_params(axis='x', colors=g)
        plt.tick_params(axis='y', colors=g)
        key_id = plt.connect('key_press_event', self.handleKeyEvent)
        mouseMove_id = plt.connect('motion_notify_event', self.handleMouseEvent)
        mouseClick_id = plt.connect('button_press_event', self.handleMouseEvent)

        self.cbar.update_normal(self.show())

        plt.show()
