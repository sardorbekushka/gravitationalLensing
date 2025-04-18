import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backend_bases import MouseButton
from matplotlib import colors, rc
from scipy.ndimage import gaussian_filter
from solver import Solver
from settings import g, flux, model


class Renderer:
    def __init__(self, solver: Solver) -> None:
        self.solver = solver
        self._setup_plot_config()
        self._initialize_plot_components()
        self._setup_ui_state()
        self._create_info_panel()

    def _setup_plot_config(self):
        """Configure matplotlib settings and figure appearance"""
        rc('text', usetex=True)
        rc('text.latex', preamble=(
            r'\usepackage[utf8]{inputenc}'
            r'\usepackage[russian]{babel}'
        ))
        rc('axes', edgecolor=g)

        self.limits = [[-2, 2], [-5, 2]]
        aspect_ratio = (self.limits[0][1] - self.limits[0][0]) / (self.limits[1][1] - self.limits[1][0])
        figsize = (aspect_ratio * 7 * 1.24, 7)

        plt.figure(figsize=figsize, facecolor='black')
        plt.style.use('dark_background')
        self.ax = plt.axes()
        self.ax.set_xlim(self.limits[0])
        self.ax.set_ylim(self.limits[1])
        self.ax.set_xlabel('mas')
        self.ax.set_ylabel('mas')

    def _initialize_plot_components(self):
        """Initialize all visualization components"""
        # Lens visualization
        self.lens_pos = self.ax.scatter([], [], marker='x', c='r')
        self.lens_circle, = self.ax.plot([], [], color=g, linewidth=1, linestyle=':')

        # Image and source visualization
        self.scatter_image = self.ax.scatter(
            [], [], c=[], cmap='inferno', s=5,
            norm=colors.LogNorm(vmin=0.1 * flux, vmax=100 * flux)
        )
        self.scatter_source = self.ax.scatter(
            [], [], c=[], cmap='viridis', s=0.5, alpha=0.2
        )

        # Color bar setup
        self.cbar = plt.colorbar(self.scatter_image, ax=self.ax)
        self.cbar.set_label('flux, Jy')

        # Layer ordering
        self.scatter_source.set_zorder(1)
        self.scatter_image.set_zorder(2)

    def _setup_ui_state(self):
        """Initialize UI state flags"""
        self.show_source = False
        self.show_image = True
        self.show_magnification = True
        self.show_data = True
        self.order = None

    def _create_info_panel(self):
        """Create the information text panel"""
        solver = self.solver
        self.data = (
            f'mass: {solver.get_lens_mass:.3e} kg\n'
            f'$\\theta_E$: {solver.get_einstein_radius:.3e} mas\n'
            f'$z_s$: {solver.source.z}\n'
            f'$D_{{ls}}:$ {round(solver.lens.D_ls, 2)} kpc\n'
            f'jet direction: {round(solver.get_source_direction, 5)}$^\\circ$'
        )

        self.title = (
            f'$H_0$: {model.H0}, '
            f'$\\Omega_M$: {model.Om0}, '
            f'$\\Omega_0$: {model.Ode0}'
        )
        self.ax.set_title(self.title, fontsize=12)

        bbox = dict(boxstyle='round', fc=g, ec=g / 2, alpha=0.3)
        self.text_block = self.ax.text(
            0.64, 0.99, self.data, fontsize=8, bbox=bbox,
            color=g, horizontalalignment='left',
            verticalalignment='top', transform=self.ax.transAxes
        )

    def update_order(self):
        """Update the point ordering for visualization"""
        self.order = np.arange(len(self.solver.get_source_points))

    def _update_visualization_state(self, magnification):
        """Update visualization based on current state flags"""
        if self.show_magnification:
            self.scatter_image.set_array(magnification)
            self.scatter_image.set_cmap('inferno')
            self.scatter_image.norm = colors.LogNorm(vmin=0.1 * flux, vmax=100 * flux)
            self.cbar.set_label('flux, Jy', color=g)
        else:
            self.scatter_image.set_array(self.order)
            self.scatter_image.set_cmap('viridis')
            self.scatter_image.norm = colors.Normalize(
                vmin=self.order[0], vmax=self.order[-1]
            )
            self.cbar.set_label('point order', color=g)

        self.scatter_source.set_visible(self.show_source)
        self.text_block.set_visible(self.show_data)

    def update_data(self):
        """Update the information panel data"""
        solver = self.solver
        self.data = (
            f'mass: {solver.get_lens_mass:.3e} kg\n'
            f'$\\theta_E$: {solver.get_einstein_radius:.3e} mas\n'
            f'$z_s$: {solver.source.z}\n'
            f'$D_{{ls}}:$ {round(solver.lens.D_ls, 2)} kpc\n'
            f'jet direction: {round(solver.get_source_direction, 5)}$^\\circ$'
        )
        self.text_block.set_text(self.data)

    def show(self):
        """Update and redraw the visualization"""
        positions, magnification = self.solver.process_image()
        magnification *= flux

        # Update point ordering
        self.update_order()

        # Update visualization elements
        self.scatter_image.set_offsets(positions)
        self.scatter_source.set_offsets(self.solver.get_source_points)
        self._update_visualization_state(magnification)

        # Update lens visualization
        center = self.solver.get_lens_center
        radius = self.solver.get_einstein_radius
        theta = np.linspace(0, 2 * np.pi, 100)
        x = radius * np.cos(theta) + center[0]
        y = radius * np.sin(theta) + center[1]

        self.lens_circle.set_data(x, y)
        self.lens_pos.set_offsets([center])

        plt.draw()
        return self.scatter_image

    def _setup_event_handlers(self):
        """Connect UI event handlers"""
        plt.connect('key_press_event', self._handle_key_event)
        plt.connect('motion_notify_event', self._handle_mouse_event)
        plt.connect('button_press_event', self._handle_mouse_event)

    def _handle_key_event(self, event):
        """Handle keyboard input events"""
        # Movement controls
        movement = {
            'right': [0.5, 0], 'left': [-0.5, 0],
            'up': [0, 0.5], 'down': [0, -0.5]
        }.get(event.key)

        if movement:
            self.solver.move_lens(movement)

        # Toggle controls
        toggles = {
            'm': 'show_magnification', 'ь': 'show_magnification',
            'h': 'show_source', 'р': 'show_source',
            'd': 'show_data', 'в': 'show_data',
            'i': 'show_image', 'ш': 'show_image'
        }.get(event.key)

        if toggles:
            setattr(self, toggles, not getattr(self, toggles))

        # Mass controls
        if event.key in {'+', '_'}:
            self.solver.change_mass(1 if event.key == '+' else -1)

        # Source controls
        if event.key in {'w', 'z', 'ц', 'я'}:
            direction = 1 if event.key in {'w', 'ц'} else -1
            self.solver.decline_source(direction)

        # Distance controls
        if event.key in {'W', 'Z', 'Ц', 'Я'}:
            direction = -1 if event.key in {'W', 'Ц'} else 1
            self.solver.move_Dls(direction)

        self.update_data()
        self.show()

    def _handle_mouse_event(self, event):
        """Handle mouse events"""
        if event.inaxes and event.button is MouseButton.LEFT:
            self.solver.set_lens_center([event.xdata, event.ydata])
            self.show()

    def start(self):
        """Start the interactive visualization"""
        self.ax.set_facecolor(g / 10)
        plt.grid(linestyle=':', color=g / 2)
        plt.tick_params(axis='x', colors=g)
        plt.tick_params(axis='y', colors=g)

        self._setup_event_handlers()
        self.cbar.update_normal(self.show())
        plt.show()

    def show_blur(self):
        """Show blurred version of the image"""
        self.ax.clear()
        positions, magnification = self.solver.process_image()
        magnification *= flux

        # Image dimensions and scaling
        width, height = 400, 700
        x_scale = width / (self.limits[0][1] - self.limits[0][0])
        y_scale = height / (self.limits[1][1] - self.limits[1][0])

        # Create and blur image
        image = np.zeros((height, width), dtype=np.float64)
        x_coords = ((positions[:, 0] - self.limits[0][0]) * x_scale).astype(int)
        y_coords = ((positions[:, 1] - self.limits[1][0]) * y_scale).astype(int)

        image[y_coords, x_coords] = magnification
        image = np.nan_to_num(image)

        # Gaussian blur parameters
        sigma_x = 0.48 / 4 * width / 2.35482
        sigma_y = 1.06 / 7 * height / 2.35482
        blurred = gaussian_filter(image, sigma=[sigma_y, sigma_x])
        blurred[blurred <= 0] = 1e-10

        # Display blurred image
        self.scatter_image = self.ax.imshow(
            blurred, cmap='hot', extent=self.limits[0] + self.limits[1],
            origin='lower', vmin=0, vmax=1e-3
        )

        # Add lens indicator
        center = self.solver.get_lens_center
        self.ax.add_patch(plt.Circle(
            center, self.solver.get_einstein_radius(),
            color='gray', linestyle=':', fill=False
        ))
        self.ax.scatter(center[0], center[1], marker='x', c='r')

        plt.draw()
        return self.scatter_image