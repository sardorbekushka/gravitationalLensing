from settings import *


class Renderer:
    def __init__(self, solver, ax) -> None:
        self.solver = solver
        self.ax = ax

    def handleEvent(self, event):
        pass

    def show(self):
        p = self.solver.processImage()
        self.ax.scatter(p[0], p[1], s=5, alpha=0.2)
        e_an = self.solver.getEinsteinRadius()
        theta = np.linspace(0, 2 * np.pi, 100)
        self.ax.plot(e_an * np.cos(theta), e_an * np.sin(theta), color=g, linestyle=':')
        pp = self.solver.source.points.T
        self.ax.scatter(pp[0], pp[1], color=[0.9, 0.9, 0.9], alpha=0.1, s=5)

        data = (rf'mass: {"{:.3e}".format(self.solver.getLensMass())} kg' + '\n' +
                rf'$\theta_E$: {"{:.3e}".format(self.solver.getEinsteinRadius())} arcsec' + '\n' +
                r'$D_{l}$: ' + f'{"{:.3e}".format(self.solver.getLensDistance())} kpc \n' +
                rf'$D_s$: {"{:.3e}".format(self.solver.getSourceDistance())} kpc' + '\n' +
                rf'jet direction: {self.solver.getSourceDirection()}$^\circ$')
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
        lim = np.array([-e_an, e_an]) * 2
        self.ax.set_xlim(self.solver.getLensCenter()[0] + lim)
        self.ax.set_ylim(self.solver.getLensCenter()[1] + lim)
        plt.grid(linestyle='--', color=g/2)
        plt.tick_params(axis='x', colors=g)
        plt.tick_params(axis='y', colors=g)
        # plt.legend(loc='best')

        # cid = plt.connect('button_press_event', self.handleEvent)
        plt.show()
        