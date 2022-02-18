import numpy as np

from dimension_reduction.figures.style import colors


class fig_delta_phi_z(object):

    def __init__(self, ax, settings):

        z = settings.z

        indices_restr = np.abs(z) < 50

        z_restr = z[indices_restr]

        xlabel_phase_difference_z = r'$z \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        ylabel_phase_difference_z = r'$\Delta \, \phi\; / \; \pi$'
        
        # self.line_phase_sum_z,        = ax.plot(settings.z, np.zeros_like(settings.z), linewidth=1, linestyle='-', color=colors.wet_asphalt, label=r'$(\varphi_1 + \varphi_2)/2$')
        # self.line_phase_difference_z, = ax.plot(settings.z, np.zeros_like(settings.z), linewidth=1, linestyle='-', color=colors.peter_river, label=r'$(\varphi_1 - \varphi_2)/2$')

        self.line_delta_phi_selected_of_z_restr,  = ax.plot(z_restr, np.zeros_like(z_restr), linewidth=1.00, linestyle='--', color=colors.peter_river, label=r'single run')
        self.line_delta_phi_circ_mean_of_z_restr, = ax.plot(z_restr, np.zeros_like(z_restr), linewidth=1.00, linestyle='-',  color=colors.wet_asphalt, label=r'circ mean')

        self.indices_restr = indices_restr


        ax.set_xlabel(xlabel_phase_difference_z)
        ax.set_ylabel(ylabel_phase_difference_z)
        
        ax.set_xticks(settings.z_ticks)

        ax.set_yticks([-0.5, 0, 0.5], minor=False)
        ax.set_yticks([-0.75, -0.25, 0.25, 0.75], minor=True)

        ax.set_xlim(settings.z_min, settings.z_max)
        ax.set_ylim(-0.8, 0.8)
        
        # ax.grid(b=True, which='major', color=colors.color_gridlines_major, linestyle='-', linewidth=0.5)

        ax.grid(b=True, which='major', color=colors.color_gridlines_major, linestyle='-', linewidth=0.5)
        ax.grid(b=True, which='minor', color=colors.color_gridlines_minor, linestyle='-', linewidth=0.5, alpha=0.2)
        
        ax.legend(loc='upper right', bbox_to_anchor=(1.05, 1.225), fancybox=settings.fancybox, framealpha=settings.framealpha, ncol=2)
    
    
    def update(self, delta_phi_circ_mean_of_z, delta_phi_x1_x2_of_z_psf_selected):

        indices_restr = self.indices_restr

        self.line_delta_phi_circ_mean_of_z_restr.set_ydata(delta_phi_circ_mean_of_z[indices_restr]/np.pi)
        self.line_delta_phi_selected_of_z_restr.set_ydata(delta_phi_x1_x2_of_z_psf_selected[indices_restr]/np.pi)