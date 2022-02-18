import numpy as np

from numpy import zeros_like

from numpy import pi

from .. style import colors


class fig_phase_z_eff(object):

    def __init__(self, ax, settings):


        self.hbar = settings.hbar

        self.m_atom = settings.m_atom
    

        self.line_phase_z_eff, = ax.plot(settings.z, zeros_like(settings.z), linewidth=1, linestyle='-', color=colors.wet_asphalt, label='effective')
        self.line_phase_z, = ax.plot(settings.z, zeros_like(settings.z), linewidth=1, linestyle='--', color=colors.peter_river, label='$x=y=0$')

        ax.set_xlim(settings.z_min, settings.z_max)
        ax.set_ylim(-1.1, 1.1)
        
        
        ax.set_xlabel(settings.xlabel_density_z)
        
        ax.set_xticks(settings.z_ticks)
        
        ax.grid(b=True, which='major', color=settings.color_gridlines_major, linestyle='-', linewidth=0.5)
        
        # ax.set_ylabel(r'$\varphi \, / \, \pi$')
        ax.set_ylabel(r'$\cos(\varphi)$')

        ax.legend(loc='lower right', bbox_to_anchor=(1.0, 0.0), fancybox=settings.fancybox, framealpha=settings.framealpha, ncol=1)
        
        
        
        ax_V_z = ax.twinx()
    
        self.line_V_z, = ax_V_z.plot(settings.z, zeros_like(settings.z), linewidth=1, linestyle='-', color=colors.sun_flower)
        
        
        # -----------------------------------------------------------------------------------------
        # V_z_comparison_1 = 0.5 * self.m_atom * settings.vec_omega_z_comparison[0]**2 * (settings.z*1e-6)**2
        # V_z_comparison_2 = 0.5 * self.m_atom * settings.vec_omega_z_comparison[1]**2 * (settings.z*1e-6)**2
        #
        # label_V_z_comparison = '$({0:d}, {1:d})\,$Hz'.format(settings.vec_nu_z_comparison[0], settings.vec_nu_z_comparison[1])
        #
        # ax_V_z.plot(settings.z, V_z_comparison_1 / (self.hbar * 2.0 * np.pi * 1000), linewidth=1, linestyle='--', color=colors.sun_flower, label=label_V_z_comparison)
        # ax_V_z.plot(settings.z, V_z_comparison_2 / (self.hbar * 2.0 * np.pi * 1000), linewidth=1, linestyle='--', color=colors.sun_flower)
        # -----------------------------------------------------------------------------------------

        ax_V_z.set_xlim(settings.z_min, settings.z_max)
        ax_V_z.set_ylim(settings.V_min, settings.V_max)
        
        ax_V_z.set_ylabel(settings.ylabel_V_x_y_z)
        
        # ax_V_z.legend(loc='upper right', bbox_to_anchor=(1.0, 1.2), fancybox=settings.fancybox, framealpha=settings.framealpha, ncol=2)
    
    
    def update(self, phase_z_eff, phase_z, V_z):
        
        # print(np.min(density_z))
        # print(np.max(density_z))
        
        scaling_V = self.hbar * 2 * pi * 1000
        
        V_z = V_z / scaling_V

        # self.line_phase_z_eff.set_ydata(phase_z_eff / np.pi)
        # self.line_phase_z.set_ydata(phase_z / np.pi)

        self.line_phase_z_eff.set_ydata(np.cos(phase_z_eff))
        self.line_phase_z.set_ydata(np.cos(phase_z))

        self.line_V_z.set_ydata(V_z)
