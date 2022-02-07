from numpy import zeros_like

from numpy import pi

from ..style import colors


class fig_density_y_eff(object):
    
    
    def __init__(self, ax, settings_graphics):

        self.hbar = settings_graphics.hbar
        self.m_atom = settings_graphics.m_atom

        self.line_density_y_eff, = ax.plot(settings_graphics.y, zeros_like(settings_graphics.y), linewidth=1, linestyle='-', color=colors.wet_asphalt)
        # self.line_density_y_x1, = ax.plot(settings_graphics.y, zeros_like(settings_graphics.y), linewidth=1, linestyle='-', color=colors.wet_asphalt)

        ax.set_xlim(settings_graphics.z_min, settings_graphics.z_max)
        ax.set_ylim(0, settings_graphics.density_y_effective_max)
        
        ax.set_xlabel(settings_graphics.xlabel_density_y)
        
        ax.set_xticks(settings_graphics.y_ticks)
        
        ax.grid(b=True, which='major', color=settings_graphics.color_gridlines_major, linestyle='-', linewidth=0.5)
        
        ax.set_ylabel(settings_graphics.ylabel_density_x_y_z_effective)
        
        ax.set_anchor('W')
        
        
        
        ax_V_y = ax.twinx()
        
        self.line_V_y, = ax_V_y.plot(settings_graphics.y, zeros_like(settings_graphics.y), linewidth=1, linestyle='-', color=colors.sun_flower)
        
        
        
        # -----------------------------------------------------------------------------------------
        # V_y_comparison_1 = 0.5 * self.m_atom * settings_graphics.vec_omega_y_comparison[0]**2 * (settings_graphics.y*1e-6)**2
        # V_y_comparison_2 = 0.5 * self.m_atom * settings_graphics.vec_omega_y_comparison[1]**2 * (settings_graphics.y*1e-6)**2
        #
        # label_V_y_comparison = '$({0:1.1f}, {1:1.1f})\,$kHz'.format(settings_graphics.vec_nu_y_comparison[0]/1e3, settings_graphics.vec_nu_y_comparison[1]/1e3)
        #
        # ax_V_y.plot(settings_graphics.y, V_y_comparison_1 / (self.hbar * 2.0 * pi * 1000), linewidth=1, linestyle='--', color=colors.sun_flower, label=label_V_y_comparison)
        # ax_V_y.plot(settings_graphics.y, V_y_comparison_2 / (self.hbar * 2.0 * pi * 1000), linewidth=1, linestyle='--', color=colors.sun_flower)
        # -----------------------------------------------------------------------------------------

        ax_V_y.set_xlim(settings_graphics.y_min, settings_graphics.y_max)
        ax_V_y.set_ylim(settings_graphics.potential_min, settings_graphics.potential_max)
        
        ax_V_y.set_ylabel(settings_graphics.ylabel_V_x_y_z)
        
        # ax_V_y.legend(loc='upper right', bbox_to_anchor=(1.025, 1.265), fancybox=False, framealpha=1.0)
        # ax_V_y.legend(loc='upper right', bbox_to_anchor=(1.0, 1.2), fancybox=settings_graphics.fancybox, framealpha=settings_graphics.framealpha)


    def update(self, density_y_effective, V_y):
        
        
        scaling_V = (self.hbar * 2 * pi * 1000)
        
        V_y = V_y / scaling_V
        
        
        self.line_density_y_eff.set_ydata(density_y_effective / 1e06)

        self.line_V_y.set_ydata(V_y)
