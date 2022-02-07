from numpy import zeros_like

from numpy import pi

from .. style import colors


class fig_density_y(object):
    
    
    def __init__(self, ax, settings):

        self.hbar = settings.hbar
        self.m_atom = settings.m_atom

        self.line_density_y, = ax.plot(settings.y, zeros_like(settings.y), linewidth=1, linestyle='-', color=colors.wet_asphalt)
        # self.line_density_y_x1, = ax.plot(settings.y, zeros_like(settings.y), linewidth=1, linestyle='-', color=colors.wet_asphalt)
        
        ax.set_xlim(settings.y_min, settings.y_max)
        
        ax.set_ylim(0, settings.density_max)
        
        ax.set_xlabel(settings.xlabel_density_y)
        
        ax.set_xticks(settings.y_ticks)
        
        ax.grid(b=True, which='major', color=settings.color_gridlines_major, linestyle='-', linewidth=0.5)
        
        ax.set_ylabel(settings.ylabel_density_x_y_z)
        
        ax.set_anchor('W')
        
        
        
        ax_V_y = ax.twinx()
        
        self.line_V_y, = ax_V_y.plot(settings.y, zeros_like(settings.y), linewidth=1, linestyle='-', color=colors.sun_flower)
        
        
        
        # -----------------------------------------------------------------------------------------
        # V_y_comparison_1 = 0.5 * self.m_atom * settings.vec_omega_y_comparison[0]**2 * (settings.y*1e-6)**2
        # V_y_comparison_2 = 0.5 * self.m_atom * settings.vec_omega_y_comparison[1]**2 * (settings.y*1e-6)**2
        #
        # label_V_y_comparison = '$({0:1.1f}, {1:1.1f})\,$kHz'.format(settings.vec_nu_y_comparison[0]/1e3, settings.vec_nu_y_comparison[1]/1e3)
        #
        # ax_V_y.plot(settings.y, V_y_comparison_1 / (self.hbar * 2.0 * pi * 1000), linewidth=1, linestyle='--', color=colors.sun_flower, label=label_V_y_comparison)
        # ax_V_y.plot(settings.y, V_y_comparison_2 / (self.hbar * 2.0 * pi * 1000), linewidth=1, linestyle='--', color=colors.sun_flower)
        # -----------------------------------------------------------------------------------------

        ax_V_y.set_xlim(settings.y_min, settings.y_max)
        ax_V_y.set_ylim(settings.potential_min, settings.potential_max)
        
        ax_V_y.set_ylabel(settings.ylabel_V_x_y_z)
        
        # ax_V_y.legend(loc='upper right', bbox_to_anchor=(1.025, 1.265), fancybox=False, framealpha=1.0)
        # ax_V_y.legend(loc='upper right', bbox_to_anchor=(1.0, 1.2), fancybox=settings.fancybox, framealpha=settings.framealpha)


    def update(self, density_y, V_y):
        
        
        scaling_V = (self.hbar * 2 * pi * 1000)
        
        V_y = V_y / scaling_V
        
        
        self.line_density_y.set_ydata(density_y)

        self.line_V_y.set_ydata(V_y)
