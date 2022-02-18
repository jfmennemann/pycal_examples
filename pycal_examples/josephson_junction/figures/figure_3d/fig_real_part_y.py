from numpy import zeros_like

from numpy import pi

from .. style import colors


class fig_real_part_y(object):
    
    
    def __init__(self, ax, settings):

        self.hbar = settings.hbar
        self.m_atom = settings.m_atom

        self.line_real_part_y, = ax.plot(settings.y, zeros_like(settings.y), linewidth=1, linestyle='-', color=colors.wet_asphalt)
        self.line_imag_part_y, = ax.plot(settings.y, zeros_like(settings.y), linewidth=1, linestyle='-', color=colors.peter_river)
        
        ax.set_xlim(settings.y_min, settings.y_max)
        
        # ax.set_ylim(-1.1, +1.1)
        ax.set_ylim(settings.real_part_min, settings.real_part_max)

        ax.set_xlabel(settings.xlabel_density_y)
        
        ax.set_xticks(settings.y_ticks)
        
        ax.grid(b=True, which='major', color=settings.color_gridlines_major, linestyle='-', linewidth=0.5)
        
        # ax.set_ylabel(r'$\mathrm{real} \;\, \mathrm{part} \;\, \mathrm{in} \;\, \mathrm{a. u.}$')
        ax.set_ylabel(r'arbitrary units')

        ax.set_anchor('W')
        
        
        
        ax_V_y = ax.twinx()
        
        self.line_V_y, = ax_V_y.plot(settings.y, zeros_like(settings.y), linewidth=1, linestyle='-', color=colors.sun_flower)
        


        ax_V_y.set_xlim(settings.y_min, settings.y_max)
        ax_V_y.set_ylim(settings.potential_min, settings.potential_max)
        
        ax_V_y.set_ylabel(settings.ylabel_V_x_y_z)
        
        # ax_V_y.legend(loc='upper right', bbox_to_anchor=(1.025, 1.265), fancybox=False, framealpha=1.0)
        # ax_V_y.legend(loc='upper right', bbox_to_anchor=(1.0, 1.2), fancybox=settings.fancybox, framealpha=settings.framealpha)


    def update(self, real_part_y, imag_part_y, V_y):
        
        
        scaling_V = (self.hbar * 2 * pi * 1000)
        
        V_y = V_y / scaling_V
        
        
        self.line_real_part_y.set_ydata(real_part_y)
        self.line_imag_part_y.set_ydata(imag_part_y)

        self.line_V_y.set_ydata(V_y)
