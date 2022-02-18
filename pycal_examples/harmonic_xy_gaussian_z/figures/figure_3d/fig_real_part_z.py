from numpy import zeros_like

from numpy import pi

from .. style import colors


class fig_real_part_z(object):

    def __init__(self, ax, settings_graphics):

        self.hbar = settings_graphics.hbar

        self.m_atom = settings_graphics.m_atom

        self.line_real_part_z, = ax.plot(settings_graphics.z, zeros_like(settings_graphics.z), linewidth=1, linestyle='-', color=colors.wet_asphalt, label=r'$\Re\, \psi$')
        self.line_imag_part_z, = ax.plot(settings_graphics.z, zeros_like(settings_graphics.z), linewidth=1, linestyle='-', color=colors.peter_river, label=r'$\Im\, \psi$')

        ax.set_xlim(settings_graphics.z_min, settings_graphics.z_max)
        
        # ax.set_ylim(-1, +1)
        ax.set_ylim(settings_graphics.real_part_min, settings_graphics.real_part_max)

        ax.set_xlabel(settings_graphics.xlabel_density_z)
        
        ax.set_xticks(settings_graphics.z_ticks)
        
        ax.grid(b=True, which='major', color=settings_graphics.color_gridlines_major, linestyle='-', linewidth=0.5)
        
        # ax.set_ylabel(r'$\mathrm{real} \;\, \mathrm{part} \;\, \mathrm{in} \;\, \mathrm{a. u.}$')
        ax.set_ylabel(r'arbitrary units')

        ax.legend(loc='lower right', bbox_to_anchor=(1.0, 0.0), fancybox=settings_graphics.fancybox, framealpha=settings_graphics.framealpha, ncol=1)
        
        
        ax_V_z = ax.twinx()
    
        self.line_V_z, = ax_V_z.plot(settings_graphics.z, zeros_like(settings_graphics.z), linewidth=1, linestyle='-', color=colors.sun_flower)


        ax_V_z.set_xlim(settings_graphics.z_min, settings_graphics.z_max)
        ax_V_z.set_ylim(settings_graphics.V_min, settings_graphics.V_max)
        
        ax_V_z.set_ylabel(settings_graphics.ylabel_V_x_y_z)
        
        # ax_V_z.legend(loc='upper right', bbox_to_anchor=(1.0, 1.2), fancybox=settings_graphics.fancybox, framealpha=settings_graphics.framealpha, ncol=2)
    
    
    def update(self, real_part_z, imag_part_z, V_z):
        
        scaling_V = self.hbar * 2 * pi * 1000
        
        V_z = V_z / scaling_V
        
        self.line_real_part_z.set_ydata(real_part_z)
        self.line_imag_part_z.set_ydata(imag_part_z)

        self.line_V_z.set_ydata(V_z)
