from numpy import zeros_like

import numpy as np

from .. style import colors


class fig_real_part_x(object):
    
    
    def __init__(self, ax, settings_graphics):

        self.hbar = settings_graphics.hbar
        self.m_atom = settings_graphics.m_atom
    
        self.line_real_part_x, = ax.plot(settings_graphics.x, zeros_like(settings_graphics.x), linewidth=1, linestyle='-', color=colors.wet_asphalt)
        self.line_imag_part_x, = ax.plot(settings_graphics.x, zeros_like(settings_graphics.x), linewidth=1, linestyle='-', color=colors.peter_river)

        # ax.set_ylim(-1.1, +1.1)
        ax.set_ylim(settings_graphics.real_part_min, settings_graphics.real_part_max)
        
        ax.set_xlabel(settings_graphics.xlabel_density_x)
        
        ax.grid(b=True, which='major', color=settings_graphics.color_gridlines_major, linestyle='-', linewidth=0.5)
        
        ax.set_xticks(settings_graphics.x_ticks)
        
        ax.grid(b=True, which='minor', color=settings_graphics.color_gridlines_minor, linestyle='-', linewidth=0.5, alpha=0.2)
        
        # ax.set_ylabel(r'$\Re\, \Im \;\, \mathrm{in} \;\, \mathrm{a. u.}$')
        # ax.set_ylabel(r'$\Re\, \Im \;\, \mathrm{part} \;\, \mathrm{in} \;\, \mathrm{a. u.}$')
        ax.set_ylabel(r'arbitrary units')
        
        
        
        ax_V_x = ax.twinx()
    
        self.line_V_x, = ax_V_x.plot(settings_graphics.x, zeros_like(settings_graphics.x), linewidth=1, linestyle='-', color=colors.sun_flower)

        ax_V_x.set_xlim(settings_graphics.x_min, settings_graphics.x_max)
        ax_V_x.set_ylim(settings_graphics.V_min, settings_graphics.V_max)
        
        ax_V_x.set_ylabel(settings_graphics.ylabel_V_x_y_z)
        
        # ax_V_x.legend(loc='upper right', bbox_to_anchor=(1.025, 1.25), fancybox=False, framealpha=1.0)
        # ax_V_x.legend(loc='upper right', bbox_to_anchor=(1.0, 1.2), fancybox=settings_graphics.fancybox, framealpha=settings_graphics.framealpha)
    
    
    
    def update(self, real_part_x, imag_part_x, V_x):
        
        scaling_V = self.hbar * 2 * np.pi * 1000
        
        V_x = V_x / scaling_V

        self.line_real_part_x.set_ydata(real_part_x)
        self.line_imag_part_x.set_ydata(imag_part_x)
        
        self.line_V_x.set_ydata(V_x)
