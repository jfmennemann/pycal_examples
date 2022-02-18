import numpy as np

from numpy import zeros_like

from numpy import pi

from ..style import colors


class fig_phase_y(object):
    
    
    def __init__(self, ax, settings_graphics):

        self.hbar = settings_graphics.hbar
        self.m_atom = settings_graphics.m_atom

        self.line_phase_y, = ax.plot(settings_graphics.y, zeros_like(settings_graphics.y), linewidth=1, linestyle='-', color=colors.wet_asphalt)

        
        ax.set_xlim(settings_graphics.y_min, settings_graphics.y_max)
        
        ax.set_ylim(-1.1, +1.1)
        # ax.set_ylim(settings_graphics.real_part_min, settings_graphics.real_part_max)

        ax.set_xlabel(settings_graphics.xlabel_density_y)
        
        ax.set_xticks(settings_graphics.y_ticks)
        
        ax.grid(b=True, which='major', color=settings_graphics.color_gridlines_major, linestyle='-', linewidth=0.5)
        
        # ax.set_ylabel(r'$\mathrm{real} \;\, \mathrm{part} \;\, \mathrm{in} \;\, \mathrm{a. u.}$')
        ax.set_ylabel(r'$\cos(\varphi)$')

        ax.set_anchor('W')
        
        
        
        ax_V_y = ax.twinx()
        
        self.line_V_y, = ax_V_y.plot(settings_graphics.y, zeros_like(settings_graphics.y), linewidth=1, linestyle='-', color=colors.sun_flower)
        


        ax_V_y.set_xlim(settings_graphics.y_min, settings_graphics.y_max)
        ax_V_y.set_ylim(settings_graphics.potential_min, settings_graphics.potential_max)
        
        ax_V_y.set_ylabel(settings_graphics.ylabel_V_x_y_z)
        
        # ax_V_y.legend(loc='upper right', bbox_to_anchor=(1.025, 1.265), fancybox=False, framealpha=1.0)
        # ax_V_y.legend(loc='upper right', bbox_to_anchor=(1.0, 1.2), fancybox=settings_graphics.fancybox, framealpha=settings_graphics.framealpha)


    def update(self, phase_y, V_y):

        scaling_V = (self.hbar * 2 * pi * 1000)
        
        V_y = V_y / scaling_V

        self.line_phase_y.set_ydata(np.cos(phase_y))

        self.line_V_y.set_ydata(V_y)
