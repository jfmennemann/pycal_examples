from numpy import zeros_like

import numpy as np

from ..style import colors


class fig_phase_x(object):
    
    
    def __init__(self, ax, settings_graphics):

        self.hbar = settings_graphics.hbar
        self.m_atom = settings_graphics.m_atom
    
        self.line_phase_x, = ax.plot(settings_graphics.x, zeros_like(settings_graphics.x), linewidth=1, linestyle='-', color=colors.wet_asphalt)
        

        ax.set_ylim(-1.1, +1.1)
        # ax.set_ylim(settings_graphics.real_part_min, settings_graphics.real_part_max)
        
        ax.set_xlabel(settings_graphics.xlabel_density_x)
        
        ax.grid(b=True, which='major', color=settings_graphics.color_gridlines_major, linestyle='-', linewidth=0.5)
        
        ax.set_xticks(settings_graphics.x_ticks)
        
        ax.grid(b=True, which='minor', color=settings_graphics.color_gridlines_minor, linestyle='-', linewidth=0.5, alpha=0.2)

        ax.set_ylabel(r'$\cos(\varphi)$')
        
        
        ax_V_x = ax.twinx()
    
        self.line_V_x, = ax_V_x.plot(settings_graphics.x, zeros_like(settings_graphics.x), linewidth=1, linestyle='-', color=colors.sun_flower)

        ax_V_x.set_xlim(settings_graphics.x_min, settings_graphics.x_max)
        ax_V_x.set_ylim(settings_graphics.potential_min, settings_graphics.potential_max)
        
        ax_V_x.set_ylabel(settings_graphics.ylabel_V_x_y_z)
        
        # ax_V_x.legend(loc='upper right', bbox_to_anchor=(1.025, 1.25), fancybox=False, framealpha=1.0)
        # ax_V_x.legend(loc='upper right', bbox_to_anchor=(1.0, 1.2), fancybox=settings_graphics.fancybox, framealpha=settings_graphics.framealpha)
    
    
    
    def update(self, phase_x, V_x):
        
        scaling_V = (self.hbar * 2 * np.pi * 1000)
        
        V_x = V_x / scaling_V

        self.line_phase_x.set_ydata(np.cos(phase_x))
        
        self.line_V_x.set_ydata(V_x)
