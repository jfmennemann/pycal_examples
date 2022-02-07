from numpy import zeros_like

import numpy as np

from ..style import colors


class fig_density_x_eff(object):
    
    
    def __init__(self, ax, settings_graphics):

        self.hbar = settings_graphics.hbar
        self.m_atom = settings_graphics.m_atom
    
        self.line_density_x_eff, = ax.plot(settings_graphics.x, zeros_like(settings_graphics.x), linewidth=1, linestyle='-', color=colors.wet_asphalt)

        ax.set_xlim(settings_graphics.z_min, settings_graphics.z_max)
        ax.set_ylim(0, settings_graphics.density_x_effective_max)
        
        ax.set_xlabel(settings_graphics.xlabel_density_x)
        
        ax.grid(b=True, which='major', color=settings_graphics.color_gridlines_major, linestyle='-', linewidth=0.5)
        
        ax.set_xticks(settings_graphics.x_ticks)
        
        ax.grid(b=True, which='minor', color=settings_graphics.color_gridlines_minor, linestyle='-', linewidth=0.5, alpha=0.2)
        
        ax.set_ylabel(settings_graphics.ylabel_density_x_y_z_effective)
        
        
        
        ax_V_x = ax.twinx()
    
        self.line_V_x, = ax_V_x.plot(settings_graphics.x, zeros_like(settings_graphics.x), linewidth=1, linestyle='-', color=colors.sun_flower)

        # -----------------------------------------------------------------------------------------
        # V_x_comparison_1 = 0.5 * settings_graphics.m_atom * settings_graphics.vec_omega_x_comparison[0]**2 * (settings_graphics.x*1e-6)**2
        # V_x_comparison_2 = 0.5 * settings_graphics.m_atom * settings_graphics.vec_omega_x_comparison[1]**2 * (settings_graphics.x*1e-6)**2
        #
        # label_V_x_comparison = '$({0:1.1f}, {1:1.1f})\,$kHz'.format(settings_graphics.vec_nu_x_comparison[0]/1e3, settings_graphics.vec_nu_x_comparison[1]/1e3)
        #
        # ax_V_x.plot(settings_graphics.x, V_x_comparison_1 / (self.hbar * 2.0 * np.pi * 1000), linewidth=1, linestyle='--', color=colors.sun_flower, label=label_V_x_comparison)
        # ax_V_x.plot(settings_graphics.x, V_x_comparison_2 / (self.hbar * 2.0 * np.pi * 1000), linewidth=1, linestyle='--', color=colors.sun_flower)
        # -----------------------------------------------------------------------------------------
        
        ax_V_x.set_xlim(settings_graphics.x_min, settings_graphics.x_max)
        ax_V_x.set_ylim(settings_graphics.potential_min, settings_graphics.potential_max)
        
        ax_V_x.set_ylabel(settings_graphics.ylabel_V_x_y_z)
        
        # ax_V_x.legend(loc='upper right', bbox_to_anchor=(1.025, 1.25), fancybox=False, framealpha=1.0)
        
        # ax_V_x.legend(loc='upper right', bbox_to_anchor=(1.0, 1.2), fancybox=settings_graphics.fancybox, framealpha=settings_graphics.framealpha)
    
    
    
    def update(self, density_x, V_x):
        
        scaling_V = (self.hbar * 2 * np.pi * 1000)
        
        V_x = V_x / scaling_V

        self.line_density_x_eff.set_ydata(density_x / 1e06)
        
        self.line_V_x.set_ydata(V_x)
