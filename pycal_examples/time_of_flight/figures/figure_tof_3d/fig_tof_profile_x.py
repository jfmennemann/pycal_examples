import numpy as np

from numpy import zeros_like

from ..style import colors


class fig_tof_profile_x(object):

    def __init__(self, ax, settings):
    
        self.line_profile_tof_x, = ax.plot(settings.x, zeros_like(settings.x), linewidth=1, linestyle='-', color=colors.wet_asphalt)
        
        ax.set_xlim(settings.x_min, settings.x_max)
        ax.set_ylim(0, 1.5)
        
        ax.set_xlabel(settings.xlabel_profile_tof_x)
        ax.set_ylabel(settings.ylabel_profile_tof_x)
        
        ax.grid(b=True, which='major', color=settings.color_gridlines_major, linestyle='-', linewidth=0.5)
        
        ax.set_xticks(settings.x_ticks)
        
        self.ax = ax
    
    
    def update(self, profile_tof_x, x_tof):
        
        x_tof = x_tof / 1e-6
        
        dx_tof = x_tof[1] - x_tof[0]
        
        x_tof_min = x_tof[0]
        x_tof_max = x_tof[-1] + dx_tof
        
        profile_tof_x = profile_tof_x / np.max(profile_tof_x)
        
        self.line_profile_tof_x.set_xdata(x_tof)
        self.line_profile_tof_x.set_ydata(profile_tof_x)
        
        self.ax.set_xlim(x_tof_min, x_tof_max)
        
        self.ax.set_xticks(np.linspace(x_tof_min, x_tof_max, 5, True))
