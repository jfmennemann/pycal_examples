import numpy as np


class fig_tof_xy(object):

    def __init__(self, ax, settings):
        
        ax.set_xlabel(settings.xlabel_xy)
        ax.set_ylabel(settings.ylabel_xy)
        
        # ax.set_xticks([y_min, 0, y_max])
        # ax.set_yticks([x_min, 0, x_max])
        
        ax.set_xticks(settings.x_ticks)
        ax.set_yticks(settings.y_ticks)
        
        ax.set_anchor('W')
        
        image_tof_xy = np.zeros((settings.Jx, settings.Jy))
        
        self.im_tof_xy = ax.imshow(np.transpose(image_tof_xy),
                                   extent=[settings.x_min, settings.x_max, settings.y_min, settings.y_max],
                                   cmap=settings.colormap_density,
                                   aspect='equal',
                                   interpolation='bilinear',
                                   vmin=0,
                                   vmax=1,
                                   origin='lower')
        
        
        self.ax = ax
        
        
    def update(self, density_tof_xy, x_tof, y_tof):
        
        x_tof = x_tof / 1e-6
        y_tof = y_tof / 1e-6
        
        dx_tof = x_tof[1] - x_tof[0]
        dy_tof = y_tof[1] - y_tof[0]
        
        x_tof_min = x_tof[0]
        y_tof_min = y_tof[0]
        
        x_tof_max = x_tof[-1] + dx_tof
        y_tof_max = y_tof[-1] + dy_tof
        
        
        density_tof_xy = density_tof_xy / np.max(density_tof_xy)
        
        
        self.im_tof_xy.set_data(np.transpose(density_tof_xy))
        self.im_tof_xy.set_extent([x_tof_min, x_tof_max, y_tof_min, y_tof_max])
        
        self.ax.set_xticks([x_tof_min, 0, x_tof_max])
        self.ax.set_yticks([y_tof_min, 0, y_tof_max])
