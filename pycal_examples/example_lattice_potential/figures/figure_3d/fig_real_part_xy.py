import numpy as np

# from ..style import colors


class fig_real_part_xy(object):
    
    
    def __init__(self, ax, settings):
        
        Jx = settings.Jx
        Jy = settings.Jy
        
        real_part_xy = np.zeros((Jx, Jy))
        
        
        ax.set_xlabel(settings.xlabel_xy)
        ax.set_ylabel(settings.ylabel_xy)
        
        ax.set_xticks(settings.y_ticks)
        ax.set_yticks(settings.x_ticks)
        
        ax.set_anchor('W')
        
        self.image_real_part_xy = ax.imshow(real_part_xy,
                                            extent=[settings.y_min, settings.y_max, settings.x_min, settings.x_max],
                                            cmap='RdBu',
                                            aspect='equal',
                                            interpolation='bilinear',
                                            vmin=-1,
                                            vmax=1,
                                            origin='lower')


    def update(self, real_part_xy, density_xy):

        s = np.max(np.sqrt(density_xy))

        real_part_xy = real_part_xy / s

        self.image_real_part_xy.set_data(real_part_xy)
