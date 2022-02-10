import numpy as np


class fig_density_xy(object):
    
    
    def __init__(self, ax, settings):
        
        Jx = settings.Jx
        Jy = settings.Jy
        
        density_xy = np.zeros((Jx, Jy))
        
        
        ax.set_xlabel(settings.xlabel_xy)
        ax.set_ylabel(settings.ylabel_xy)
        
        ax.set_xticks(settings.y_ticks)
        ax.set_yticks(settings.x_ticks)
        
        ax.set_anchor('W')
        
        self.image_density_xy = ax.imshow(density_xy,
                                          extent=[settings.y_min, settings.y_max, settings.x_min, settings.x_max],
                                          cmap=settings.colormap_density,
                                          aspect='equal',
                                          interpolation='bilinear',
                                          vmin=0,
                                          vmax=1,
                                          origin='lower')

        self.density_max = settings.density_max


    def update(self, density_xy):

        # if np.max(density_xy) > 0:
        #
        #     image_density_xy = density_xy / np.max(density_xy)
        #
        # else:
        #
        #     image_density_xy = density_xy

        image_density_xy = density_xy / self.density_max

        self.image_density_xy.set_data(image_density_xy)
