import numpy as np


class fig_density_xz(object):

    def __init__(self, ax, settings):
    
        Jx = settings.Jx
        Jz = settings.Jz
        
        density_xz = np.zeros((Jx, Jz))
    
        ax.set_xlabel(settings.xlabel_xz)
        ax.set_ylabel(settings.ylabel_xz)
        
        ax.set_xticks(settings.z_ticks)
        ax.set_yticks(settings.x_ticks)
        
        self.image_density_xz = ax.imshow(
                                        density_xz,
                                        extent=[settings.z_min, settings.z_max, settings.x_min, settings.x_max],
                                        cmap=settings.colormap_density,
                                        aspect='auto',
                                        interpolation='bilinear',
                                        vmin=0,
                                        vmax=1,
                                        origin='lower')
    
       
    def update(self, density_xz):

        image_density_xz = np.squeeze(density_xz) / np.max(density_xz)
        
        self.image_density_xz.set_data(image_density_xz)
