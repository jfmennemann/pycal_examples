import numpy as np


class fig_tof_xz(object):

    def __init__(self, ax, settings):

        ax.set_xlabel(settings.xlabel_xz)
        ax.set_ylabel(settings.ylabel_xz)
        
        ax.set_xticks(settings.x_ticks)
        ax.set_yticks(settings.z_ticks)
        
        image_tof_xz = np.zeros((settings.Jx, settings.Jz))
        
        self.im_tof_xz = ax.imshow(np.transpose(image_tof_xz),
                                   extent=[settings.x_min, settings.x_max, settings.z_min, settings.z_max],
                                   cmap=settings.colormap_density,
                                   aspect='auto',
                                   interpolation='bilinear',
                                   vmin=0,
                                   vmax=1,
                                   origin='lower')
    
        self.ax = ax


    def update(self, density_tof_xz, x_tof, z_tof):
        
        x_tof = x_tof / 1e-6
        z_tof = z_tof / 1e-6
        
        dx_tof = x_tof[1] - x_tof[0]
        dz_tof = z_tof[1] - z_tof[0]
        
        x_tof_min = x_tof[0]
        z_tof_min = z_tof[0]
        
        x_tof_max = x_tof[-1] + dx_tof
        z_tof_max = z_tof[-1] + dz_tof
        
        density_tof_xz = density_tof_xz / np.max(density_tof_xz)
        
        self.im_tof_xz.set_data(np.transpose(density_tof_xz))
        self.im_tof_xz.set_extent([x_tof_min, x_tof_max, z_tof_min, z_tof_max])
        
        self.ax.set_xticks([x_tof_min, 0, x_tof_max])
        self.ax.set_yticks([z_tof_min, 0, z_tof_max])
