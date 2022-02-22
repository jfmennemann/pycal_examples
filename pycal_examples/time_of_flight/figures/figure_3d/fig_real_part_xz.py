import numpy as np


class fig_real_part_xz(object):

    def __init__(self, ax, settings):

        real_part_xz = np.zeros((settings.Jx, settings.Jz))

        ax.set_xlabel(settings.xlabel_xz)
        ax.set_ylabel(settings.ylabel_xz)

        ax.set_xticks(settings.z_ticks)
        ax.set_yticks(settings.x_ticks)

        self.image_real_part_xz = ax.imshow(real_part_xz,
                                            extent=[settings.z_min, settings.z_max, settings.x_min, settings.x_max],
                                            cmap='RdBu',
                                            aspect='auto',
                                            interpolation='bilinear',
                                            vmin=-1,
                                            vmax=+1,
                                            origin='lower')

    # def update(self, delta_phi_xz_selected, density_xz_selected):
    #
    #     activation_function = lambda u: 1 / (1 + np.exp(-50 * (u - 0.15)))
    #
    #     filter_density_xz_selected = activation_function(density_xz_selected / np.max(density_xz_selected))
    #
    #     delta_phi_xz_selected = filter_density_xz_selected * delta_phi_xz_selected
    #
    #     image_delta_phi_xz = np.squeeze(delta_phi_xz_selected)
    #
    #     self.image_delta_phi_xz.set_data(image_delta_phi_xz)

    def update(self, real_part_xz, density_xz):

        s = np.max(np.sqrt(density_xz))

        real_part_xz = real_part_xz / s

        self.image_real_part_xz.set_data(real_part_xz)