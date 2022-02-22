import numpy as np


def activation_function(u):

    return 1.0 / (1.0 + np.exp(-50.0 * (u - 0.15)))


class fig_phase_difference_xz(object):

    def __init__(self, ax, settings):

        Jx = settings.Jx
        Jz = settings.Jz

        phase_difference_xz = np.zeros((Jx, Jz))

        ax.set_xlabel(settings.xlabel_xz)
        ax.set_ylabel(settings.ylabel_xz)

        ax.set_xticks(settings.z_ticks)
        ax.set_yticks(settings.x_ticks)

        phase_difference_xz_min = -0.5*np.pi
        phase_difference_xz_max = +0.5*np.pi

        self.image_phase_difference_xz = ax.imshow(
                                                phase_difference_xz,
                                                extent=[settings.z_min, settings.z_max, settings.x_min, settings.x_max],
                                                cmap='RdBu',
                                                aspect='auto',
                                                interpolation='bilinear',
                                                vmin=phase_difference_xz_min,
                                                vmax=phase_difference_xz_max,
                                                origin='lower')

    def update(self, phase_difference_xz, density_xz):

        filter_density_xz = activation_function(density_xz / np.max(density_xz))

        image_phase_difference_xz = filter_density_xz * phase_difference_xz

        self.image_phase_difference_xz.set_data(image_phase_difference_xz)
