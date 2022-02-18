import numpy as np


class fig_delta_phi_xy(object):

    def __init__(self, ax, settings):

        delta_phi_xy = np.zeros((settings.Jx, settings.Jy))

        ax.set_xlabel(settings.xlabel_xy)
        ax.set_ylabel(settings.ylabel_xy)

        ax.set_xticks(settings.y_ticks)
        ax.set_yticks(settings.x_ticks)

        ax.set_anchor('W')

        delta_phi_xy_min = -0.5 * np.pi
        delta_phi_xy_max = +0.5 * np.pi

        self.image_delta_phi_xy = ax.imshow(delta_phi_xy,
                                            extent=[settings.y_min, settings.y_max, settings.x_min, settings.x_max],
                                            cmap='RdBu',
                                            aspect='equal',
                                            interpolation='bilinear',
                                            vmin=delta_phi_xy_min,
                                            vmax=delta_phi_xy_max,
                                            origin='lower')

    def update(self, delta_phi_xy_selected, density_xy_selected):

        activation_function = lambda u: 1 / (1 + np.exp(-50 * (u - 0.15)))

        filter_density_xy_selected = activation_function(density_xy_selected / np.max(density_xy_selected))

        delta_phi_xy_selected = filter_density_xy_selected * delta_phi_xy_selected

        image_delta_phi_xy = np.squeeze(delta_phi_xy_selected)

        self.image_delta_phi_xy.set_data(image_delta_phi_xy)
