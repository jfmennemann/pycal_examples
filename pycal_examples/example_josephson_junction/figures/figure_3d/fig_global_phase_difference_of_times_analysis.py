import numpy as np

from .. style import colors


class fig_global_phase_difference_of_times_analysis(object):

    def __init__(self, ax, settings):

        self.line_global_phase_difference_of_times_analysis, = ax.plot([], [], linewidth=1.0, linestyle='-', color=colors.wet_asphalt, label=r'')

        ax.set_xlim(settings.t_min, settings.t_max)
        ax.set_ylim(-0.55, 0.55)

        ax.set_xlabel(settings.xlabel_phase_difference_of_times_analysis)
        ax.set_ylabel(r'$\Delta \, \Phi\; / \; \pi$')

        ax.set_xticks(settings.t_ticks_major, minor=False)
        ax.set_xticks(settings.t_ticks_minor, minor=True)

        ax.set_yticks([-0.5, 0, 0.5], minor=False)
        ax.set_yticks([-0.25, 0.25], minor=True)

        ax.grid(b=True, which='major', color=colors.color_gridlines_major, linestyle='-', linewidth=0.5)
        ax.grid(b=True, which='minor', color=colors.color_gridlines_minor, linestyle='-', linewidth=0.5, alpha=0.2)
        
        # ax.legend(loc='upper right', bbox_to_anchor=(1.05, 1.25), ncol=3, fancybox=settings.fancybox, framealpha=settings.framealpha)


    def update(self,
               global_phase_difference_of_times_analysis,
               times_analysis,
               nr_times_analysis):

        self.line_global_phase_difference_of_times_analysis.set_ydata(global_phase_difference_of_times_analysis[0:nr_times_analysis] / np.pi)
        self.line_global_phase_difference_of_times_analysis.set_xdata(times_analysis[0:nr_times_analysis] / 1e-3)
