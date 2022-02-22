from .. style import colors


class fig_number_imbalance_of_times_analysis(object):

    def __init__(self, ax, settings):

        self.line_number_imbalance_of_times_analysis, = ax.plot([], [], linewidth=1.0, linestyle='-', color=colors.wet_asphalt, label='single run')

        # self.line_number_imbalance_mean_of_times_analysis,       = ax.plot([], [], linewidth=1.0, linestyle='-', color=colors.wet_asphalt, label='mean')

        # self.line_number_imbalance_std_top_of_times_analysis,    = ax.plot([], [], linewidth=0.5, linestyle='--', color=colors.wet_asphalt, label='mean $\pm$ std')

        # self.line_number_imbalance_std_top_of_times_analysis,    = ax.plot([], [], linewidth=0.5, linestyle='--', color=colors.wet_asphalt)
        # self.line_number_imbalance_std_bottom_of_times_analysis, = ax.plot([], [], linewidth=0.5, linestyle='--', color=colors.wet_asphalt)


        ax.set_xlim(settings.t_min, settings.t_max)
        ax.set_ylim(-0.35, 0.35)
        
        ax.set_xlabel(settings.xlabel_number_imbalance_of_times_analysis)
        ax.set_ylabel(settings.ylabel_number_imbalance_of_times_analysis)
        
        ax.set_xticks(settings.t_ticks_major, minor=False)
        ax.set_xticks(settings.t_ticks_minor, minor=True)
        
        ax.set_yticks([-0.2, 0.0, 0.2], minor=False)
        ax.set_yticks([-0.3, -0.1, 0.1, 0.3], minor=True)
        
        ax.grid(b=True, which='major', color=colors.color_gridlines_major, linestyle='-', linewidth=0.5)
        ax.grid(b=True, which='minor', color=colors.color_gridlines_minor, linestyle='-', linewidth=0.5, alpha=0.2)
        
        # ax.legend(loc='upper right', bbox_to_anchor=(1.05, 1.2), ncol=3, fancybox=settings.fancybox, framealpha=settings.framealpha)
        
        
    def update(self,
               # number_imbalance_mean_of_times_analysis,
               # number_imbalance_std_of_times_analysis,
               number_imbalance_of_times_analysis,
               times_analysis,
               nr_times_analysis):
        
        # number_imbalance_mean_temp = number_imbalance_mean_of_times_analysis[0:nr_times_analysis]
        # number_imbalance_std_temp = number_imbalance_std_of_times_analysis[0:nr_times_analysis]

        # self.line_number_imbalance_mean_of_times_analysis.set_ydata(number_imbalance_mean_temp)
        # self.line_number_imbalance_mean_of_times_analysis.set_xdata(times_analysis_temp)
        
        # self.line_number_imbalance_std_top_of_times_analysis.set_ydata(number_imbalance_mean_temp + number_imbalance_std_temp)
        # self.line_number_imbalance_std_top_of_times_analysis.set_xdata(times_analysis_temp)
    
        # self.line_number_imbalance_std_bottom_of_times_analysis.set_ydata(number_imbalance_mean_temp - number_imbalance_std_temp)
        # self.line_number_imbalance_std_bottom_of_times_analysis.set_xdata(times_analysis_temp)

        self.line_number_imbalance_of_times_analysis.set_ydata(number_imbalance_of_times_analysis[0:nr_times_analysis])
        self.line_number_imbalance_of_times_analysis.set_xdata(times_analysis[0:nr_times_analysis]/1e-3)
