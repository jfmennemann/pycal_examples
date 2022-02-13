import matplotlib.pyplot as plt

from PyQt5 import QtWidgets 


import numpy as np



from .fig_density_xz import fig_density_xz
from .fig_density_xy import fig_density_xy

from .fig_density_x import fig_density_x
from .fig_density_y import fig_density_y
from .fig_density_z import fig_density_z



from .fig_phase_difference_xz import fig_phase_difference_xz

from .fig_phase_difference_z_x1_x2 import fig_phase_difference_z_x1_x2

from .fig_control_inputs_of_times import fig_control_inputs_of_times

from .fig_global_phase_difference_of_times_analysis import fig_global_phase_difference_of_times_analysis

from .fig_number_imbalance_of_times_analysis import fig_number_imbalance_of_times_analysis


from .. style import colors



import scipy.constants

hbar = scipy.constants.hbar



class Figure3d(object):


    def __init__(self,
                 x,
                 y,
                 z,
                 times,
                 settings_figure_3d):


        m_atom = settings_figure_3d['m_atom']

        density_max = settings_figure_3d['density_max']

        # density_z_eff_max = settings_figure_3d['density_z_eff_max']

        V_max = settings_figure_3d['V_max']

        abs_z_restr = settings_figure_3d['abs_z_restr']


        x = x / 1e-6
        y = y / 1e-6
        z = z / 1e-6

        times = times / 1e-3


        x_min = x[0]
        y_min = y[0]
        z_min = z[0]
                
        x_max = -x_min
        y_max = -y_min
        z_max = -z_min
        
        t_min = times[0]
        t_max = times[-1]

        Jx = x.size
        Jy = y.size
        Jz = z.size

        abs_z_restr = abs_z_restr / 1e-6

        indices_z_restr = np.abs(z) < abs_z_restr
        
        

        x_ticks = np.array([-2, -1, 0, 1, 2])
        y_ticks = np.array([-1, 0, 1])


        if np.round(z_max) == 5:
            z_ticks = np.array([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])

        elif np.round(z_max) == 10:
            z_ticks = np.array([-10, -5, 0, 5, 10])

        elif np.round(z_max) == 20:
            z_ticks = np.array([-20, -10, 0, 10, 20])

        elif np.round(z_max) == 40:
            z_ticks = np.array([-40, -30, -20, -10, 0, 10, 20, 30, 40])

        elif np.round(z_max) == 50:
            z_ticks = np.array([-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50])

        elif np.round(z_max) == 54:
            z_ticks = np.array([-54, 0, 54])

        elif np.round(z_max) == 60:
            z_ticks = np.array([-60, -40, -20, 0, 20, 40, 60])

        elif np.round(z_max) == 80:
            z_ticks = np.array([-80, -60, -40, -20, 0, 20, 40, 60, 80])

        elif np.round(z_max) == 100:
            z_ticks = np.array([-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100])

        elif np.round(z_max) == 200:
            z_ticks = np.array([-200, -160, -120, -80, -40, 0, 40, 80, 120, 160, 200])

        else:
            z_ticks = np.array([z_min, z_max])



        if t_max == 2:
            t_ticks_major = np.array([0, 1, 2])

        elif t_max == 2.5:
            t_ticks_major = np.array([0, 0.5, 1, 1.5, 2, 2.5])

        elif t_max == 4:
            t_ticks_major = np.array([0, 1, 2, 3, 4])

        elif t_max == 5:
            t_ticks_major = np.array([0, 1, 2, 3, 4, 5])

        elif t_max == 8:
            t_ticks_major = np.array([0, 2, 4, 6, 8])

        elif t_max == 10:
            t_ticks_major = np.array([0, 2, 4, 6, 8, 10])

        elif t_max == 20:
            t_ticks_major = np.array([0, 4, 8, 12, 16, 20])

        elif t_max == 40:
            t_ticks_major = np.array([0, 10, 20, 30, 40])

        elif t_max == 80:
            t_ticks_major = np.array([0, 20, 40, 60, 80])

        elif t_max == 160:
            t_ticks_major = np.array([0, 40, 80, 120, 160])

        elif t_max == 200:
            t_ticks_major = np.array([0, 40, 80, 120, 160, 200])

        else:
            t_ticks_major = np.array([0, t_max])



        t_ticks_minor = 0.5 * (t_ticks_major[0:-1] + t_ticks_major[1:])



        settings_graphics = type('', (), {})()

        settings_graphics.density_max = density_max

        # settings_graphics.density_z_eff_max = density_z_eff_max


        settings_graphics.hbar = hbar
        settings_graphics.m_atom = m_atom

        settings_graphics.abs_z_restr = abs_z_restr

        settings_graphics.indices_z_restr = indices_z_restr




        settings_graphics.ylabel_density_x_y_z = r'$\mathrm{density} \;\, \mathrm{in} \;\, \mathrm{m}^{-3}$'
        settings_graphics.ylabel_density_x_y_z_effective = r'$\mathrm{density} \;\, \mathrm{in} \;\, \mu \mathrm{m}^{-1}$'




        settings_graphics.real_part_min = -np.sqrt(settings_graphics.density_max)
        settings_graphics.real_part_max = +np.sqrt(settings_graphics.density_max)

        settings_graphics.V_max = V_max

        settings_graphics.phase_difference_of_times_analysis_min = -1.6
        settings_graphics.phase_difference_of_times_analysis_max = +1.6

        settings_graphics.number_imbalance_of_times_analysis_min = -0.4
        settings_graphics.number_imbalance_of_times_analysis_max = +0.4



        settings_graphics.times = times
        
        settings_graphics.t_min = t_min
        settings_graphics.t_max = t_max



        settings_graphics.t_ticks_major = t_ticks_major
        settings_graphics.t_ticks_minor = t_ticks_minor
        
        
        settings_graphics.Jx = Jx
        settings_graphics.Jy = Jy
        settings_graphics.Jz = Jz
        
        settings_graphics.x = x
        settings_graphics.y = y
        settings_graphics.z = z
        
        settings_graphics.x_min = x_min
        settings_graphics.y_min = y_min
        settings_graphics.z_min = z_min
        
        settings_graphics.x_max = x_max
        settings_graphics.y_max = y_max
        settings_graphics.z_max = z_max
        
        settings_graphics.x_ticks = x_ticks
        settings_graphics.y_ticks = y_ticks
        settings_graphics.z_ticks = z_ticks
        
        
        
        # settings_graphics.vec_nu_x_comparison = np.array([1.4e3, 2.1e3])
        # settings_graphics.vec_omega_x_comparison = 2 * np.pi * settings_graphics.vec_nu_x_comparison
        
        # settings_graphics.vec_nu_y_comparison = np.array([1.4e3, 2.1e3])
        # settings_graphics.vec_omega_y_comparison = 2 * np.pi * settings_graphics.vec_nu_y_comparison
        
        # settings_graphics.vec_nu_z_comparison = np.array([7, 11])
        # settings_graphics.vec_omega_z_comparison = 2 * np.pi * settings_graphics.vec_nu_z_comparison
        
        
        settings_graphics.xlabel_phase_difference_of_times_analysis = r'$t             \;\, \mathrm{in} \;\, \mathrm{ms}$'
        settings_graphics.ylabel_phase_difference_of_times_analysis = r'$\bar{\varphi} \;\, \mathrm{in} \;\, \mathrm{rad}$'
        
        settings_graphics.xlabel_number_imbalance_of_times_analysis = r'$t             \;\, \mathrm{in} \;\, \mathrm{ms}$'
        settings_graphics.ylabel_number_imbalance_of_times_analysis = r'$\Delta N$'
        
        
        settings_graphics.phase_difference_of_times_analysis_min = -1.6
        settings_graphics.phase_difference_of_times_analysis_max = +1.6
        
        settings_graphics.number_imbalance_of_times_analysis_min = -0.4
        settings_graphics.number_imbalance_of_times_analysis_max = +0.4

        settings_graphics.xlabel_profile_tof_x = r'$x \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings_graphics.ylabel_profile_tof_x = r'$\mathrm{density} \;\, \mathrm{in} \;\, \mathrm{a.u.}$'
        
        
        
        settings_graphics.xlabel_xz = r'$z \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings_graphics.ylabel_xz = r'$x \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        
        settings_graphics.xlabel_xy = r'$y \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings_graphics.ylabel_xy = r'$x \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        
        settings_graphics.ylabel_V_x_y_z = r'$V \;\, \mathrm{in} \;\, h \times \mathrm{kHz}$'

        
        settings_graphics.xlabel_density_x = r'$x \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings_graphics.xlabel_density_y = r'$y \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings_graphics.xlabel_density_z = r'$z \;\, \mathrm{in} \;\, \mu \mathrm{m}$'



        
        
        
        
        settings_graphics.colormap_density = colors.colormap_1
        
        settings_graphics.color_gridlines_major = colors.color_gridlines_major
        settings_graphics.color_gridlines_minor = colors.color_gridlines_minor

        settings_graphics.framealpha = 1.0
        settings_graphics.fancybox = False
        

        self.fig_name = "figure_main"
                
        self.fig = plt.figure(self.fig_name, facecolor="white", constrained_layout=False)
        
        
        
        window = self.fig.canvas.window()
        
        window.findChild(QtWidgets.QToolBar).setVisible(False)
        window.statusBar().setVisible(False)
        
        
        
        resolution = '2560x1440'
        # resolution = '1920x1080'
        
        if resolution == '2560x1440':
            
            n_pixels_x = 1800
            n_pixels_y = 900
            
            pos_x = 2560 - n_pixels_x
            pos_y = 0
            
            plt.rcParams.update({'font.size': 10})
            
        
        if resolution == '1920x1080':
            
            n_pixels_x = 800
            n_pixels_y = 400
            
            pos_x = 1920 - n_pixels_x
            pos_y = 0
            
            plt.rcParams.update({'font.size': 6})
        




        window.setGeometry(pos_x, pos_y, n_pixels_x, n_pixels_y)

        self.gridspec = self.fig.add_gridspec(ncols=4, nrows=5, left=0.055, bottom=0.065, right=0.985, top=0.965, wspace=0.50, hspace=0.65, width_ratios=[2.5, 1, 1, 2], height_ratios=[1, 1, 1, 1, 1])
        
        
        
        # =========================================================================================
        ax_00 = self.fig.add_subplot(self.gridspec[0, 0])
        ax_01 = self.fig.add_subplot(self.gridspec[0, 1])
        ax_03 = self.fig.add_subplot(self.gridspec[0, 3])

        ax_10 = self.fig.add_subplot(self.gridspec[1, 0])
        ax_11 = self.fig.add_subplot(self.gridspec[1, 1])
        ax_12 = self.fig.add_subplot(self.gridspec[1, 2])
        ax_13 = self.fig.add_subplot(self.gridspec[1, 3])
        
        ax_20 = self.fig.add_subplot(self.gridspec[2, 0])
        # ax_21 = self.fig.add_subplot(self.gridspec[2, 1])
        # ax_22 = self.fig.add_subplot(self.gridspec[2, 2])
        ax_23 = self.fig.add_subplot(self.gridspec[2, 3])

        ax_30 = self.fig.add_subplot(self.gridspec[3, 0])
        # ax_40 = self.fig.add_subplot(self.gridspec[4, 0])

        self.fig_density_xz = fig_density_xz(ax_00, settings_graphics)
        self.fig_density_xy = fig_density_xy(ax_01, settings_graphics)

        self.fig_density_z = fig_density_z(ax_10, settings_graphics)
        self.fig_density_y = fig_density_y(ax_11, settings_graphics)
        self.fig_density_x = fig_density_x(ax_12, settings_graphics)

        self.fig_phase_difference_xz = fig_phase_difference_xz(ax_20, settings_graphics)

        # self.fig_real_part_z = fig_real_part_z(ax_20, settings_graphics)
        # self.fig_real_part_y = fig_real_part_y(ax_21, settings_graphics)
        # self.fig_real_part_x = fig_real_part_x(ax_22, settings_graphics)

        # self.fig_density_z_eff = fig_density_z_eff(ax_30, settings_graphics)
        self.fig_phase_difference_z_x1_x2 = fig_phase_difference_z_x1_x2(ax_30, settings_graphics)

        self.fig_control_inputs_of_times = fig_control_inputs_of_times(ax_03, settings_graphics)
        self.fig_global_phase_difference_of_times_analysis = fig_global_phase_difference_of_times_analysis(ax_13, settings_graphics)
        self.fig_number_imbalance_of_times_analysis = fig_number_imbalance_of_times_analysis(ax_23, settings_graphics)
        # =========================================================================================
        
            
        plt.ion()
        
        plt.draw()
        plt.pause(0.001)





    def update_data(self, data):

        V_x = data.V_x

        V_y_x1 = data.V_y_x1

        V_z_x1 = data.V_z_x1

        density_xz = data.density_xz
        density_xy = data.density_xy

        density_x = data.density_x

        density_y_x1 = data.density_y_x1

        density_z_x1 = data.density_z_x1
        density_z_x2 = data.density_z_x2


        phase_difference_xz = data.phase_difference_xz

        phase_difference_z_x1_x2 = data.phase_difference_z_x1_x2


        self.fig_density_xz.update(density_xz)
        self.fig_density_xy.update(density_xy)

        self.fig_density_x.update(density_x, V_x)

        self.fig_density_y.update(density_y_x1, V_y_x1)

        self.fig_density_z.update(density_z_x1, density_z_x2, V_z_x1)

        self.fig_phase_difference_xz.update(phase_difference_xz, density_xz)

        self.fig_phase_difference_z_x1_x2.update(phase_difference_z_x1_x2)



    # def update_data_tof_selected(self, data_tof_selected):
    #
    #     self.fig_tof_xz.update(data_tof_selected.density_tof_xz)
    #
    #     self.fig_tof_xy.update(data_tof_selected.density_tof_xy)
    #
    #     self.fig_tof_profile_x.update(data_tof_selected.profile_tof_x)



    def redraw(self):

        plt.figure(self.fig_name)

        plt.draw()

        self.fig.canvas.start_event_loop(0.001)



    def export(self, filepath):

        plt.figure(self.fig_name)

        plt.draw()

        self.fig.canvas.start_event_loop(0.001)

        plt.savefig(filepath,
                    dpi=None,
                    facecolor='w',
                    edgecolor='w',
                    format='png',
                    transparent=False,
                    bbox_inches=None,
                    pad_inches=0,
                    metadata=None)
