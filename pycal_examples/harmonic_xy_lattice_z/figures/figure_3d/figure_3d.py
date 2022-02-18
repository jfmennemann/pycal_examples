import matplotlib.pyplot as plt

from PyQt5 import QtWidgets 


import numpy as np


from .fig_density_xz import fig_density_xz
from .fig_density_xy import fig_density_xy

from .fig_density_x import fig_density_x
from .fig_density_y import fig_density_y
from .fig_density_z import fig_density_z

from .fig_real_part_x import fig_real_part_x
from .fig_real_part_y import fig_real_part_y
from .fig_real_part_z import fig_real_part_z

from .fig_density_z_eff import fig_density_z_eff
from .fig_phase_z_eff import fig_phase_z_eff

from .fig_control_inputs_of_times import fig_control_inputs_of_times

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

        m_atom = settings_figure_3d["m_atom"]

        density_min = settings_figure_3d["density_min"]
        density_max = settings_figure_3d["density_max"]

        density_z_eff_min = settings_figure_3d["density_z_eff_min"]
        density_z_eff_max = settings_figure_3d["density_z_eff_max"]

        V_min = settings_figure_3d["V_min"]
        V_max = settings_figure_3d["V_max"]


        Jx = x.size
        Jy = y.size
        Jz = z.size
        
        assert(Jx % 2 == 0)
        assert(Jy % 2 == 0)
        assert(Jz % 2 == 0)
        
        # integer division
        index_center_x = Jx // 2
        index_center_y = Jy // 2
        index_center_z = Jz // 2
        
        
        assert(np.abs(x[index_center_x])<1e-15)
        assert(np.abs(y[index_center_y])<1e-15)
        assert(np.abs(z[index_center_z])<1e-15)
        
        
        
        times = times / 1e-3
        
        x = x / 1e-6
        y = y / 1e-6
        z = z / 1e-6

        x_min = x[0]
        y_min = y[0]
        z_min = z[0]
                
        x_max = -x_min
        y_max = -y_min
        z_max = -z_min
        
        t_min = times[0]
        t_max = times[-1]
        
        

        x_ticks = np.array([-2, -1, 0, 1, 2])
        y_ticks = np.array([-1, 0, 1])


        if np.round(z_max) == 5:
            z_ticks = np.array([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])

        if np.round(z_max) == 10:
            z_ticks = np.array([-10, -5, 0, 5, 10])

        if np.round(z_max) == 20:
            z_ticks = np.array([-20, -10, 0, 10, 20])

        if np.round(z_max) == 40:
            z_ticks = np.array([-40, -30, -20, -10, 0, 10, 20, 30, 40])

        if np.round(z_max) == 50:
            z_ticks = np.array([-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50])

        if np.round(z_max) == 54:
            z_ticks = np.array([-54, 0, 54])

        if np.round(z_max) == 60:
            z_ticks = np.array([-60, -40, -20, 0, 20, 40, 60])

        if np.round(z_max) == 80:
            z_ticks = np.array([-80, -60, -40, -20, 0, 20, 40, 60, 80])

        if np.round(z_max) == 100:
            z_ticks = np.array([-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100])

        if np.round(z_max) == 200:
            z_ticks = np.array([-200, -160, -120, -80, -40, 0, 40, 80, 120, 160, 200])



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

        elif t_max == 200:
            t_ticks_major = np.array([0, 40, 80, 120, 160, 200])

        else:
            t_ticks_major = np.array([0, t_max])



        t_ticks_minor = 0.5 * (t_ticks_major[0:-1] + t_ticks_major[1:])



        settings = type('', (), {})()

        settings.density_min = density_min
        settings.density_max = density_max

        settings.density_z_eff_min = density_z_eff_min
        settings.density_z_eff_max = density_z_eff_max


        settings.hbar = hbar
        settings.m_atom = m_atom





        settings.ylabel_density_x_y_z = r'$\mathrm{density} \;\, \mathrm{in} \;\, \mathrm{m}^{-3}$'
        settings.ylabel_density_x_y_z_effective = r'$\mathrm{density} \;\, \mathrm{in} \;\, \mu \mathrm{m}^{-1}$'




        settings.real_part_min = -np.sqrt(settings.density_max)
        settings.real_part_max = +np.sqrt(settings.density_max)

        settings.V_min = V_min
        settings.V_max = V_max

        settings.phase_difference_of_times_analysis_min = -1.6
        settings.phase_difference_of_times_analysis_max = +1.6

        settings.number_imbalance_of_times_analysis_min = -0.4
        settings.number_imbalance_of_times_analysis_max = +0.4



        settings.times = times
        
        settings.t_min = t_min
        settings.t_max = t_max



        settings.t_ticks_major = t_ticks_major
        settings.t_ticks_minor = t_ticks_minor
        
        
        settings.Jx = Jx
        settings.Jy = Jy
        settings.Jz = Jz
        
        settings.x = x
        settings.y = y
        settings.z = z
        
        settings.x_min = x_min
        settings.y_min = y_min
        settings.z_min = z_min
        
        settings.x_max = x_max
        settings.y_max = y_max
        settings.z_max = z_max
        
        settings.x_ticks = x_ticks
        settings.y_ticks = y_ticks
        settings.z_ticks = z_ticks
        
        
        
        # settings.vec_nu_x_comparison = np.array([1.4e3, 2.1e3])
        # settings.vec_omega_x_comparison = 2 * np.pi * settings.vec_nu_x_comparison
        
        # settings.vec_nu_y_comparison = np.array([1.4e3, 2.1e3])
        # settings.vec_omega_y_comparison = 2 * np.pi * settings.vec_nu_y_comparison
        
        # settings.vec_nu_z_comparison = np.array([7, 11])
        # settings.vec_omega_z_comparison = 2 * np.pi * settings.vec_nu_z_comparison
        
        
        settings.xlabel_phase_difference_of_times_analysis = r'$t             \;\, \mathrm{in} \;\, \mathrm{ms}$'
        settings.ylabel_phase_difference_of_times_analysis = r'$\bar{\varphi} \;\, \mathrm{in} \;\, \mathrm{rad}$'
        
        settings.xlabel_number_imbalance_of_times_analysis = r'$t             \;\, \mathrm{in} \;\, \mathrm{ms}$'
        settings.ylabel_number_imbalance_of_times_analysis = r'$\Delta N$'
        
        
        settings.phase_difference_of_times_analysis_min = -1.6
        settings.phase_difference_of_times_analysis_max = +1.6
        
        settings.number_imbalance_of_times_analysis_min = -0.4
        settings.number_imbalance_of_times_analysis_max = +0.4

        settings.xlabel_profile_tof_x = r'$x \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings.ylabel_profile_tof_x = r'$\mathrm{density} \;\, \mathrm{in} \;\, \mathrm{a.u.}$'
        
        
        
        settings.xlabel_xz = r'$z \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings.ylabel_xz = r'$x \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        
        settings.xlabel_xy = r'$y \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings.ylabel_xy = r'$x \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        
        settings.ylabel_V_x_y_z = r'$V \;\, \mathrm{in} \;\, h \times \mathrm{kHz}$'

        
        settings.xlabel_density_x = r'$x \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings.xlabel_density_y = r'$y \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings.xlabel_density_z = r'$z \;\, \mathrm{in} \;\, \mu \mathrm{m}$'



        
        
        
        
        settings.colormap_density = colors.colormap_1
        
        settings.color_gridlines_major = colors.color_gridlines_major
        settings.color_gridlines_minor = colors.color_gridlines_minor

        settings.framealpha = 1.0
        settings.fancybox = False
        

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
        
        ax_20 = self.fig.add_subplot(self.gridspec[2, 0])
        ax_21 = self.fig.add_subplot(self.gridspec[2, 1])
        ax_22 = self.fig.add_subplot(self.gridspec[2, 2])

        ax_30 = self.fig.add_subplot(self.gridspec[3, 0])
        ax_40 = self.fig.add_subplot(self.gridspec[4, 0])

        self.fig_density_xz = fig_density_xz(ax_00, settings)
        self.fig_density_xy = fig_density_xy(ax_01, settings)

        self.fig_density_z = fig_density_z(ax_10, settings)
        self.fig_density_y = fig_density_y(ax_11, settings)
        self.fig_density_x = fig_density_x(ax_12, settings)

        self.fig_real_part_z = fig_real_part_z(ax_20, settings)
        self.fig_real_part_y = fig_real_part_y(ax_21, settings)
        self.fig_real_part_x = fig_real_part_x(ax_22, settings)

        self.fig_density_z_eff = fig_density_z_eff(ax_30, settings)
        self.fig_phase_z_eff = fig_phase_z_eff(ax_40, settings)

        self.fig_control_inputs_of_times = fig_control_inputs_of_times(ax_03, settings)
        # =========================================================================================
        
            
        plt.ion()
        
        plt.draw()
        plt.pause(0.001)





    def update_data(self, data):

        V_x = data.V_x
        V_y = data.V_y
        V_z = data.V_z

        density_xz = data.density_xz
        density_xy = data.density_xy

        density_x = data.density_x
        density_y = data.density_y
        density_z = data.density_z

        real_part_x = data.real_part_x
        real_part_y = data.real_part_y
        real_part_z = data.real_part_z

        imag_part_x = data.imag_part_x
        imag_part_y = data.imag_part_y
        imag_part_z = data.imag_part_z

        density_z_eff = data.density_z_eff

        phase_z_eff = data.phase_z_eff
        phase_z = data.phase_z

        self.fig_density_xz.update(density_xz)
        self.fig_density_xy.update(density_xy)

        self.fig_density_x.update(density_x, V_x)
        self.fig_density_y.update(density_y, V_y)
        self.fig_density_z.update(density_z, V_z)

        self.fig_real_part_x.update(real_part_x, imag_part_x, V_x)
        self.fig_real_part_y.update(real_part_y, imag_part_y, V_y)
        self.fig_real_part_z.update(real_part_z, imag_part_z, V_z)

        self.fig_density_z_eff.update(density_z_eff, V_z)

        self.fig_phase_z_eff.update(phase_z_eff, phase_z, V_z)

    # def update_data_time_evolution(self,
    #                                data_time_evolution,
    #                                nr_times_analysis):
    #
    #     self.fig_delta_phi_of_times_analysis.update(
    #         data_time_evolution.delta_phi_x1_x2_of_z_psf_circ_mean_00_um_of_times_analysis,
    #         data_time_evolution.delta_phi_x1_x2_of_z_psf_circ_mean_20_um_of_times_analysis,
    #         data_time_evolution.times_analysis,
    #         nr_times_analysis)
    #
    #     self.fig_delta_N_of_times_analysis.update(
    #         data_time_evolution.number_imbalance_mean_of_times_analysis,
    #         data_time_evolution.number_imbalance_std_of_times_analysis,
    #         data_time_evolution.number_imbalance_selected_of_times_analysis,
    #         data_time_evolution.times_analysis,
    #         nr_times_analysis)


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
