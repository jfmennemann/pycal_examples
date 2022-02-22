import matplotlib.pyplot as plt

from PyQt5 import QtWidgets 

import numpy as np

from .. style import colors

from .fig_tof_xz import fig_tof_xz
from .fig_tof_xy import fig_tof_xy

from .fig_tof_profile_x import fig_tof_profile_x


class FigureTof3d(object):

    def __init__(self, x, y, z):
                
        Jx = x.size
        Jy = y.size
        Jz = z.size
        
        assert(Jx % 2 == 0)
        assert(Jy % 2 == 0)
        assert(Jz % 2 == 0)
        
        # integer division
        # index_center_x = Jx // 2
        # index_center_y = Jy // 2
        # index_center_z = Jz // 2
        
        # assert(np.abs(x[index_center_x])<1e-15)
        # assert(np.abs(y[index_center_y])<1e-15)
        # assert(np.abs(z[index_center_z])<1e-15)
        
        
        x = x / 1e-6
        y = y / 1e-6
        z = z / 1e-6
        
        x_min = x[0]
        y_min = y[0]
        z_min = z[0]
        
        x_max = -x_min
        y_max = -y_min
        z_max = -z_min
        
        
        x_ticks = np.array([-2, -1, 0, 1, 2])
        y_ticks = np.array([-1, 0, 1])
        z_ticks = np.array([-60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60])
        
        
        settings = type('', (), {})()
        
        
        settings.x = x
        settings.y = y
        settings.z = z
        
        settings.Jx = Jx
        settings.Jy = Jy
        settings.Jz = Jz
        
        settings.x_min = x_min
        settings.y_min = y_min
        settings.z_min = z_min
                
        settings.x_max = x_max
        settings.y_max = y_max
        settings.z_max = z_max
        
        settings.x_ticks = x_ticks
        settings.y_ticks = y_ticks
        settings.z_ticks = z_ticks
        
        
            
        settings.density_xz = np.zeros((Jx, Jz))
        settings.density_xy = np.zeros((Jx, Jy))
        
        
        
        settings.xlabel_xz            = r'$x                \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings.ylabel_xz            = r'$z                \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        
        settings.xlabel_xy            = r'$x                \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings.ylabel_xy            = r'$y                \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        
        settings.xlabel_profile_tof_x = r'$x \;\, \mathrm{in} \;\, \mu \mathrm{m}$'
        settings.ylabel_profile_tof_x = r'$\mathrm{density} \;\, \mathrm{in} \;\, \mathrm{a.u.}$'
        
        
        
        
        settings.colormap_density = colors.colormap_1
        
        settings.color_gridlines_major = colors.color_gridlines_major
        settings.color_gridlines_minor = colors.color_gridlines_minor
        


        self.fig_name = "figure_tof"
                
        self.fig = plt.figure(self.fig_name, facecolor="white", constrained_layout=False)
        
        self.gridspec = self.fig.add_gridspec(ncols=2, nrows=2, left=0.1, bottom=0.1, right=0.975, top=0.95, wspace=0.5, hspace=0.5, width_ratios=[1, 1], height_ratios=[1, 1])
        
        
        
        
        
        window = self.fig.canvas.window()
        
        window.findChild(QtWidgets.QToolBar).setVisible(False)
        window.statusBar().setVisible(False)
        
        
        
        resolution = '2560x1440'
        # resolution = '1920x1080'
        
        if resolution == '2560x1440':
            
            n_pixels_x = 800
            n_pixels_y = 800
            
            pos_x = 2560 - n_pixels_x
            pos_y = 0
            
            plt.rcParams.update({'font.size': 10})
            
        
        if resolution == '1920x1080':
            
            n_pixels_x = 400
            n_pixels_y = 400
            
            pos_x = 1920 - n_pixels_x
            pos_y = 0
            
            plt.rcParams.update({'font.size': 6})
        
        
        window.setGeometry(pos_x, pos_y, n_pixels_x, n_pixels_y)
        
        
        
        # =========================================================================================
        ax_00 = self.fig.add_subplot(self.gridspec[0, 0])
        ax_10 = self.fig.add_subplot(self.gridspec[1, 0])
        
        ax_01 = self.fig.add_subplot(self.gridspec[0, 1])
        
        
        self.fig_tof_xz = fig_tof_xz(ax_00, settings)
        self.fig_tof_xy = fig_tof_xy(ax_10, settings)
        
        self.fig_tof_profile_x = fig_tof_profile_x(ax_01, settings)
        # =========================================================================================
        
        
        
        plt.ion()
        
        plt.draw()
        plt.pause(0.001)
    
    
    
    def update_data(self, data_tof):

        self.fig_tof_xz.update(data_tof.density_tof_xz, data_tof.x_tof, data_tof.z_tof)
        self.fig_tof_xy.update(data_tof.density_tof_xy, data_tof.x_tof, data_tof.y_tof)

        # self.fig_tof_profile_x.update(data_tof.profile_tof_x, data_tof.x_tof)

        # -----------------------------------------------------------------------------------------
        plt.figure(self.fig_name)
        
        plt.draw()
        
        self.fig.canvas.start_event_loop(0.001)
        # -----------------------------------------------------------------------------------------
