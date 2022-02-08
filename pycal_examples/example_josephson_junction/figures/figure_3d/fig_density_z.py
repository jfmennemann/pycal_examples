from numpy import zeros_like

from numpy import pi

from .. style import colors


class fig_density_z(object):

    def __init__(self, ax, settings):


        self.hbar = settings.hbar

        self.m_atom = settings.m_atom

        self.line_density_z_x1, = ax.plot(settings.z, zeros_like(settings.z), linewidth=1, linestyle='-', color=colors.peter_river, label=r'$|\psi(x_1, 0, z)|^2$')
        self.line_density_z_x2, = ax.plot(settings.z, zeros_like(settings.z), linewidth=1, linestyle='-', color=colors.wet_asphalt, label=r'$|\psi(x_2, 0, z)|^2$')
        
        ax.set_xlim(settings.z_min, settings.z_max)
        
        ax.set_ylim(0, settings.density_max)


        ax.set_xlabel(settings.xlabel_density_z)
        
        ax.set_xticks(settings.z_ticks)
        
        ax.grid(b=True, which='major', color=settings.color_gridlines_major, linestyle='-', linewidth=0.5)
        
        ax.set_ylabel(settings.ylabel_density_x_y_z)

        ax.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0), fancybox=settings.fancybox, framealpha=settings.framealpha, ncol=1)
        
        
        
        ax_V_z_x1 = ax.twinx()
    
        self.line_V_z_x1, = ax_V_z_x1.plot(settings.z, zeros_like(settings.z), linewidth=1, linestyle='-', color=colors.sun_flower, label=r'$V(x_1, 0, z)$')
        
        
        # -----------------------------------------------------------------------------------------
        # V_z_comparison_1 = 0.5 * self.m_atom * settings.vec_omega_z_comparison[0]**2 * (settings.z*1e-6)**2
        # V_z_comparison_2 = 0.5 * self.m_atom * settings.vec_omega_z_comparison[1]**2 * (settings.z*1e-6)**2
        #
        # label_V_z_comparison = '$({0:d}, {1:d})\,$Hz'.format(settings.vec_nu_z_comparison[0], settings.vec_nu_z_comparison[1])
        #
        # ax_V_z.plot(settings.z, V_z_comparison_1 / (self.hbar * 2.0 * np.pi * 1000), linewidth=1, linestyle='--', color=colors.sun_flower, label=label_V_z_comparison)
        # ax_V_z.plot(settings.z, V_z_comparison_2 / (self.hbar * 2.0 * np.pi * 1000), linewidth=1, linestyle='--', color=colors.sun_flower)
        # -----------------------------------------------------------------------------------------

        ax_V_z_x1.set_xlim(settings.z_min, settings.z_max)
        ax_V_z_x1.set_ylim(settings.potential_min, settings.potential_max)
        
        ax_V_z_x1.set_ylabel(settings.ylabel_V_x_y_z)
        
        ax_V_z_x1.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), fancybox=settings.fancybox, framealpha=settings.framealpha, ncol=1)
    
    
    def update(self, density_z_x1, density_z_x2, V_z_x1):

        # print(np.min(density_z))
        # print(np.max(density_z))
        
        scaling_V = (self.hbar * 2 * pi * 1000)
        
        V_z_x1 = V_z_x1 / scaling_V
        
        self.line_density_z_x1.set_ydata(density_z_x1)
        self.line_density_z_x2.set_ydata(density_z_x2)

        self.line_V_z_x1.set_ydata(V_z_x1)
