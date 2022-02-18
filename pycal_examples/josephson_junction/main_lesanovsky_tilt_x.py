from pycal import Solver3D

import os

import h5py

from time import time

import scipy

import numpy as np

import matplotlib.pyplot as plt

from figures.figure_3d.figure_3d import Figure3d

from eval import my_eval


# -------------------------------------------------------------------------------------------------
pi = scipy.constants.pi

hbar = scipy.constants.hbar

amu = scipy.constants.physical_constants["atomic mass constant"][0]  # atomic mass unit

mu_B = scipy.constants.physical_constants["Bohr magneton"][0]

k_B = scipy.constants.Boltzmann
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
export_hdf5 = False

export_frames = False

export_psi_of_times_analysis = False
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
# close figures from previous simulation

plt.close('all')
# -------------------------------------------------------------------------------------------------


# =================================================================================================
name_potential = 'lesanovsky_tilt_x'

quickstart = False

N = 3500

u1_final = 0.56

gamma_tilt = 4.1 * 1e-26

t_final = 80e-3

T_des = 20e-9

m_Rb_87 = 87 * amu

m_atom = m_Rb_87

a_s = 5.24e-9

g_F = -1/2
m_F = -1
m_F_prime = -1

omega_perp = 2 * np.pi * 3e3
omega_para = 2 * np.pi * 22.5
omega_delta_detuning = -2 * np.pi * 50e3
omega_trap_bottom = 2 * np.pi * 1216e3
omega_rabi_max = 2 * np.pi * 575e3

x_min = -2.8e-6
x_max = +2.8e-6

y_min = -1.2e-6
y_max = +1.2e-6

z_min = -60e-6
z_max = +60e-6

Jx = 2*28
Jy = 2*12
Jz = 4*60

# dt = 0.0025e-3
dt = 0.0020e-3

n_mod_times_analysis = 100

device = 'gpu'
precision = 'double'

seed = 1

visualization = True

settings_figure_3d = {
    'm_atom': m_atom,
    'density_min': -0.2e20,
    'density_max': +2.2e20,
    'V_min': -1.0,
    'V_max': 11.0,
    'abs_z_restr': 30e-6
}
# =================================================================================================

# =================================================================================================
simulation_id = name_potential

simulation_id = simulation_id.replace(".", "_")
# =================================================================================================


# -------------------------------------------------------------------------------------------------
# hdf5

path_f_hdf5 = "./data_hdf5/"

filepath_f_hdf5 = path_f_hdf5 + simulation_id + ".hdf5"
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
# frames

path_frames_figure_3d = "./movies/frames_figure_3d/" + simulation_id + "/"

nr_frame = 0

if export_frames:

    if not os.path.exists(path_frames_figure_3d):
        os.makedirs(path_frames_figure_3d)
# -------------------------------------------------------------------------------------------------


solver = Solver3D(x_min=x_min,
                  x_max=x_max,
                  y_min=y_min,
                  y_max=y_max,
                  z_min=z_min,
                  z_max=z_max,
                  Jx=Jx,
                  Jy=Jy,
                  Jz=Jz,
                  m_atom=m_atom,
                  a_s=a_s,
                  device=device,
                  precision=precision,
                  seed=seed)


solver.init_potential(name=name_potential,
                      g_F=g_F,
                      m_F=m_F,
                      m_F_prime=m_F_prime,
                      omega_perp=omega_perp,
                      omega_para=omega_para,
                      omega_delta_detuning=omega_delta_detuning,
                      omega_trap_bottom=omega_trap_bottom,
                      omega_rabi_max=omega_rabi_max,
                      gamma_tilt=gamma_tilt)


print('solver: seed = {0:d}'.format(solver.get('seed')))
print()
print('solver: pi = {:1.16}'.format(solver.get('pi')))
print('python: pi = {:1.16}'.format(pi))
print()
print('solver: hbar = {:1.16e}'.format(solver.get('hbar')))
print('python: hbar = {:1.16e}'.format(hbar))
print('solver: hbar = {:1.16e} (scaled)'.format(solver.get('hbar', units='solver_units')))
print()
print('solver: mu_B = {:1.16e}'.format(solver.get('mu_B')))
print('python: mu_B = {:1.16e}'.format(mu_B))
print('solver: mu_B = {:1.16e} (scaled)'.format(solver.get('mu_B', units='solver_units')))
print()
print('solver: k_B = {:1.16e}'.format(solver.get('k_B')))
print('python: k_B = {:1.16e}'.format(k_B))
print('solver: k_B = {:1.16e} (scaled)'.format(solver.get('k_B', units='solver_units')))
print()
print('solver: a_s = {:1.16e}'.format(solver.get('a_s')))
print('python: a_s = {:1.16e}'.format(a_s))
print('solver: a_s = {:1.16e} (scaled)'.format(solver.get('a_s', units='solver_units')))
print()
print('m_atom = {:1.16e} (solver: si_units)'.format(solver.get('m_atom')))
print('m_atom = {:1.16e} (python: si_units)'.format(m_atom))
print('m_atom = {:1.16e} (solver: solver_units)'.format(solver.get('m_atom', units='solver_units')))
print()
print('solver: g = {:1.16e}'.format(solver.get('g')))
print('solver: g = {:1.16e} (scaled)'.format(solver.get('g', units='solver_units')))
print()

assert(solver.get('pi') == pi)
assert(solver.get('hbar') == hbar)
assert(solver.get('mu_B') == mu_B)


x = solver.get('x')
y = solver.get('y')
z = solver.get('z')

Jx = solver.get('Jx')
Jy = solver.get('Jy')
Jz = solver.get('Jz')

dx = solver.get('dx')
dy = solver.get('dy')
dz = solver.get('dz')

print('solver: dx = {:1.16e}'.format(solver.get('dx')))
print('solver: dy = {:1.16e}'.format(solver.get('dy')))
print('solver: dz = {:1.16e}'.format(solver.get('dz')))
print()


time_1 = time()

# -------------------------------------------------------------------------------------------------
solver.init_time_evolution(t_final=t_final, dt=dt)

times = solver.get('times')

n_times = times.size
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
times_analysis = times[0::n_mod_times_analysis]

n_times_analysis = times_analysis.size

assert (times_analysis[-1] == t_final)
# -------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
# init control inputs

t_ramp_up = 21.5e-3

t_phase_imprint_part_1 = 1.5e-3
t_phase_imprint_part_2 = 1.5e-3
t_ramp_down = 3.0e-3
t_help = 10.0e-3

t0 = 0.0
t1 = t0 + t_ramp_up
t2 = t1 + t_phase_imprint_part_1
t3 = t2 + t_phase_imprint_part_2
t4 = t3 + t_ramp_down
t5 = t4 + t_help

if quickstart:

    u1_0 = u1_final
    u1_1 = u1_final
    u1_2 = u1_final
    u1_3 = u1_final
    u1_4 = u1_final
    u1_5 = u1_final

else:

    u1_0 = 0.0
    u1_1 = 0.65
    u1_2 = 0.65
    u1_3 = 0.65
    u1_4 = u1_final
    u1_5 = u1_final

vec_t = np.array([t0, t1, t2, t3, t4, t5])

vec_u1 = np.array([u1_0, u1_1, u1_2, u1_3, u1_4, u1_5])
vec_u2 = np.array([0, 0, 1, 0, 0, 0])

u1_of_times = np.interp(times, vec_t, vec_u1)
u2_of_times = np.interp(times, vec_t, vec_u2)
# -------------------------------------------------------------------------------------------------


# =================================================================================================
# compute ground state solution phi
# =================================================================================================

# -------------------------------------------------------------------------------------------------
u1_0 = u1_of_times[0]
u2_0 = u2_of_times[0]

solver.set_V(u1_0, u2_0)
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
solver.compute_ground_state_solution(N=N, n_iter=5000)

phi = solver.get('phi')

density_0 = np.abs(phi)**2

density_0_max = np.max(density_0)

mu_phi = solver.get('mu_phi')

print('mu_phi / h: {0:1.4} kHz'.format(mu_phi / (1e3 * (2*pi*hbar))))
# -------------------------------------------------------------------------------------------------


# =================================================================================================
# set wave function psi to ground state solution
# =================================================================================================

solver.init_psi('ground_state_solution')



# =================================================================================================
# init figure
# =================================================================================================

if visualization:

    # ---------------------------------------------------------------------------------------------
    figure_3d = Figure3d(x, y, z, times, settings_figure_3d)

    figure_3d.fig_control_inputs_of_times.update_u(u1_of_times, u2_of_times)

    figure_3d.fig_control_inputs_of_times.update_t(0.0)
    # ---------------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------
    data = my_eval(solver)

    figure_3d.update_data(data)

    figure_3d.redraw()
    # ---------------------------------------------------------------------------------------------


# =================================================================================================
# thermal state sampling
# =================================================================================================

n_sgpe_max = 10000

n_sgpe_inc = 1000

n_sgpe = 0

while n_sgpe < n_sgpe_max:

    data = my_eval(solver)

    print('----------------------------------------------------------------------------------------')
    print('n_sgpe: {0:4d} / {1:4d}'.format(n_sgpe, n_sgpe_max))
    print()
    print('N:      {0:1.4f}'.format(data.N))
    print('----------------------------------------------------------------------------------------')
    print()

    if visualization:

        # -----------------------------------------------------------------------------------------
        figure_3d.update_data(data)

        figure_3d.redraw()
        # -----------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------
    # apply thermal state sampling process via sgpe for n_sgpe_inc time steps

    solver.iter_sgpe(T_des=T_des, gamma=1e-1, dt=dt, n_inc=n_sgpe_inc)
    # ---------------------------------------------------------------------------------------------

    n_sgpe = n_sgpe + n_sgpe_inc


# =================================================================================================
# compute time evolution
# =================================================================================================

# -------------------------------------------------------------------------------------------------
# set first and second control input

solver.set_u_of_times(u1_of_times, 1)
solver.set_u_of_times(u2_of_times, 2)
# -------------------------------------------------------------------------------------------------

if export_psi_of_times_analysis:

    psi_of_times_analysis = np.zeros((n_times_analysis, Jx, Jy, Jz), dtype=np.complex128)

else:

    psi_of_times_analysis = None



global_phase_difference_of_times_analysis = np.zeros((n_times_analysis,), dtype=np.float64)

number_imbalance_of_times_analysis = np.zeros((n_times_analysis,), dtype=np.float64)



n_inc = n_mod_times_analysis

nr_times_analysis = 0

stop = False

n = 0


while n < n_times-1:

    t = times[n]

    data = my_eval(solver)

    if export_psi_of_times_analysis:

        psi_of_times_analysis[nr_times_analysis, :] = data.psi


    global_phase_difference_of_times_analysis[nr_times_analysis] = data.global_phase_difference

    number_imbalance_of_times_analysis[nr_times_analysis] = data.number_imbalance

    print('----------------------------------------------------------------------------------------')
    print('t: {0:1.2f} / {1:1.2f}'.format(t / 1e-3, times[-1] / 1e-3))
    print('n: {0:4d} / {1:4d}'.format(n, n_times))
    print()
    print('N: {0:1.4f}'.format(data.N))
    print('----------------------------------------------------------------------------------------')
    print()

    if visualization:

        # -----------------------------------------------------------------------------------------
        figure_3d.update_data(data)

        figure_3d.fig_control_inputs_of_times.update_t(t)

        figure_3d.fig_global_phase_difference_of_times_analysis.update(global_phase_difference_of_times_analysis, times_analysis, nr_times_analysis)
        figure_3d.fig_number_imbalance_of_times_analysis.update(number_imbalance_of_times_analysis, times_analysis, nr_times_analysis)

        figure_3d.redraw()

        if export_frames:

            filepath = path_frames_figure_3d + 'frame_' + str(nr_frame).zfill(5) + '.png'

            # figure_3d.export(filepath)

            nr_frame = nr_frame + 1
        # -----------------------------------------------------------------------------------------

    nr_times_analysis = nr_times_analysis + 1


    # ---------------------------------------------------------------------------------------------
    # propagate psi for n_inc time steps

    solver.propagate(n, n_inc)
    # ---------------------------------------------------------------------------------------------

    n = n + n_inc


# -------------------------------------------------------------------------------------------------
# remember: times[n_times - 1] = t_final
# last evaluation of simulation data happens if n == n_times-1

assert(n == n_times-1)

t = times[n]

data = my_eval(solver)

if export_psi_of_times_analysis:

    psi_of_times_analysis[nr_times_analysis, :] = data.psi

print('--------------------------------------------------------------------------------')
print('t: {0:1.2f} / {1:1.2f}'.format(t / 1e-3, times[-1] / 1e-3))
print('n: {0:4d} / {1:4d}'.format(n, n_times))
print()
print('N: {0:1.4f}'.format(data.N))
print('--------------------------------------------------------------------------------')
print()

if visualization:

    # ---------------------------------------------------------------------------------------------
    figure_3d.update_data(data)

    figure_3d.fig_control_inputs_of_times.update_t(t)

    figure_3d.redraw()

    if export_frames:

        filepath = path_frames_figure_3d + 'frame_' + str(nr_frame).zfill(5) + '.png'

        figure_3d.export(filepath)

        nr_frame = nr_frame + 1
    # ---------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------

time_2 = time()

elapsed_time = time_2 - time_1

print("elapsed_time: {0:f}".format(elapsed_time))


if export_hdf5:

    # ---------------------------------------------------------------------------------------------
    # Create file
    f_hdf5 = h5py.File(filepath_f_hdf5, "w")

    # Create file, fail if exists
    # f_hdf5 = h5py.File(filepath_f_hdf5, "x")

    # f_hdf5.create_dataset("trap_id", data=trap_id)
    # ---------------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------
    # f_hdf5.create_dataset("hbar", data=hbar)

    # f_hdf5.create_dataset("N", data=N)
    # f_hdf5.create_dataset("m_atom", data=m_atom)
    # f_hdf5.create_dataset("a_s", data=a_s)

    # f_hdf5.create_dataset("nu_perp", data=nu_perp)
    # f_hdf5.create_dataset("omega_perp", data=omega_perp)

    # f_hdf5.create_dataset("pert_a", data=pert_a)
    # f_hdf5.create_dataset("pert_b", data=pert_b)

    # f_hdf5.create_dataset("x", data=x)
    # f_hdf5.create_dataset("y", data=y)
    # f_hdf5.create_dataset("z", data=z)

    # f_hdf5.create_dataset("z_min", data=z_min)
    # f_hdf5.create_dataset("z_max", data=z_max)

    # f_hdf5.create_dataset("Jz", data=Jz)
    # ---------------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------
    tmp = f_hdf5.create_group("time_evolution")

    if export_psi_of_times_analysis:

        tmp.create_dataset("psi_of_times_analysis", data=psi_of_times_analysis, dtype=np.complex128)

    # tmp.create_dataset("density_z_eff_of_times_analysis", data=density_z_eff_of_times_analysis, dtype=np.float64)

    # tmp.create_dataset("phase_z_eff_of_times_analysis", data=phase_z_eff_of_times_analysis, dtype=np.float64)
    # tmp.create_dataset("phase_z_of_times_analysis", data=phase_z_of_times_analysis, dtype=np.float64)

    tmp.create_dataset("times", data=times)
    tmp.create_dataset("dt", data=dt)
    tmp.create_dataset("n_mod_times_analysis", data=n_mod_times_analysis)
    tmp.create_dataset("times_analysis", data=times_analysis)
    # ---------------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------
    f_hdf5.close()
    # ---------------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------
    f_hdf5 = h5py.File(filepath_f_hdf5, 'r')

    list_all_items = True

    if list_all_items:

        def print_attrs(name, obj):

            print(name)
            for key, val in obj.attrs.items():
                print("    %s: %s" % (key, val))


        f_hdf5.visititems(print_attrs)

        print()
        print()
    # ---------------------------------------------------------------------------------------------

plt.ioff()
plt.show()
