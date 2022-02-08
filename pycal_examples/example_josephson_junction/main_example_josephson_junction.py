from pycal import Solver3d

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


# =================================================================================================
simulation_id = 'lesanovsky_3d'

simulation_id = simulation_id.replace(".", "_")

N = 3500

m_Rb_87 = 87 * amu  # kg

m_atom = m_Rb_87

a_s = 5.24e-9  # m

g_F = -1/2
m_F = -1
m_F_prime = -1

omega_perp = 2 * np.pi * 3e3
omega_para = 2 * np.pi * 22.5
omega_delta_detuning = -2 * np.pi * 50e3
omega_trap_bottom = 2 * np.pi * 1216e3
omega_rabi_max = 2 * np.pi * 575e3
gamma_tilt = 4.1 * 1e-26


x_min = -2.8e-6
x_max = +2.8e-6

y_min = -1.2e-6
y_max = +1.2e-6

z_min = -60e-6
z_max = +60e-6

Jx = 2*28
Jy = 2*12
Jz = 4*60

t_final = 80e-3
dt = 0.002e-3
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

# -------------------------------------------------------------------------------------------------
solver = Solver3d(m_atom=m_atom,
                  a_s=a_s,
                  device='gpu',
                  precision='single',
                  seed=1)

solver.init_grid(x_min=x_min,
                 x_max=x_max,
                 y_min=y_min,
                 y_max=y_max,
                 z_min=z_min,
                 z_max=z_max,
                 Jx=Jx,
                 Jy=Jy,
                 Jz=Jz)

solver.init_V(name='lesanovsky',
              g_F=g_F,
              m_F=m_F,
              m_F_prime=m_F_prime,
              omega_perp=omega_perp,
              omega_para=omega_para,
              omega_delta_detuning=omega_delta_detuning,
              omega_trap_bottom=omega_trap_bottom,
              omega_rabi_max=omega_rabi_max,
              gamma_tilt=gamma_tilt)

solver.init_psi(N)
# -------------------------------------------------------------------------------------------------

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
print('solver: N = {:1.16e}'.format(solver.get('N')))
print('python: N = {:1.16e}'.format(N))
print()
print('solver: a_s = {:1.16e}'.format(solver.get('a_s')))
print('python: a_s = {:1.16e}'.format(a_s))
print('solver: a_s = {:1.16e} (scaled)'.format(solver.get('a_s', units='solver_units')))
print()
print('solver: m_atom = {:1.16e}'.format(solver.get('m_atom')))
print('python: m_atom = {:1.16e}'.format(m_atom))
print('solver: m_atom = {:1.16e} (scaled)'.format(solver.get('m_atom', units='solver_units')))
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




# -------------------------------------------------------------------------------------------------
solver.init_time_evolution(t_final=t_final, dt=dt)

times = solver.get('times')

n_times = times.size
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
n_mod_times_analysis = 100

times_analysis = times[0::n_mod_times_analysis]

n_times_analysis = times_analysis.size

assert (times_analysis[-1] == t_final)
# -------------------------------------------------------------------------------------------------




# -------------------------------------------------------------------------------------------------
# init control inputs

quickstart = False

u1_final = 0.56

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


# -------------------------------------------------------------------------------------------------
plt.close('all')

settings_visualization = {
    "n_mod_times_analysis": n_mod_times_analysis,
    "density_min": 0.0,
    "density_max": 4e20,
    "density_z_eff_min": 0.0,
    "density_z_eff_max": 400.0,
    "V_min": 0.0,
    "V_max": 20.0,
    "sigma_z_min": 0.2,
    "sigma_z_max": 0.6,
    "m_atom": m_Rb_87
}

figure_3d = Figure3d(x, y, z, times, settings_visualization)

figure_3d.fig_control_inputs_of_times.update_u(u1_of_times, u2_of_times)

figure_3d.fig_control_inputs_of_times.update_t(0.0)
# -------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
solver.set_V(0.0, 0.0)
# -------------------------------------------------------------------------------------------------


time_total_start = time()

# =================================================================================================
# compute ground state solution using imaginary time propagation
# =================================================================================================

tau = dt

mu_of_iterations = []
iterations = []

iter_ground_state_inc = 100

iter_ground_state = 0

elapsed_time = 0.0

n_iter_ground_state = 5000

while iter_ground_state < n_iter_ground_state:

    # =============================================================================================
    data = my_eval(solver)

    mu_of_iterations = np.append(mu_of_iterations, data.mu)
    iterations = np.append(iterations, iter_ground_state)

    figure_3d.update_data(data)

    figure_3d.redraw()

    # ---------------------------------------------------------------------------------------------
    fig_temp = plt.figure("figure_ground_state_convergence", figsize=(6, 3), facecolor="white")

    plt.clf()

    ax_1 = fig_temp.add_subplot(111)

    ax_1.plot(iterations[1:], mu_of_iterations[1:] / (2*pi*hbar), linewidth=1.0, linestyle='-', color='k')

    ax_1.set_xlabel('iteration')
    ax_1.set_ylabel('mue / h in Hz')

    plt.tight_layout()

    plt.draw()
    fig_temp.canvas.start_event_loop(0.001)
    # ----------------------------------------------------------------------------------------------

    print('----------------------------------------------------------------------------------------')
    print('ground_state_u_0')
    print()
    print('iter_ground_state: {0:04d}/{1:04d}'.format(iter_ground_state, n_iter_ground_state))
    print()
    print('tau:               {0:1.4e} ms'.format(tau / 1e-3))
    print()
    print('N:                 {0:d}'.format(int(np.round(data.N))))
    print()
    print('mu / h:           {0:1.4} kHz'.format(data.mu / (1e3 * (2*pi*hbar))))
    print()
    print('elapsed_time:      {0:1.4f}'.format(elapsed_time))
    print('----------------------------------------------------------------------------------------')
    print()
    # =============================================================================================

    # =============================================================================================
    time_1 = time()

    solver.iter_gpe_imaginary_time(tau, iter_ground_state_inc)

    time_2 = time()

    elapsed_time = time_2 - time_1

    iter_ground_state = iter_ground_state + iter_ground_state_inc
    # =============================================================================================

    if iter_ground_state == 500 or iter_ground_state == 1000 or iter_ground_state == 1500:

        tau = tau / 2.0

plt.close("figure_ground_state_convergence")





mu_initial_state = solver.get('mu')

print('mu_initial_state / h: {0:1.4} kHz'.format(mu_initial_state / (1e3 * (2*pi*hbar))))




solver.set_u_of_times(u1_of_times, 1)
solver.set_u_of_times(u2_of_times, 2)

# =================================================================================================
# compute time evolution
# =================================================================================================

if export_psi_of_times_analysis:

    psi_of_times_analysis = np.zeros((n_times_analysis, Jx, Jy, Jz), dtype=np.complex128)

else:

    export_psi_of_times_analysis = None



delta_phi_of_times_analysis = np.zeros((n_times_analysis,), dtype=np.float64)

delta_N_of_times_analysis = np.zeros((n_times_analysis,), dtype=np.float64)

density_z_eff_of_times_analysis = np.zeros((n_times_analysis, Jz), dtype=np.float64)

phase_z_eff_of_times_analysis = np.zeros((n_times_analysis, Jz), dtype=np.float64)

phase_z_of_times_analysis = np.zeros((n_times_analysis, Jz), dtype=np.float64)


n_inc = n_mod_times_analysis

nr_times_analysis = 0

stop = False

n = 0

elapsed_time = 0.0


while n < n_times-1:

    t = times[n]

    data = my_eval(solver)

    if export_psi_of_times_analysis:

        psi_of_times_analysis[nr_times_analysis, :] = data.psi


    delta_phi_of_times_analysis[nr_times_analysis] = data.delta_phi

    delta_N_of_times_analysis[nr_times_analysis] = data.delta_N

    density_z_eff_of_times_analysis[nr_times_analysis, :] = data.density_z_eff

    phase_z_eff_of_times_analysis[nr_times_analysis, :] = data.phase_z_eff

    phase_z_of_times_analysis[nr_times_analysis, :] = data.phase_z


    print('----------------------------------------------------------------------------------------')
    print('t:            {0:1.2f} / {1:1.2f}'.format(t / 1e-3, times[-1] / 1e-3))
    print('n:            {0:4d} / {1:4d}'.format(n, n_times))
    print()
    print('N:            {0:1.4f}'.format(data.N))
    print()
    print('elapsed_time: {0:1.2f} s'.format(elapsed_time))
    print('----------------------------------------------------------------------------------------')
    print()

    # ---------------------------------------------------------------------------------------------
    figure_3d.update_data(data)

    figure_3d.fig_control_inputs_of_times.update_t(t)

    figure_3d.fig_delta_phi_of_times_analysis.update(delta_phi_of_times_analysis, times_analysis, nr_times_analysis)
    figure_3d.fig_delta_N_of_times_analysis.update(delta_N_of_times_analysis, times_analysis, nr_times_analysis)

    figure_3d.redraw()

    if export_frames:

        filepath = path_frames_figure_3d + 'frame_' + str(nr_frame).zfill(5) + '.png'

        figure_3d.export(filepath)

        nr_frame = nr_frame + 1
    # ---------------------------------------------------------------------------------------------

    nr_times_analysis = nr_times_analysis + 1



    # ---------------------------------------------------------------------------------------------
    # propagate psi for n_inc time steps

    time_1 = time()

    solver.propagate(n, n_inc)

    time_2 = time()

    elapsed_time = time_2 - time_1
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

density_z_eff_of_times_analysis[nr_times_analysis, :] = data.density_z_eff

phase_z_eff_of_times_analysis[nr_times_analysis, :] = data.phase_z_eff
phase_z_of_times_analysis[nr_times_analysis, :] = data.phase_z

print('--------------------------------------------------------------------------------')
print('t:             {0:1.2f} / {1:1.2f}'.format(t / 1e-3, times[-1] / 1e-3))
print('n:             {0:4d} / {1:4d}'.format(n, n_times))
print()
print('N:             {0:1.4f}'.format(data.N))
print()
print('elapsed_time:  {0:1.2f}'.format(elapsed_time))
print('--------------------------------------------------------------------------------')
print()

# -------------------------------------------------------------------------------------------------
figure_3d.update_data(data)

figure_3d.fig_control_inputs_of_times.update_t(t)

figure_3d.redraw()

if export_frames:

    filepath = path_frames_figure_3d + 'frame_' + str(nr_frame).zfill(5) + '.png'

    figure_3d.export(filepath)

    nr_frame = nr_frame + 1
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


time_total_end = time()

elapsed_time = time_total_end - time_total_start

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

    tmp.create_dataset("density_z_eff_of_times_analysis", data=density_z_eff_of_times_analysis, dtype=np.float64)

    tmp.create_dataset("phase_z_eff_of_times_analysis", data=phase_z_eff_of_times_analysis, dtype=np.float64)
    tmp.create_dataset("phase_z_of_times_analysis", data=phase_z_of_times_analysis, dtype=np.float64)

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

# input("done!")
