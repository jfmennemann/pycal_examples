from pycal import Solver3D

import os

import h5py

from time import time

import scipy

from scipy import interpolate

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
simulation_id = 'lattice_3d'

simulation_id = simulation_id.replace(".", "_")

N = 2000

m_Rb_87 = 87 * amu  # kg

m_atom = m_Rb_87

a_s = 5.24e-9  # m

omega_perp = 2 * np.pi * 1e3

x_min = -2e-6
x_max = +2e-6

y_min = -2e-6
y_max = +2e-6

z_min = -20e-6
z_max = +20e-6

Jx = 64
Jy = 64
Jz = 320

t_final = 8e-3
dt = 0.00025e-3

device = 'gpu'
precision = 'double'

seed = 1
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



# =================================================================================================
# init solver
# =================================================================================================

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

# =================================================================================================
# init potential
# =================================================================================================

solver.init_potential(name='lattice',
                      omega_perp=omega_perp,
                      V_z_0_min=0.0,
                      V_z_0_max=15.0 * omega_perp * hbar,
                      k=16)



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
print()



time_1 = time()


# =================================================================================================
# init time evolution
# =================================================================================================

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



# =================================================================================================
# init control inputs
# =================================================================================================

# -------------------------------------------------------------------------------------------------
t0 = 0e-3
t1 = 1e-3
t2 = 2e-3
t3 = 3e-3

u0 = 1.0
u1 = 1.0
u2 = 0.0
u3 = 0.0

vec_t = np.array([t0, t1, t2, t3])
vec_u = np.array([u0, u1, u2, u3])

f = interpolate.PchipInterpolator(vec_t, vec_u)

u1_of_times = f(times)
# -------------------------------------------------------------------------------------------------
# =================================================================================================


# =================================================================================================
# compute ground state solution phi
# =================================================================================================

# -------------------------------------------------------------------------------------------------
u1_0 = u1_of_times[0]

solver.set_V(u1_0)
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
solver.compute_ground_state_solution(N=N, n_iter=5000)

phi = solver.get('phi')

density_0 = np.abs(phi)**2

density_0_max = np.max(density_0)

mu_initial_state = solver.get('mu_phi')

print('mu_phi / h: {0:1.4} kHz'.format(mu_initial_state / (1e3 * (2*pi*hbar))))
# -------------------------------------------------------------------------------------------------


# =================================================================================================
# set wave function psi to ground state solution
# =================================================================================================

solver.init_psi('ground_state_solution')



# =================================================================================================
# init figure
# =================================================================================================

# -------------------------------------------------------------------------------------------------
settings_figure_3d = {
    "density_max": 1.1 * density_0_max,
    "density_z_eff_max": 400.0,
    "V_min": 0.0,
    "V_max": 20.0,
    "sigma_z_min": 0.2,
    "sigma_z_max": 0.6,
    "m_atom": m_Rb_87
}

figure_3d = Figure3d(x, y, z, times, settings_figure_3d)

figure_3d.fig_control_inputs_of_times.update_u(u1_of_times)

figure_3d.fig_control_inputs_of_times.update_t(0.0)
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
data = my_eval(solver)

figure_3d.update_data(data)

figure_3d.redraw()
# -------------------------------------------------------------------------------------------------



# =================================================================================================
# compute time evolution
# =================================================================================================

# -------------------------------------------------------------------------------------------------
# set control input

solver.set_u_of_times(u1_of_times, 1)
# -------------------------------------------------------------------------------------------------

if export_psi_of_times_analysis:

    psi_of_times_analysis = np.zeros((n_times_analysis, Jx, Jy, Jz), dtype=np.complex128)

else:

    psi_of_times_analysis = None


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

    density_z_eff_of_times_analysis[nr_times_analysis, :] = data.density_z_eff

    phase_z_eff_of_times_analysis[nr_times_analysis, :] = data.phase_z_eff
    phase_z_of_times_analysis[nr_times_analysis, :] = data.phase_z


    print('----------------------------------------------------------------------------------------')
    print('t:             {0:1.2f} / {1:1.2f}'.format(t / 1e-3, times[-1] / 1e-3))
    print('n:             {0:4d} / {1:4d}'.format(n, n_times))
    print()
    print('N:             {0:1.4f}'.format(data.N))
    print()
    print('elapsed_time:  {0:1.2f}'.format(elapsed_time))
    print('----------------------------------------------------------------------------------------')
    print()

    # ---------------------------------------------------------------------------------------------
    figure_3d.update_data(data)

    figure_3d.fig_control_inputs_of_times.update_t(t)

    figure_3d.redraw()

    if export_frames:

        filepath = path_frames_figure_3d + 'frame_' + str(nr_frame).zfill(5) + '.png'

        figure_3d.export(filepath)

        nr_frame = nr_frame + 1
    # ---------------------------------------------------------------------------------------------

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
