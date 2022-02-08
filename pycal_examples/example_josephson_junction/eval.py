import numpy as np


def get_indices_x1_x2(V_x, Jx, index_center_x):

    index_x1 = -1

    for jx in range(1, index_center_x + 1):

        if V_x[jx - 1] >= V_x[jx] and V_x[jx] <= V_x[jx + 1]:
            index_x1 = jx
            break

    assert (index_x1 > 0)

    index_x2 = -1

    for jx in range(index_center_x, Jx):

        if V_x[jx - 1] >= V_x[jx] and V_x[jx] <= V_x[jx + 1]:
            index_x2 = jx
            break

    assert (index_x2 > 0)

    return index_x1, index_x2


def compute_psi_complete(psi, fill_boundaries=False):

    Jx = psi.shape[0]
    Jy = psi.shape[1]
    Jz = psi.shape[2]

    psi_complete = np.zeros((Jx+1, Jy+1, Jz+1), dtype=np.complex128)

    psi_complete[:Jx, :Jy, :Jz] = psi

    if fill_boundaries:

        psi_complete[-1, :, :] = psi_complete[0, :, :]
        psi_complete[:, -1, :] = psi_complete[:, 0, :]
        psi_complete[:, :, -1] = psi_complete[:, :, 0]

    return psi_complete


def compute_phase_difference(psi):

    psi_complete = compute_psi_complete(psi)

    psi_complete_flip_x = np.flip(psi_complete, 0)

    tmp = psi_complete * np.conj(psi_complete_flip_x)

    phase_difference_complete = np.angle(tmp)

    phase_difference = phase_difference_complete[:-1, :-1, :-1]

    return phase_difference



def compute_phase_difference_z_x1_x2(psi_z_x1, psi_z_x2):

    phase_difference_z_x1_x2 = np.angle(psi_z_x1 * np.conj(psi_z_x2))

    return phase_difference_z_x1_x2


def compute_global_phase_difference(psi, index_center_x):

    psi_complete = compute_psi_complete(psi)

    psi_complete_flip_x = np.flip(psi_complete, 0)

    tmp = psi_complete * np.conj(psi_complete_flip_x)

    tmp = np.sum(tmp[:index_center_x, :, :])

    delta_phi = np.angle(tmp)

    return delta_phi


def compute_number_imbalance(psi, dx, dy, dz, index_center_x):

    density = np.real(psi * np.conj(psi))

    N = (dx * dy * dz) * np.sum(density)

    density_1 = density[:index_center_x+1, :, :]
    density_2 = density[index_center_x:, :, :]

    N_1 = (dx * dy * dz) * np.sum(density_1)
    N_2 = (dx * dy * dz) * np.sum(density_2)

    number_imbalance = (N_2 - N_1) / N

    return number_imbalance



def my_eval(solver):

    dx = solver.get('dx')
    dy = solver.get('dy')
    dz = solver.get('dz')

    Jx = solver.get('Jx')

    index_center_x = solver.get('index_center_x')
    index_center_y = solver.get('index_center_y')
    index_center_z = solver.get('index_center_z')

    data = type('', (), {})()

    # -----------------------------------------------------------------------------------------------------------------
    V = solver.get('V')

    V_x = np.squeeze(V[:, index_center_y, index_center_z])

    index_x1, index_x2 = get_indices_x1_x2(V_x, Jx, index_center_x)

    data.V_x = V_x

    data.V_y_x1 = np.squeeze(V[index_x1, :, index_center_z])
    data.V_y_x2 = np.squeeze(V[index_x2, :, index_center_z])

    data.V_z_x1 = np.squeeze(V[index_x1, index_center_y, :])
    data.V_z_x2 = np.squeeze(V[index_x2, index_center_y, :])
    # -----------------------------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------------------------
    psi = solver.get('psi')

    psi_x = np.squeeze(psi[:, index_center_y, index_center_z])

    psi_y_x1 = np.squeeze(psi[index_x1, :, index_center_z])
    psi_y_x2 = np.squeeze(psi[index_x2, :, index_center_z])

    psi_z_x1 = np.squeeze(psi[index_x1, index_center_y, :])
    psi_z_x2 = np.squeeze(psi[index_x2, index_center_y, :])

    psi_xz = np.squeeze(psi[:, index_center_y, :])
    psi_xy = np.squeeze(psi[:, :, index_center_z])

    data.psi = psi

    data.psi_x = psi_x

    data.psi_y_x1 = psi_y_x1
    data.psi_y_x2 = psi_y_x2

    data.psi_z_x1 = psi_z_x1
    data.psi_z_x2 = psi_z_x2

    data.psi_xz = psi_xz
    data.psi_xy = psi_xy
    # -----------------------------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------------------------
    data.density_x = np.abs(psi_x) ** 2

    data.density_y_x1 = np.abs(psi_y_x1) ** 2
    data.density_y_x2 = np.abs(psi_y_x2) ** 2

    data.density_z_x1 = np.abs(psi_z_x1) ** 2
    data.density_z_x2 = np.abs(psi_z_x2) ** 2

    data.density_xz = np.abs(psi_xz) ** 2
    data.density_xy = np.abs(psi_xy) ** 2
    # -----------------------------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------------------------
    data.N = solver.get('N')
    # -----------------------------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------------------------
    data.mu = solver.get('mu')
    # -----------------------------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------------------------
    data.global_phase_difference = compute_global_phase_difference(psi, index_center_x)
    # -----------------------------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------------------------
    data.number_imbalance = compute_number_imbalance(psi, dx, dy, dz, index_center_x)
    # -----------------------------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------------------------
    phase_difference = compute_phase_difference(psi)

    phase_difference_xz = np.squeeze(phase_difference[:, index_center_y, :])

    data.phase_difference_xz = phase_difference_xz
    # -----------------------------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------------------------
    data.phase_difference_z_x1_x2 = compute_phase_difference_z_x1_x2(psi_z_x1, psi_z_x2)
    # -----------------------------------------------------------------------------------------------------------------

    return data
