import numpy as np


# def get_indices_x1_x2(V_x, Jx, index_center_x):
#
#     index_x1 = -1
#
#     for jx in range(1, index_center_x + 1):
#
#         if V_x[jx - 1] >= V_x[jx] and V_x[jx] <= V_x[jx + 1]:
#             index_x1 = jx
#             break
#
#     assert (index_x1 > 0)
#
#     index_x2 = -1
#
#     for jx in range(index_center_x, Jx):
#
#         if V_x[jx - 1] >= V_x[jx] and V_x[jx] <= V_x[jx + 1]:
#             index_x2 = jx
#             break
#
#     assert (index_x2 > 0)
#
#     return index_x1, index_x2


# def compute_psi_complete(psi, fill_boundaries=False):
#
#     Jx = psi.shape[0]
#     Jy = psi.shape[1]
#     Jz = psi.shape[2]
#
#     psi_complete = np.zeros((Jx+1, Jy+1, Jz+1), dtype=np.complex128)
#
#     psi_complete[:Jx, :Jy, :Jz] = psi
#
#     if fill_boundaries:
#
#         psi_complete[-1, :, :] = psi_complete[0, :, :]
#         psi_complete[:, -1, :] = psi_complete[:, 0, :]
#         psi_complete[:, :, -1] = psi_complete[:, :, 0]
#
#     return psi_complete


# def compute_phase_difference(psi):
#
#     psi_complete = compute_psi_complete(psi)
#
#     psi_complete_flip_x = np.flip(psi_complete, 0)
#
#     tmp = psi_complete * np.conj(psi_complete_flip_x)
#
#     phase_difference_complete = np.angle(tmp)
#
#     phase_difference = phase_difference_complete[:-1, :-1, :-1]
#
#     return phase_difference



# def compute_phase_difference_z_x1_x2(psi_z_x1, psi_z_x2):
#
#     phase_difference_z_x1_x2 = np.angle(psi_z_x1 * np.conj(psi_z_x2))
#
#     return phase_difference_z_x1_x2


# def compute_global_phase_difference(psi, index_center_x):
#
#     psi_complete = compute_psi_complete(psi)
#
#     psi_complete_flip_x = np.flip(psi_complete, 0)
#
#     tmp = psi_complete * np.conj(psi_complete_flip_x)
#
#     tmp = np.sum(tmp[:index_center_x, :, :])
#
#     delta_phi = np.angle(tmp)
#
#     return delta_phi


# def compute_number_imbalance(psi, dx, dy, dz, index_center_x):
#
#     density = np.real(psi * np.conj(psi))
#
#     N = (dx * dy * dz) * np.sum(density)
#
#     density_1 = density[:index_center_x+1, :, :]
#     density_2 = density[index_center_x:, :, :]
#
#     N_1 = (dx * dy * dz) * np.sum(density_1)
#     N_2 = (dx * dy * dz) * np.sum(density_2)
#
#     number_imbalance = (N_2 - N_1) / N
#
#     return number_imbalance



def my_eval_tof(solver):

    data_tof = type('', (), {})()

    # -----------------------------------------------------------------------------------------------------------------
    x_tof_stage_1 = solver.get('x_tof_stage_1')
    y_tof_stage_1 = solver.get('y_tof_stage_1')
    z_tof_stage_1 = solver.get('z_tof_stage_1')

    Jx_tof_stage_1 = x_tof_stage_1.size
    Jy_tof_stage_1 = y_tof_stage_1.size
    Jz_tof_stage_1 = z_tof_stage_1.size

    psi_tof_stage_1 = solver.get('psi_tof_stage_1')

    psi_xz_tof_stage_1 = np.squeeze(psi_tof_stage_1[:, Jy_tof_stage_1//2, :])
    psi_xy_tof_stage_1 = np.squeeze(psi_tof_stage_1[:, :, Jz_tof_stage_1//2])

    data_tof.density_xz_tof_stage_1 = np.abs(psi_xz_tof_stage_1)**2
    data_tof.density_xy_tof_stage_1 = np.abs(psi_xy_tof_stage_1)**2

    data_tof.x_tof_stage_1 = x_tof_stage_1
    data_tof.y_tof_stage_1 = y_tof_stage_1
    data_tof.z_tof_stage_1 = z_tof_stage_1

    data_tof.Jx_tof_stage_1 = Jx_tof_stage_1
    data_tof.Jy_tof_stage_1 = Jy_tof_stage_1
    data_tof.Jz_tof_stage_1 = Jz_tof_stage_1
    # -----------------------------------------------------------------------------------------------------------------

    return data_tof
