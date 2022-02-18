"""
tau = dt

mu_of_iterations = []
iterations = []

# iter_ground_state_inc = n_mod_times_analysis
iter_ground_state_inc = 100

iter_ground_state = 0

elapsed_time = 0.0

n_iter_ground_state = 1000

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

    print('------------------------------------------------')
    print('imaginary time step propagation:')
    print()
    print('iter:         {0:04d}/{1:04d}'.format(iter_ground_state, n_iter_ground_state))
    print()
    print('tau:          {0:1.4e} ms'.format(tau / 1e-3))
    print()
    print('mu/h:         {0:1.4} kHz'.format(data.mu / (1e3 * (2*pi*hbar))))
    print()
    print('N:            {0:d}'.format(int(np.round(data.N))))
    print()
    print('elapsed time: {0:1.4f} s'.format(elapsed_time))
    print('------------------------------------------------')
    print()
    # =============================================================================================

    # =============================================================================================
    time_1 = time()

    solver.imaginary_time_step_iteration(tau, iter_ground_state_inc)

    time_2 = time()

    elapsed_time = time_2 - time_1

    iter_ground_state = iter_ground_state + iter_ground_state_inc
    # =============================================================================================

    if iter_ground_state == 500 or iter_ground_state == 1000 or iter_ground_state == 1500:

        tau = tau / 2.0

plt.close("figure_ground_state_convergence")
"""