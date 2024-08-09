'''
Author: Edoardo Caldarelli
Affiliation: Institut de Robòtica i Informàtica Industrial, CSIC-UPC
email: ecaldarelli@iri.upc.edu
July 2024
'''
import time

import scipy.linalg
import scipy.signal
import pathlib
import random
import control
import pickle
import scipy.linalg
import scipy.signal
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV

from regressors import *
from dynamical_systems import *

def validate_dyn_sys(regressor: KoopmanRegressor, true_trajectory, test_controls):
    n_samples_traj = true_trajectory.shape[1]
    A = regressor.A
    B = regressor.B
    C = regressor.C
    xcurr_true = true_trajectory[:, 0].reshape([-1, 1])
    phi_xcurr = regressor.lift(xcurr_true)

    xcurr = phi_xcurr
    lifted_traj = xcurr
    simulated_traj = C @ xcurr
    for i in range(0, n_samples_traj - 1):
        xcurr = A @ xcurr + B @ test_controls[:, i].reshape([-1, 1])
        simulated_traj = np.hstack((simulated_traj, C @ xcurr))
        lifted_traj = np.hstack((lifted_traj, xcurr))
    # Compute RMSE
    plt.plot(true_trajectory[0, :])
    plt.plot(simulated_traj[0, :])
    plt.show()
    rmse = np.sqrt(np.sum(np.square(true_trajectory - simulated_traj))) / np.sqrt(np.sum(np.square(simulated_traj))) * 100

    return rmse


def learn_hyperparams(X, Y, kapprox, path_to_data, state_bounds=None):
    if kapprox == 'nystrom':
        cv_regr = KoopmanNystromRegressor(n_inputs)
        # kernels = [LinearKernelWrapper(10 * i) for i in range(-4, 4)]
        kernels = []
        for i in range(-2, 1):
                    kernels.append(KernelWrapper([10**i]))

        clf = GridSearchCV(cv_regr, {'kernel': kernels,
                                    'gamma': np.power(10.0, np.arange(-6, -2, 0.25)),
                                    'm': [500]},
                                    scoring='neg_root_mean_squared_error',
                           verbose=3,
                           n_jobs=-1)
    else:
        cv_regr = KoopmanSplineRegressor(n_inputs, state_bounds)
        clf = GridSearchCV(cv_regr, {'gamma': np.power(10.0, np.arange(-6, -2, 0.25)),
                                                'm': [500]},
                           scoring='neg_root_mean_squared_error',
                           verbose=3,
                           n_jobs=-1)

    clf.fit(X.T, Y.T)
    with open(f"{path_to_data}/cross_validated_kern_params_{kapprox}.npy", 'wb') as f:
        pickle.dump(clf, f)


def lqr_control(num_steps, reference, initial_state, regressor, K):
    A = regressor.A
    B = regressor.B
    C = regressor.C

    phi_new = regressor.lift(initial_state)
    phi_reference = regressor.lift(reference)

    visited_states = np.empty((n_states, 0))
    visited_states = np.hstack((initial_state, visited_states))
    u_s = np.empty((n_inputs, 0))
    x_new = initial_state
    for i in range(0, num_steps):
        if i % 20 == 0:
            print(i)
        u_op = K @ (phi_reference - phi_new)
        u_s = np.hstack((u_s, u_op.reshape(B.shape[1], 1)))
        new_state = C @ phi_new
        visited_states = np.hstack((visited_states, x_new))
        # phi_new = A @ phi_new + B @ u_op
        x_new = dynamical_system.update_SOM(x_new, u_op)
        phi_new = regressor.lift(x_new)
        x_s = visited_states[0, :].T
    return x_s, u_s

def open_loop_control(dynamical_system, initial_state, controls):
    Jcurr = initial_state ** 2
    state = initial_state
    states = state.reshape([-1, 1])
    for i in range(0, controls.shape[1]):
        state = dynamical_system.update_SOM(state, controls[:, i])
        Jcurr = Jcurr + state ** 2 + controls[:, i] ** 2
        states = np.hstack((states, state))
    return states, Jcurr


def generate_dataset(dynamical_system: DynamicalSystem, n_trajs: int, n_samples_traj: int):
    X = np.zeros((dynamical_system.n_states + 1, n_trajs * n_samples_traj))
    Y = np.zeros((dynamical_system.n_states, n_trajs * n_samples_traj))
    indx = 0
    for i in range(0, n_trajs):
        x_curr = np.random.uniform(dynamical_system.state_lb, dynamical_system.state_ub)
        curr_traj = np.array([x_curr]).reshape([-1, 1])
        for j in range(0, n_samples_traj):
            u_curr = np.random.uniform(dynamical_system.input_lb, dynamical_system.input_ub).reshape([dynamical_system.n_inputs, 1])
            curr_augm_state = np.vstack((x_curr, u_curr))
            X[:, indx] = np.squeeze(curr_augm_state)
            x_curr = dynamical_system.update_SOM(x_curr, u_curr).reshape([-1, 1])
            curr_traj = np.hstack((curr_traj, x_curr))
            Y[:, indx] = np.squeeze(x_curr)
            indx = indx + 1
    return X, Y


def simulate_true_system(dynamical_system: DynamicalSystem, T):
    initial_state = np.random.uniform(dynamical_system.state_lb, dynamical_system.state_ub)
    times = np.linspace(0, T, int(1 / dynamical_system.Ts))
    u_s = 2 * times
    state = np.array(initial_state).reshape([-1, 1])
    visited_states = state.reshape([-1, 1])
    for u in u_s:
        state = dynamical_system.update_SOM(state, u)
        visited_states = np.hstack((visited_states, state.reshape([-1, 1])))
    # plt.plot(visited_states[0, :])
    # plt.show()
    return visited_states, u_s.reshape([dynamical_system.n_inputs, -1])


if __name__ == '__main__':
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rcParams.update({'font.size': 22})
    system = "hjb"
    gen_dataset = False
    cross_validate = False
    validate_sys_id = False  # Set to false to test the LQR, true to test the open loop forecasts
    test_exact_regr = True
    dynamical_system = None
    if system == 'hjb':
        model_params = {'Ts': 0.01,
                        'name': system,
                        'n_states': 1,
                        'n_inputs': 1,
                        'state_lb': -1.0,
                        'state_ub': 1.0,
                        'input_lb': [-1],
                        'input_ub': [1]}
        dynamical_system = HJB(**model_params)
    else:
        print("Invalid system ID!")
        exit(1)

    n_trajs = 20
    simulation_horizon = 2  # s
    n_samples_traj = int(simulation_horizon // dynamical_system.Ts)

    np.random.seed(0)  # Fix seed for dataset creation
    random.seed(0)
    path_to_data = pathlib.Path(f"./{system}")
    path_to_data.mkdir(exist_ok=True)

    if gen_dataset:
        X, Y = generate_dataset(dynamical_system, n_trajs, n_samples_traj)
        with open(f"{path_to_data}/dataset.npy", 'wb') as f:
            pickle.dump((X, Y), f)
    with open(f"{path_to_data}/dataset.npy", 'rb') as f:
        X, Y = pickle.load(f)

    ms = np.around(np.logspace(1, 2.3, num=20)).astype(int)  # np.arange(10, 200, 100)
    # ms = ms[:-3]
    ms = np.array([100])
    # ms = np.logspace(1, 2.6, num=5, dtype=int)

    # ms = np.array([20])
    labels = ['nystrom']
    n_inputs = dynamical_system.n_inputs
    n_states = dynamical_system.n_states
    test_trajectories = []
    test_controls = []
    n_seeds = 200
    for seed in range(0, n_seeds):
        np.random.seed(seed)  # Fix seed
        random.seed(seed)
        test_trajectory, test_control = simulate_true_system(dynamical_system, simulation_horizon)

        test_trajectories.append(test_trajectory)
        test_controls.append(test_control)

    plt.show()

    for k, kapprox in enumerate(labels):
        if cross_validate:
            n_val_trajs = 20
            n_samples_val_trajs = int(simulation_horizon // dynamical_system.Ts)
            Xval, Yval = generate_dataset(dynamical_system, n_val_trajs, n_samples_val_trajs)
            if kapprox == 'nystrom':
                learn_hyperparams(Xval, Yval, kapprox, path_to_data)
            else:
                learn_hyperparams(Xval, Yval, kapprox, path_to_data, state_bounds_params)


    if validate_sys_id:
        for k, kapprox in enumerate(labels):
            all_rmse_across_seeds = np.empty((0, ms.shape[0]))
            with open(f"{path_to_data}/cross_validated_kern_params_nystrom.npy", 'rb') as f:
                clf = pickle.load(f)
                kernel_params = clf.best_params_
            print(kernel_params)
            for seed in range(0, n_seeds):
                np.random.seed(seed)  # Fix seed
                random.seed(seed)

                all_rmses = []  # Collect RMSEs across all testing trajectories
                all_ls = [1, 1]
                rmses = []

                for m_indx, m in enumerate(ms):
                    print('m ', m, " seed ", seed)

                    # print("Curr number of features: ", m)
                    regressor = None
                    if kapprox == 'nystrom':
                        # kernel_params['kernel'] = KernelWrapper([1] * n_states)
                        kernel_params['m'] = m
                        # print(kernel_params)
                        regressor = KoopmanNystromRegressor(n_inputs, **kernel_params)
                    elif kapprox == 'splines':
                        kernel_params['m'] = m
                        regressor = KoopmanSplineRegressor(n_inputs, state_bounds_params=state_bounds_params,
                                                           **kernel_params)
                    else:
                        del kernel_params['m']
                        regressor = KoopmanKernelRegressor(n_inputs,
                                                           **kernel_params)
                    regressor.fit(X.T,
                                  Y.T)  # Careful with transposition (shape required by sklearn estimator API, used for CV)
                    # Evaluate the prediction accuracy in open loop
                    test_trajectory = test_trajectories[seed]
                    test_control = test_controls[seed]
                    curr_rmse = validate_dyn_sys(regressor, test_trajectory, test_control)
                    # plt.scatter(test_trajectory[0, 0], test_trajectory[1, 0], color='C3', s=curr_rmse)
                    rmses.append(curr_rmse)
                print(rmses)
                all_rmses.append(rmses)
                all_rmses_array = np.array(all_rmses)
                all_rmse_across_seeds = np.vstack((all_rmse_across_seeds, all_rmses_array))
            print(np.median(all_rmse_across_seeds, axis=0))
            print(np.percentile(all_rmse_across_seeds, axis=0, q=15))
            print(np.percentile(all_rmse_across_seeds, axis=0, q=85))
            # np.savetxt(f"{path_to_data}/all_rmses_{kapprox}_double_dataset.csv", all_rmse_across_seeds)
            plt.show()
    else:
        first_states_nys = []
        for k, kapprox in enumerate(labels):
            with open(f"{path_to_data}/cross_validated_kern_params_nystrom.npy", 'rb') as f:
                clf = pickle.load(f)
                kernel_params = clf.best_params_
            all_rmses_control_across_ms = []
            if kapprox == 'nystrom':
                for m in ms:
                    print("m:", m)
                    rmses_control = []
                    for seed in range(0, n_seeds):
                        np.random.seed(seed)  # Fix seed
                        random.seed(seed)
                        print("Seed ", seed)
                        kernel_params['m'] = m
                        regressor = KoopmanNystromRegressor(n_inputs, **kernel_params)
                        regressor.fit(X.T,
                                      Y.T)
                        A = regressor.A
                        B = regressor.B
                        C = regressor.C
                        initial_state = np.array([0.9]).reshape([-1, 1])

                        R = np.eye(n_inputs)
                        Q = C.T @ C
                        reference = np.array([0.0]).reshape([-1, 1])

                        K, _, E = control.dlqr(A, B, Q, R)
                        end = int(10 * simulation_horizon / dynamical_system.Ts)
                        xs, us = lqr_control(end, reference, initial_state, regressor, K)

                        states, Jcurr = open_loop_control(dynamical_system, initial_state, us)
                        np.savetxt(f"{path_to_data}/Jcurr_m_{m}_seed_{seed}.csv", Jcurr)

                        if kapprox == 'nystrom':
                            first_states_nys.append(states[0, :])
                        # plt.plot(states[0, :], states[1, :])
                        # plt.show()
                        Jtrue = initial_state ** 2
                        xnew_true = initial_state
                        u_opt_true = []
                        for i in range(0, int(end)):
                            u_opt_curr = xnew_true ** 3 - xnew_true * np.sqrt(1 + xnew_true ** 4)
                            u_opt_true.append(u_opt_curr.squeeze())
                            xnew_true = dynamical_system.update_SOM(xnew_true, u_opt_curr)
                            Jtrue = Jtrue + xnew_true ** 2 + u_opt_curr ** 2
                        np.savetxt(f"{path_to_data}/Jtrue.csv", Jtrue)
                        rmse_control = np.sqrt(np.sum(np.square(us.squeeze() - u_opt_true))) / np.sqrt(np.sum(np.square(u_opt_true))) * 100
                        rmses_control.append(rmse_control)
                        if m == 100 and seed == 0:
                            np.savetxt("hjb/control_nystrom_kernel.csv", us)
                            np.savetxt("hjb/control_optimal.csv", u_opt_true)

                        # plt.figure(figsize=[10, 4])
                        # plt.plot(np.arange(0, 10 * simulation_horizon, dynamical_system.Ts), u_opt_true, linewidth=2)
                        # plt.plot(np.arange(0, 10 * simulation_horizon, dynamical_system.Ts), us.squeeze(), linewidth=3, linestyle='-.')
                        # plt.legend(["True opt.\ control", "Nys.\ Matérn-5/2"], bbox_to_anchor=(0.0, 1.02, 1.0, 0.2),
                        #            loc='lower left',
                        #            mode='expand',
                        #            borderaxespad=0, ncol=2, handlelength=1.0)
                        # plt.xlabel("t [s]")
                        # plt.grid(visible=True, which='both')
                        # plt.tight_layout()
                        # plt.savefig(f"{path_to_data}/control_comparison_hjb.png", dpi=300, bbox_inches='tight',
                        #             pad_inches=0)
                        # plt.show()
                    print(rmses_control)
                    all_rmses_control_across_ms.append(rmses_control)
            else:
                np.random.seed(0)  # Fix seed
                random.seed(0)
                print("Seed ", 0)
                del kernel_params['m']
                print(kapprox)
                print('starting fit ...')
                regressor = KoopmanKernelRegressor(n_inputs, **kernel_params)
                regressor.fit(X.T,
                              Y.T)
                print('... fit completed')

                A = regressor.A
                B = regressor.B
                C = regressor.C
                initial_state = np.array([0.9]).reshape([-1, 1])

                R = np.eye(n_inputs)
                Q = C.T @ C
                reference = np.array([0.0]).reshape([-1, 1])
                print('computing lqr gain...')
                now = time.time()
                K, _, E = control.dlqr(A, B, Q, R)
                print('... gain computed', time.time() - now)
                np.savetxt('Kexact.csv', K)
                end = int(10 * simulation_horizon / dynamical_system.Ts)
                xs, us = lqr_control(end, reference, initial_state, regressor, K)

                states, Jcurr = open_loop_control(dynamical_system, initial_state, us)
                np.savetxt(f"{path_to_data}/Jcurr_EXACT_seed_{seed}.csv", Jcurr)

                if kapprox == 'nystrom':
                    first_states_nys.append(states[0, :])
                # plt.plot(states[0, :], states[1, :])
                # plt.show()
                Jtrue = initial_state ** 2
                xnew_true = initial_state
                u_opt_true = []
                for i in range(0, int(end)):
                    u_opt_curr = xnew_true ** 3 - xnew_true * np.sqrt(1 + xnew_true ** 4)
                    u_opt_true.append(u_opt_curr.squeeze())
                    xnew_true = dynamical_system.update_SOM(xnew_true, u_opt_curr)
                    Jtrue = Jtrue + xnew_true ** 2 + u_opt_curr ** 2
                np.savetxt(f"{path_to_data}/Jtrue.csv", Jtrue)
                rmse_control = np.sqrt(np.sum(np.square(us.squeeze() - u_opt_true))) / np.sqrt(
                    np.sum(np.square(u_opt_true))) * 100
                print(rmse_control)
                np.savetxt("hjb/control_exact_kernel.csv", us)
                # plt.figure(figsize=[10, 4])
                # plt.plot(np.arange(0, 10 * simulation_horizon, dynamical_system.Ts), u_opt_true, linewidth=2)
                # plt.plot(np.arange(0, 10 * simulation_horizon, dynamical_system.Ts), us.squeeze(), linewidth=3, linestyle='-.')
                # plt.legend(["True opt.\ control", "Nys.\ Matérn-5/2"], bbox_to_anchor=(0.0, 1.02, 1.0, 0.2),
                #            loc='lower left',
                #            mode='expand',
                #            borderaxespad=0, ncol=2, handlelength=1.0)
                # plt.xlabel("t [s]")
                # plt.grid(visible=True, which='both')
                # plt.tight_layout()
                # plt.savefig(f"{path_to_data}/control_comparison_hjb.png", dpi=300, bbox_inches='tight',
                #             pad_inches=0)
                # plt.show()

                # np.savetxt(f"{path_to_data}/control_rmse.csv", all_rmses_control_across_ms)
        # first_states_nys = np.array(first_states_nys)
        # second_states_nys = np.array(second_states_nys)
        #
        # # np.savetxt(f"{path_to_experiment}/first_state_control_splines.csv", first_states_splines)
        # # np.savetxt(f"{path_to_experiment}/first_state_control_nystrom.csv", first_states_nys)
        # #
        # # np.savetxt(f"{path_to_experiment}/second_state_control_splines.csv", second_states_splines)
        # # np.savetxt(f"{path_to_experiment}/second_state_control_nystrom.csv", second_states_nys)
        # #
        # # first_states_splines = np.loadtxt(f"{path_to_experiment}/first_state_control_splines.csv", )
        # # first_states_nys = np.loadtxt(f"{path_to_experiment}/first_state_control_nystrom.csv", )
        # #
        # # second_states_splines = np.loadtxt(f"{path_to_experiment}/second_state_control_splines.csv", )
        # # second_states_nys = np.loadtxt(f"{path_to_experiment}/second_state_control_nystrom.csv", )
        # fig = plt.figure(figsize=(8, 4))
        #
        #
        # plt.plot(np.arange(0, first_states_nys.shape[1]) * dynamical_system.Ts, np.median(first_states_nys, axis=0), color=f'C0', linewidth=2)
        # plt.fill_between(np.arange(0, first_states_nys.shape[1]) * dynamical_system.Ts, np.percentile(first_states_nys, axis=0, q=15), np.percentile(first_states_nys, axis=0, q=85), alpha=0.3, color='C0')
        # plt.xlabel('t [s]')
        # plt.ylabel('$x_1$')
        # plt.grid(visible=True, which='both')
        #
        # plt.legend(["Nyström Matérn-5/2"], bbox_to_anchor=(0.0, 1.02, 1.0, 0.2), loc='lower left',
        #            mode='expand',
        #            borderaxespad=0, ncol=3, handlelength=1.0)
        # plt.tight_layout()
        # plt.savefig(f"{path_to_data}/control_x1.png", dpi=300, bbox_inches='tight',
        #             pad_inches=0)
        # plt.show()
        # fig = plt.figure(figsize=(8, 4))
        #
        # plt.plot(np.arange(0, second_states_nys.shape[1]) * dynamical_system.Ts, np.median(second_states_nys, axis=0), color=f'C0', linewidth=2)
        # plt.fill_between(np.arange(0, second_states_nys.shape[1]) * dynamical_system.Ts, np.percentile(second_states_nys, axis=0, q=15), np.percentile(second_states_nys, axis=0, q=85), alpha=0.3, color='C0')
        # plt.xlabel('t [s]')
        # plt.ylabel('$x_2$')
        # plt.grid(visible=True, which='both')
        #
        # plt.legend(["Nyström Matérn-5/2"], bbox_to_anchor=(0.0, 1.02, 1.0, 0.2), loc='lower left',
        #            mode='expand',
        #            borderaxespad=0, ncol=3, handlelength=1.0)
        # plt.tight_layout()
        # plt.savefig(f"{path_to_data}/control_x2.png", dpi=300, bbox_inches='tight',
        #             pad_inches=0)
        # plt.show()



