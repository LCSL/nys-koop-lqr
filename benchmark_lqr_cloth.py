'''
Author: Edoardo Caldarelli
Affiliation: Institut de Robòtica i Informàtica Industrial, CSIC-UPC
email: ecaldarelli@iri.upc.edu
January 2024
'''


import pathlib
import random
import matplotlib; matplotlib.use("TkAgg")
import control
import pickle
from sklearn.model_selection import GridSearchCV
from regressors import ThreeDimensionalKernel, KoopmanNystromRegressor, KoopmanRegressor, KoopmanSplineRegressor
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
    rmse = np.sqrt(np.mean(np.square(true_trajectory - simulated_traj)))

    return rmse


def learn_hyperparams(X, Y, kapprox):
    ls = np.power(10.0, np.arange(0, 3.5))
    if kapprox == 'nystrom':
        cv_regr = KoopmanNystromRegressor(n_inputs)
    elif kapprox == 'splines':
        cv_regr = KoopmanSplineRegressor(n_inputs)
    if kapprox == 'nystrom':
        kernel_candidates = []
        max_exp_lth = 3
        for i in range(0, max_exp_lth):
            for j in range(0, max_exp_lth):
                for k in range(0, max_exp_lth):
                    kernel_candidates.append(ThreeDimensionalKernel(10**i, 10**j, 10**k, Y.shape[0]))
        clf = GridSearchCV(cv_regr, {'kernel': kernel_candidates,
                                    'gamma': np.power(10.0, np.arange(-7, -4)),
                                    'm': [500]},
                                    scoring='neg_root_mean_squared_error',
                           verbose=3,
                           n_jobs=-1)
    else:
        clf = GridSearchCV(cv_regr, {
            'gamma': np.power(10.0, np.arange(-7, -4)),
            'm': [500]},
                           scoring='neg_root_mean_squared_error',
                           verbose=3,
                           n_jobs=-1)
    cv_result = clf.fit(X.T, Y.T)
    return clf


def lqr_control(num_steps, reference, initial_state, regressor, K):
    A = regressor.A
    B = regressor.B
    C = regressor.C
    phi_new = regressor.lift(initial_state)
    phi_reference = regressor.lift(reference)
    visited_states = np.empty((n_states, 0))
    visited_states = np.hstack((initial_state, visited_states))
    u_s = np.empty((n_inputs, 0))
    u_s = np.hstack((u_s, initial_state[[168, 169, 170, 189, 190, 191], :]))
    for i in range(0, num_steps):
        u_op = K @ (phi_reference - phi_new)
        u_s = np.hstack((u_s, u_s[:, -1].reshape([-1, 1]) + u_op))
        new_state = C @ phi_new
        visited_states = np.hstack((visited_states, new_state))
        phi_new = A @ phi_new + B @ u_op
    n_nodes = int(n_states / 3)
    x_s = np.zeros((n_nodes, visited_states.shape[1]))
    y_s = np.zeros((n_nodes, visited_states.shape[1]))
    z_s = np.zeros((n_nodes, visited_states.shape[1]))
    for j in range(visited_states.shape[1]):
        kx = 0
        ky = 0
        kz = 0
        for i in range(0, n_states):
            if i % 3 == 0:
                x_s[kx, j] = visited_states[i, j]
                kx += 1
            if i % 3 == 1:
                y_s[ky, j] = visited_states[i, j]
                ky += 1
            if i % 3 == 2:
                z_s[kz, j] = visited_states[i, j]
                kz += 1
    final_us = np.vstack((u_s[0, :], u_s[3, :], u_s[1, :], u_s[4, :], u_s[2, :], u_s[5, :]))
    return x_s, y_s, z_s, final_us


def update_scatter(i, traj, ax):
    ax.clear()
    ax.view_init(elev=12, azim=-6)
    # Setting the axes properties
    ax.set(xlim3d=(-1, 1), xlabel='X')
    ax.set(ylim3d=(0.5, 2.5), ylabel='Y')
    ax.set(zlim3d=(-1, 1), zlabel='Z')
    ax.scatter(traj[0, :, i].T, traj[1, :, i].T, traj[2, :, i].T)


def create_data_matrices(trajs, controls, indices):
    states = np.empty((n_states, 0))
    inputs = np.empty((n_inputs, 0))
    next_states = np.empty((n_states, 0))
    for i in indices:
        cloth_trajectory = trajs[i]
        cloth_inputs = controls[i]
        states = np.hstack((states, cloth_trajectory[:, :-1]))
        next_states = np.hstack((next_states, cloth_trajectory[:, 1:]))
        inputs = np.hstack((inputs, cloth_inputs[:, :-1]))
    augmented_states = np.vstack((states, inputs))
    X = augmented_states
    Y = next_states
    return X, Y


if __name__ == '__main__':
    ms = np.logspace(1.0, 2.6, num=20, dtype=int)  # np.arange(10, 200, 100)
    # ms = np.array([400])
    labels = ['nystrom', 'splines']
    all_trajs = []
    all_controls = []
    # swing_angle = 120
    path_to_experiment = pathlib.Path(f"./8x8_cloth_swing_xyz")
    n_states = 8 * 8 * 3
    n_inputs = 6
    n_trajs = 50
    n_training_trajs = 30
    validate_sys_id = False
    cross_validate = False
    n_val_trajs = 10  # for hyperparameter learning
    for i in range(0, n_trajs):
        traj = np.loadtxt(f'{path_to_experiment}/state_samples_cloth_swing_{i}.csv', delimiter=',').T
        controls = np.loadtxt(f'{path_to_experiment}/input_samples_cloth_swing_{i}.csv', delimiter=',')[:, :n_inputs].T
        all_trajs.append(traj)
        all_controls.append(controls)
    val_trajs = all_trajs[0:n_val_trajs]
    val_controls = all_controls[0:n_val_trajs]
    all_trajs = all_trajs[n_val_trajs:]
    all_controls = all_controls[n_val_trajs:]
    if cross_validate:
        for kapprox in labels:
            Xval, Yval = create_data_matrices(val_trajs, val_controls, range(0, n_val_trajs))
            best_clf = learn_hyperparams(Xval, Yval, kapprox)
            with open(f"{path_to_experiment}/cross_validated_kern_params_{kapprox}_cloth_swing.npy", 'wb') as f:
                pickle.dump(best_clf, f)
    if validate_sys_id:
        regressor = None
        for k, kapprox in enumerate(labels):
            all_rmse_across_seeds = np.empty((0, ms.shape[0]))
            print(kapprox)
            for seed in range(0, 20):
                np.random.seed(seed)  # Fix seed
                random.seed(seed)
                trajs_indices = np.arange(0, n_trajs - n_val_trajs)
                np.random.shuffle(trajs_indices)
                trainings_trajs_indices = trajs_indices[:n_training_trajs]
                testing_trajs_indices = trajs_indices[n_training_trajs:]
                test_trajs = [all_trajs[i] for i in testing_trajs_indices.tolist()]
                test_controlss = [all_controls[i] for i in testing_trajs_indices.tolist()]
                X, Y = create_data_matrices(all_trajs, all_controls, trainings_trajs_indices)
                with open(f"{path_to_experiment}/cross_validated_kern_params_{kapprox}_cloth_swing.npy", 'rb') as f:
                    clf = pickle.load(f)
                    kernel_params = clf.best_params_
                # kernel_params = {}
                print('----------------K PARAMS----------------', kernel_params)
                pathdata = pathlib.Path(f"{path_to_experiment}/sim_results/{kapprox}/data")
                pathplots = pathlib.Path(f"{path_to_experiment}/sim_results/{kapprox}/plots")
                pathlib.Path.mkdir(pathdata, parents=True, exist_ok=True)
                pathlib.Path.mkdir(pathplots, parents=True, exist_ok=True)
                all_rmses = []  # Collect RMSEs across all testing trajectories
                for i in range(0, len(test_trajs)):
                    test_traj = test_trajs[i]
                    test_controls = test_controlss[i]
                    rmses = []
                    for m_indx, m in enumerate(ms):
                        print("Curr seed and number of features: ", seed, " ", m)
                        kernel_params['m'] = m
                        regressor = None
                        if kapprox == 'nystrom':
                            regressor = KoopmanNystromRegressor(n_inputs, **kernel_params)
                        elif kapprox == 'splines':
                            regressor = KoopmanSplineRegressor(n_inputs, state_bounds_params=None, **kernel_params)
                        regressor.fit(X.T, Y.T)  # Careful with transposition (shape required by sklearn estimator API, used for CV)
                        # Evaluate the prediction accuracy in open loop
                        curr_rmse = validate_dyn_sys(regressor, test_traj, test_controls)
                        rmses.append(curr_rmse)
                    print(rmses)
                    all_rmses.append(rmses)
                all_rmses_array = np.array(all_rmses)
                all_rmse_across_seeds = np.vstack((all_rmse_across_seeds, all_rmses_array))
            print(np.median(all_rmse_across_seeds, axis=0))
            print(np.percentile(all_rmse_across_seeds, axis=0, q=15))
            print(np.percentile(all_rmse_across_seeds, axis=0, q=85))
            np.savetxt(f"{pathdata}/all_rmses_{kapprox}_cloth_swing_angle.csv", all_rmse_across_seeds)
    else:
        for seed in range(0, 50):
            for k, kapprox in enumerate(labels):
                print(kapprox)
                np.random.seed(seed)  # Fix seed
                random.seed(seed)
                trajs_indices = np.arange(0, n_trajs)
                trainings_trajs_indices = trajs_indices[:n_training_trajs]
                X, Y = create_data_matrices(all_trajs, all_controls, trainings_trajs_indices)

                with open(f"{path_to_experiment}/cross_validated_kern_params_{kapprox}_cloth_swing.npy", 'rb') as f:
                    clf = pickle.load(f)
                    kernel_params = clf.best_params_
                m = 100
                kernel_params['m'] = m

                regressor = None
                if kapprox == 'nystrom':
                    regressor = KoopmanNystromRegressor(n_inputs, **kernel_params)
                elif kapprox == 'splines':
                    regressor = KoopmanSplineRegressor(n_inputs, state_bounds_params=None, **kernel_params)
                regressor.fit(X.T, Y.T)
                A = regressor.A
                B = regressor.B
                C = regressor.C
                initial_state = all_trajs[0][:, 0].reshape([-1, 1])
                R = np.eye(n_inputs)
                Q = 0.075e-1 * C.T @ C
                Q = (Q + Q.T) / 2
                offset = np.zeros((n_states, 1))
                alpha = np.pi / 4
                z_top = initial_state[-1]
                vertical_shift = 0.0
                horizontal_shift = 0.0
                for i in range(0, n_states):
                    if i % 3 == 0:
                        offset[i] = 0
                    if i % 3 == 1:
                        r = abs(z_top - initial_state[i + 1])
                        offset[i] = r * np.sin(alpha) + horizontal_shift
                    if i % 3 == 2:
                        r = abs(z_top - initial_state[i])
                        offset[i] = r - r * np.cos(alpha) + vertical_shift
                reference = initial_state + offset
                pathdata = pathlib.Path(f"{path_to_experiment}/sim_results/{kapprox}/data")
                pathplots = pathlib.Path(f"{path_to_experiment}/sim_results/{kapprox}/plots")
                pathdata.mkdir(parents=True, exist_ok=True)
                pathplots.mkdir(parents=True, exist_ok=True)

                np.savetxt(f"{pathdata}/reference_lqr_m_{m}.csv", reference)
                K, _, _ = control.dlqr(A, B, Q, R)
                K_SOM = K[[0, 3, 1, 4, 2, 5], :]
                # Careful with actual input sequence used by F. Coltraro's simulator in matlab
                np.savetxt(f"{pathdata}/REBUTTAL_K_lqr_seed_{seed}_m_{m}.csv", K_SOM)
                with open(f"{pathdata}/REBUTTAL_regressor_seed_{seed}_m_{m}.npy", 'wb') as f:
                    pickle.dump(regressor, f)

                end = 60
                xs, ys, zs, us = lqr_control(end, reference, initial_state, regressor, K)
