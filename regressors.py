'''
Author: Edoardo Caldarelli
Affiliation: Institut de Robòtica i Informàtica Industrial, CSIC-UPC
email: ecaldarelli@iri.upc.edu
January 2024
'''


import scipy.linalg
import scipy.signal
from sklearn.gaussian_process.kernels import RBF, Matern
from sklearn.base import BaseEstimator
import numpy as np

class ThreeDimensionalKernel():
    def __init__(self, lx, ly, lz, n_states):
        l = [lx, ly, lz]
        all_ls = np.zeros((1, n_states))
        for i in range(0, all_ls.shape[1]):
            j = i % 3
            all_ls[:, i] = l[j]
        self.kernel = RBF(all_ls)

class KernelWrapper():
    def __init__(self, ls):
        self.kernel = Matern(ls, nu=2.5)

class KoopmanRegressor(BaseEstimator):
    def __init__(self, n_inputs, gamma, m):
        self.gamma = gamma
        self.m = m
        self.A = None
        self.B = None
        self.C = None
        self.weights = None
        self.n_inputs = n_inputs

    def lift(self, X):
        raise NotImplementedError

    def fit(self, X, Y):
        raise NotImplementedError

    def predict(self, X_aug):
        # Method takes the augmented state, this is transposed
        n_states = X_aug.shape[1] - self.n_inputs
        X = X_aug.T
        X_lifted = self.lift(X[:n_states, :])
        # X_lifted = np.vstack((X[:n_states, :], X_lifted))
        phi_X = np.vstack((X_lifted, X[n_states:, :]))
        return (self.weights @ phi_X).T


class KoopmanNystromRegressor(KoopmanRegressor):
    def __init__(self, n_inputs, kernel, gamma=None, m=None):
        super().__init__(n_inputs, gamma, m)
        self.kernel = kernel
        self.nystrom_centers_input = None
        self.nystrom_centers_output = None
        self.jitter = 1e-6

    def fit(self, X, Y):
        X = X.T
        Y = Y.T

        n_states = X.shape[0] - self.n_inputs
        gamma_n = self.gamma * X.shape[1]

        if self.nystrom_centers_output is None:
            nystrom_centers_indices = np.random.choice(np.arange(0, Y.shape[1]), size=self.m,
                                                       replace=False)
            self.nystrom_centers_output = Y[:, nystrom_centers_indices]
        if self.nystrom_centers_input is None:
            self.nystrom_centers_input = self.nystrom_centers_output

        kern = self.kernel

        # Build all kernel matrices: m means Nys. landmarks, n means training points, x means state only
        K_mm_out = kern.kernel(self.nystrom_centers_output.T, self.nystrom_centers_output.T) + self.jitter * np.eye(self.m)
        K_mm_out_sqrt = scipy.linalg.sqrtm(K_mm_out).real
        K_mn_out = kern.kernel(self.nystrom_centers_output.T, Y.T)
        K_mn_in_x = kern.kernel(self.nystrom_centers_input.T, X[:n_states, :].T)
        K_mm_in_x = kern.kernel(self.nystrom_centers_input.T, self.nystrom_centers_input.T) + self.jitter * np.eye(self.m)
        K_mm_in_x_out = kern.kernel(self.nystrom_centers_input.T, self.nystrom_centers_output.T)

        # Augment to account for control with linear kernel (S_u is u, proj. is identity)
        K_mn_in = np.vstack((K_mn_in_x, X[n_states:, :]))
        K_mm_in = scipy.linalg.block_diag(K_mm_in_x, np.eye(self.n_inputs))

        # Invert matrix to compute "full" operator G in the notes
        inner_term = K_mn_in @ K_mn_in.T + gamma_n * K_mm_in
        right_term = scipy.linalg.block_diag(scipy.linalg.solve(K_mm_out_sqrt, K_mm_in_x_out.T, assume_a='her').T, np.eye(self.n_inputs))
        left_term = scipy.linalg.solve(K_mm_out_sqrt, (K_mn_out @ K_mn_in.T), assume_a='her')

        sol = scipy.linalg.lstsq(inner_term, right_term)[0]
        G_ls = left_term @ sol

        self.A = G_ls[:, : self.m]
        self.B = G_ls[:, self.m:]

        # Kernel ridge regression to reconstruct the state from the features (output space)
        inner_term_rec = gamma_n * K_mm_out + K_mn_out @ K_mn_out.T
        right_term_rec = scipy.linalg.sqrtm(K_mm_out).real
        left_term_rec = Y @ K_mn_out.T
        sol_rec = scipy.linalg.lstsq(inner_term_rec, right_term_rec)[0]
        self.C = left_term_rec @ sol_rec
        W = self.C @ G_ls

        self.weights = W

    def lift(self, X):
        # X is not augmented, no input
        kern = self.kernel
        Kmm = kern.kernel(self.nystrom_centers_output.T, self.nystrom_centers_output.T) + self.jitter * np.eye(self.m)
        Kmm_sqrt = scipy.linalg.sqrtm(Kmm).real
        Kmn = kern.kernel(self.nystrom_centers_output.T, X.T)
        phi = scipy.linalg.solve(Kmm_sqrt, Kmn, assume_a='her')
        return phi


class KoopmanSplineRegressor(KoopmanRegressor):
    def __init__(self, n_inputs, state_bounds_params=None, m=None, gamma=None):
        super().__init__(n_inputs, gamma, m)
        self.state_bounds_params = state_bounds_params
        self.centers = None

    def compute_centers(self, X):
        # n_states = self.state_bounds_params.shape[1]
        if self.state_bounds_params is not None:
            length = np.sqrt(np.random.uniform(0, self.state_bounds_params[0], size=(1, self.m)))
            angle = np.pi * np.random.uniform(0, self.state_bounds_params[1], size=(1, self.m))
            centers = np.vstack((length * np.cos(angle), length * np.sin(angle)))
            return centers
        else:
            centers_indices = np.random.choice(np.arange(0, X.shape[1]), size=self.m,
                                               replace=False)
            return X[:, centers_indices]

    def fit(self, X, Y):
        # Exact same code as in Korda and Mezic
        X = X.T
        Y = Y.T
        gamma_n = self.gamma * X.shape[1]

        n_states = X.shape[0] - self.n_inputs

        phi_x = self.lift(X[:n_states, :])
        phi_y = self.lift(Y)
        phi_x = np.vstack((phi_x, X[n_states:, :]))  # Stack control
        phi_y = np.vstack((phi_y, X[:n_states, :]))  # Learn optimal linear state reconstruction

        cov = phi_x @ phi_x.T
        cross_cov = phi_y @ phi_x.T
        M_ls = cross_cov @ scipy.linalg.pinv(cov + gamma_n * np.eye(cov.shape[0]))
        W = M_ls[self.m:, :self.m] @ M_ls[:self.m, :]

        self.A = M_ls[:self.m, :self.m]
        self.B = M_ls[:self.m, self.m:]
        self.C = M_ls[self.m:, :self.m]

        self.weights = W

    def lift(self, X):
        # Uses same set of centers for inputs and outputs
        if self.centers is None:
            self.centers = self.compute_centers(X)
        Cbig = self.centers
        phi = np.zeros((self.m, X.shape[1]))
        for i in range(0, Cbig.shape[1]):
            Ccurr = Cbig[:, i].reshape([-1, 1])
            r_squared_x = np.sum(np.square(X - Ccurr), axis=0, keepdims=True)
            phi[i, :] = np.multiply(r_squared_x, np.log(np.sqrt(r_squared_x)))
        phi = np.nan_to_num(phi, nan=0.0)
        return phi