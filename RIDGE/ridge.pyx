# distutils: language = c++

from Ridge cimport Ridge as Ridgecpp
from libcpp.vector cimport vector
import numpy as np

cdef class Ridge:
    """
        Ridge regression class. Initiate the class via:

        Ridge(x, y)
        
    """
    cdef Ridgecpp c_ridge

    def __cinit__(self, x, y, zeroth=1, true_error=False):
        y = np.array(y).flatten()
        if len(x)!=len(y):
            raise ValueError('Values are not consistent')
        self.c_ridge.initialize_sets(x.flatten(), y.flatten(), zeroth, true_error)

    def set_number_of_cv_set(self, n_set):
        self.c_ridge.set_number_of_cv_set(n_set)

    def run_least_square(self):
        self.c_ridge.least_square()

    def run_random(self, max_step):
        self.c_ridge.random(max_step)

    def run_raw(self, max_step):
        self.c_ridge.raw(max_step)

    def run_classic(self, max_step):
        self.c_ridge.classic(max_step)

    def run_conjugate_gradient(self, max_step):
        self.c_ridge.conjugate_gradient(max_step)

    def run_gradient_descent(self, max_step, prefactor, max_descent=0):
        """
            Args:
                max_step (int): Number of steps
                prefactor (float): Prefactor to be multiplied with the gradient
                max_descent (float): Maximum displacement at each step
        """
        self.c_ridge.gradient_descent(max_step, prefactor, max_descent)

    def run_fixed_descent(self, max_step, displacement):
        self.c_ridge.fixed_descent(max_step, displacement)

    def get_coeff(self):
        return np.array(self.c_ridge.get_coeff())

    def get_lambda(self):
        return np.array(self.c_ridge.get_lambda())

    def activate_debug(self):
        self.c_ridge.activate_debug()

    def get_validation_error(self):
        return self.c_ridge.get_validation_error()

    def get_true_error(self):
        return self.c_ridge.get_true_error()

    def get_determinant(self, regularization=True):
        return np.array(self.c_ridge.get_determinant(regularization))

    def get_sum_H(self):
        return self.c_ridge.get_sum_H()

    def get_derivative(self):
        return np.array(self.c_ridge.get_derivative())

    def get_hessian(self):
        H = np.array(self.c_ridge.get_hessian())
        n = int(np.sqrt(len(H)))
        return H.reshape(n, n)

    def set_lambda(self, lambda_in):
        self.c_ridge.set_lambda(lambda_in)

    def reset_increment(self):
        self.c_ridge.reset_increment()

    def set_min_lambda(self, min_value=-10):
        L = self.get_lambda()
        L[L<min_value] = min_value
        self.c_ridge.set_lambda(L)

    def set_max_lambda(self, min_value=10):
        L = self.get_lambda()
        L[L>min_value] = min_value
        self.c_ridge.set_lambda(L)
