from libcpp.vector cimport vector

cdef extern from "Ridge.cpp":
    pass

cdef extern from "Ridge.h":
    cdef cppclass Ridge:
        Ridge() except +
        void activate_debug()
        void set_number_of_cv_set(int) except +
        void initialize_sets(vector[double], vector[double], int, int) except +
        void ridge(double) except +
        double get_validation_error() except +
        double get_true_error() except +
        void run_classic(int) except +
        void chi_training(bool) except +
        void least_square() except +
        void random(int) except +
        void raw(int) except +
        void classic(int) except +
        void reset_increment() except +
        void set_lambda(vector[double]) except +
        void gradient_descent(int, double, double) except +
        void fixed_descent(int, double) except +
        void conjugate_gradient(int) except +
        vector[double] get_coeff()
        vector[double] get_lambda()
        vector[double] get_determinant(int) except +
        vector[double] get_sum_H() except +
        vector[double] get_derivative() except +
        vector[double] get_hessian() except +
