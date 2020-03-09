from libcpp.vector cimport vector

cdef extern from "MC.cpp":
    pass

cdef extern from "MC.h":
    cdef cppclass MC:
        MC() except +
        void create_atoms(int) except +
        void set_landau_coeff(vector[double], int, int) except +
        void set_heisenberg_coeff(vector[double], vector[int], vector[int], int, int) except +
        void clear_landau_coeff(int) except +
        void clear_heisenberg_coeff(int) except +
        void run(double, int) except +
        void activate_debug()
        vector[double] get_magnetic_moments()
        void set_magnetic_moments(vector[double]) except +
        double get_acceptance_ratio()
        vector[double] get_acceptance_ratios()
        double get_energy(int)
        double get_energy_variance(int)
        double get_mean_energy(int)
        double set_lambda(double) except +
        double get_eta() except +
        void set_eta(double) except +
        void prepare_qmc(double, int) except +
        double get_steps_per_second()
        double set_magnitude(vector[double], vector[double], vector[double]) except +
        int get_number_of_atoms()
        double run_gradient_descent(int, double, double) except +
        void reset()
