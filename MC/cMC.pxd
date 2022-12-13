from libcpp.vector cimport vector

cdef extern from "cMC.cpp":
    pass

cdef extern from "cMC.h":
    cdef cppclass cMC:
        cMC() except +
        void create_atoms(int) except +
        void set_landau_coeff(vector[double], int, int) except +
        void set_heisenberg_coeff(vector[double], vector[int], vector[int], int, int) except +
        void clear_landau_coeff(int) except +
        void clear_heisenberg_coeff(int) except +
        void run(double, int, int) except +
        void activate_debug()
        vector[double] get_magnetic_moments()
        vector[double] get_magnetic_gradients()
        void set_magnetic_moments(vector[double]) except +
        double get_acceptance_ratio()
        vector[double] get_acceptance_ratios()
        double get_energy(int)
        double get_energy_variance(int)
        double get_mean_energy(int)
        double set_lambda(double) except +
        double get_steps_per_second()
        double set_magnitude(vector[double], vector[double], vector[int]) except +
        int get_number_of_atoms()
        double run_gradient_descent(int, double, double, double) except +
        void select_id(vector[int]) except +
        void set_metadynamics(double, double, double, int, double, int) except +
        void reset()
        void switch_spin_dynamics(bool, double, double, bool) except +
        vector[double] get_magnetization() except +
        vector[double] get_histogram(int) except +

