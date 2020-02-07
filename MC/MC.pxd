from libcpp.vector cimport vector

cdef extern from "MC.cpp":
    pass

cdef extern from "MC.h":
    cdef cppclass MC:
        MC() except +
        void create_atoms(vector[double], vector[double], vector[int], vector[int], vector[double]) except +
        void append_parameters(vector[double], vector[double], vector[int], vector[int], vector[double])
        void set_landau_coeff(vector[double], int, int) except +
        void set_heisenberg_coeff(vector[double], vector[int], vector[int], int, int) except +
        void run(double, int) except +
        vector[double] get_magnetic_moments()
        double get_acceptance_ratio()
        double get_energy(int)
        double get_energy_variance(int)
        double get_mean_energy(int)
        double set_lambda(double) except +
        double get_steps_per_second()
        double set_magnitude(double, double, double)
        void reset()
