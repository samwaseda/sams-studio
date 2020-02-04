from libcpp.vector cimport vector

cdef extern from "MC.cpp":
    pass

cdef extern from "MC.h":
    cdef cppclass MC:
        MC() except +
        void create_atoms(int, vector[double], vector[double], vector[int], vector[int], vector[double])
        double run(double, int) except +
        vector[double] get_magnetic_moments()
        double get_acceptance_ratio()
        double get_energy()
        double get_energy_variance()
        double get_mean_energy()
        double set_lambda(double) except +
        double get_steps_per_second()
        double set_magnitude(double, double, double)
        void reset()
