from libcpp.vector cimport vector

cdef extern from "MC.cpp":
    pass

cdef extern from "MC.h":
    cdef cppclass MC:
        MC() except +
        MC(double, bool) except +
        void create_atoms(int, int, vector[double], vector[double], vector[int], vector[int], vector[double])
        double run(double)
        double output(string, bool)
