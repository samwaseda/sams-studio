from libcpp.vector cimport vector

cdef extern from "MC.cpp":
    pass

cdef extern from "MC.h":
    cdef cppclass Energy:
        Energy() except +
        Energy(double, bool) except +
        void create_atoms(int, int, vector[double], vector[double])
        double MC(double)
        double output(string, bool)
