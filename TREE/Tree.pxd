from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "Tree.cpp":
    pass

cdef extern from "Tree.h":
    cdef cppclass Tree:
        Tree() except +
        void append(vector[double], int) except +
        int get_index() except +
        int get_jump_id() except +
        void remove() except +
        void choose_event(double) except +
        double get_kappa() except +
        string get_structure() except +
        void update_kappa(vector[double]) except +
