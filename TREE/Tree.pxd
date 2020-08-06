from libcpp.vector cimport vector

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
