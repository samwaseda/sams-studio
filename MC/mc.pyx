#distutils: language = c++

from MC cimport Energy
from libcpp.vector import vector

cdef class PyMC:
    cdef Energy c_mc

    def __cinit__(self, structure):
        n_tot = len(structure)
        num_neigh = 58
        A = structure.A
        B = structure.B
        self.c_mc.create_atoms(n_tot, num_neigh, A, B)
