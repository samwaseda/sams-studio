# distutils: language = c++

from MC cimport MC as MCcpp
from libcpp.vector cimport vector
import numpy as np

cdef class MC:
    cdef MCcpp c_mc

    def __cinit__(self, structure):
        n_tot = len(structure)
        A = structure.A
        B = structure.B
        J = np.array(structure.J)
        num_neigh = J.shape[1]
        neigh = np.array(structure.neigh)
        me = np.arange(len(structure))[:, np.newaxis]*np.ones_like(neigh)
        me = me[J!=None].flatten()
        neigh = neigh[J!=None].flatten()
        J = J[J!=None].flatten()
        self.c_mc.create_atoms(n_tot, num_neigh, A, B, me, neigh, J)
