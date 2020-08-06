# distutils: language = c++

from Tree cimport Tree as Treecpp
from libcpp.vector cimport vector
import numpy as np

cdef class Tree:
    """
        Tree class. Initiate the class via:
    """
    cdef Treecpp c_tree

    def __cinit__(self):
        self.number_of_atoms = None
        self.number_of_possible_jumps = 4

    def choose_event(self, kappa):
        self.c_tree.choose_event(kappa)

    def append(self, kappa, indices):
        kappa = np.array(kappa).reshape(-1, self.number_of_possible_jumps)
        if len(kappa) != len(indices):
            raise ValueError('Length of kappa must be the same as that of indices')

    def remove(self):
        self.c_tree.remove()

    def get_index(self):
        return self.c_tree.get_index()

    def get_jump_id(self):
        return self.c_tree.get_jump_id()
