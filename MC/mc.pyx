# distutils: language = c++

from MC cimport MC as MCcpp
from libcpp.vector cimport vector
import numpy as np

cdef class MC:
    cdef MCcpp c_mc

    def __cinit__(self, structure):
        A = structure.A
        B = structure.B
        J = np.array(structure.J)
        num_neigh = J.shape[1]
        neigh = np.array(structure.neigh)
        me = np.arange(len(structure))[:, np.newaxis]*np.ones_like(neigh)
        me = me[J!=None].flatten()
        neigh = neigh[J!=None].flatten()
        J = J[J!=None].flatten()
        if len(A)!=len(structure) or len(B)!=len(structure):
            raise ValueError('Length of A or B is not the same as the structure length')
        if len(J)>len(structure)*num_neigh:
            raise ValueError('Problem with the structure length and/or number of neighbors')
        if np.min(B)<0:
            raise ValueError('Negative B value will make the simulation explode')
        if np.max(neigh)>=len(structure) or np.min(neigh)<0:
            raise AssertionError('Problem with the neighbor ids')
        if np.max(me)>=len(structure) or np.min(me)<0:
            raise AssertionError('Problem with the ids')
        self.c_mc.create_atoms(num_neigh, A, B, me, neigh, J)

    def run(self, temperature, number_of_iterations=1):
        self.c_mc.run(temperature, number_of_iterations)

    def output(self, configuration=True, reset=True):
        _ = self.c_mc.output(configuration)
        if reset:
            self.c_mc.reset()
