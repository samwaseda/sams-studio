# distutils: language = c++

from MC cimport MC as MCcpp
from libcpp.vector cimport vector
import numpy as np

cdef class MC:
    cdef MCcpp c_mc

    def __cinit__(self, structure):
        A = np.array(structure.A)
        thermodynamic_integration = False
        if len(A.shape)==2:
            thermodynamic_integration = True
        B = np.array(structure.B)
        J = np.array(structure.J)
        neigh = np.array(structure.neigh)
        if thermodynamic_integration:
            me = np.arange(len(structure))[np.newaxis, :, np.newaxis]*np.ones_like(neigh)
            if len(A[0])!=len(structure) or len(B[0])!=len(structure):
                raise ValueError('Length of A or B is not the same as the structure length')
        else:
            me = np.arange(len(structure))[:, np.newaxis]*np.ones_like(neigh)
            if len(A)!=len(structure) or len(B)!=len(structure):
                raise ValueError('Length of A or B is not the same as the structure length')
        if np.min(B)<0:
            raise ValueError('Negative B value will make the simulation explode')
        if np.max(neigh)>=len(structure) or np.min(neigh)<0:
            raise AssertionError('Problem with the neighbor ids')
        if np.max(me)>=len(structure) or np.min(me)<0:
            raise AssertionError('Problem with the ids')
        if thermodynamic_integration:
            self.c_mc.create_atoms(A[0], B[0], me[0].flatten()[J[0].flatten()!=0], neigh[0].flatten()[J[0].flatten()!=0], J[0].flatten()[J[0].flatten()!=0])
            self.c_mc.append_parameters(A[1], B[1], me[1].flatten()[J[1].flatten()!=0], neigh[1].flatten()[J[1].flatten()!=0], J[1].flatten()[J[1].flatten()!=0])
        else:
            self.c_mc.create_atoms(A, B, me.flatten()[J.flatten()!=0], neigh.flatten()[J.flatten()!=0], J.flatten()[J.flatten()!=0])

    def get_magnetic_moments(self):
        m = self.c_mc.get_magnetic_moments()
        return np.array(m).reshape(-1, 3)

    def run(self, temperature, number_of_iterations=1, reset=True):
        if reset:
            self.c_mc.reset()
        self.c_mc.run(temperature, number_of_iterations)

    def get_acceptance_ratio(self):
        return self.c_mc.get_acceptance_ratio()

    def get_energy(self, index=0):
        return self.c_mc.get_energy(index)

    def get_mean_energy(self, index=0):
        return self.c_mc.get_mean_energy(index)

    def get_energy_variance(self, index=0):
        return self.c_mc.get_energy_variance(index)

    def set_lambda(self, val):
        self.c_mc.set_lambda(val)

    def get_steps_per_second(self):
        return self.c_mc.get_steps_per_second()

    def set_magnitude(self, dm, dphi, dtheta):
        self.c_mc.set_magnitude(dm, dphi, dtheta)

