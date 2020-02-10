# distutils: language = c++

from MC cimport MC as MCcpp
from libcpp.vector cimport vector
import numpy as np

cdef class MC:
    cdef MCcpp c_mc
    cdef readonly int number_of_atoms

    def __cinit__(self, structure):
        self.number_of_atoms = len(structure)
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
            self.c_mc.append_parameters(A[1], B[0], me[1].flatten()[J[1].flatten()!=0], neigh[1].flatten()[J[1].flatten()!=0], J[1].flatten()[J[1].flatten()!=0])
        else:
            self.c_mc.create_atoms(A, B ,me.flatten()[J.flatten()!=0], neigh.flatten()[J.flatten()!=0], J.flatten()[J.flatten()!=0])

    def set_landau_coeff(self, coeff, deg, index=0):
        coeff = np.array([coeff]).flatten()
        if len(coeff)==1:
            coeff = np.tile(coeff, self.number_of_atoms)
        self.c_mc.set_landau_coeff(coeff, deg, index)

    def set_heisenberg_coeff(self, coeff, me, neigh, deg, index=0):
        coeff = np.array([coeff]).flatten()
        me = np.array(me).flatten()
        neigh = np.array(neigh).flatten()
        if len(coeff)==1:
            coeff = np.tile(coeff, len(me))
        if len(coeff)!=len(me) or len(me)!=len(neigh):
            raise ValueError('Length of vectors not the same')
        self.c_mc.set_heisenberg_coeff(coeff, me, neigh, deg, index)

    def clear_heisenberg_coeff(self, index=0):
        self.c_mc.clear_heisenberg_coeff(index)

    def clear_landau_coeff(self, index=0):
        self.c_mc.clear_landau_coeff(index)

    def clear_all_coeff(self):
        self.clear_heisenberg_coeff()
        self.clear_heisenberg_coeff(1)
        self.clear_landau_coeff()
        self.clear_landau_coeff(1)

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
        dm = np.array([dm]).flatten()
        dphi = np.array([dphi]).flatten()
        dtheta = np.array([dtheta]).flatten()
        if len(dm)==1:
            dm = np.array(self.number_of_atoms*dm.tolist())
        if len(dphi)==1:
            dphi = np.array(self.number_of_atoms*dphi.tolist())
        if len(dtheta)==1:
            dtheta = np.array(self.number_of_atoms*dtheta.tolist())
        self.c_mc.set_magnitude(dm, dphi, dtheta)

