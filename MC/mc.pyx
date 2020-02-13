# distutils: language = c++

from MC cimport MC as MCcpp
from libcpp.vector cimport vector
import numpy as np

cdef class MC:
    cdef MCcpp c_mc

    def __cinit__(self, number_of_atoms):
        self.c_mc.create_atoms(number_of_atoms)

    def set_landau_coeff(self, coeff, deg, index=0):
        coeff = np.array([coeff]).flatten()
        if len(coeff)==1:
            coeff = np.tile(coeff, self.c_mc.get_number_of_atoms())
        self.c_mc.set_landau_coeff(coeff, deg, index)

    def set_heisenberg_coeff(self, coeff, me, neigh, deg=1, index=0):
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
            dm = np.array(self.c_mc.get_number_of_atoms()*dm.tolist())
        if len(dphi)==1:
            dphi = np.array(self.c_mc.get_number_of_atoms()*dphi.tolist())
        if len(dtheta)==1:
            dtheta = np.array(self.c_mc.get_number_of_atoms()*dtheta.tolist())
        self.c_mc.set_magnitude(dm, dphi, dtheta)

