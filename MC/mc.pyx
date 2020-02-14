# distutils: language = c++

from MC cimport MC as MCcpp
from libcpp.vector cimport vector
import numpy as np

cdef class MC:
    cdef MCcpp c_mc

    def __cinit__(self, number_of_atoms):
        """
        args:
            number_of_atoms (int): number of atoms
        """
        self.c_mc.create_atoms(number_of_atoms)

    def set_landau_coeff(self, coeff, deg, index=0):
        """
            Args:
                coeff (float/list/ndarray): Coefficient to the Landau term. If a single number
                                            is given, the same parameter is applied to all the
                                            atoms.
                deg (int): Polynomial degree (usually an even number)
                index (int): Potential index for thermodynamic integration (0 or 1; choose 0 if
                             not thermodynamic integration)

            Comment:
                Landau term is given by: sum_i coeff_i*m_i^deg
        """
        coeff = np.array([coeff]).flatten()
        if len(coeff)==1:
            coeff = np.tile(coeff, self.c_mc.get_number_of_atoms())
        if len(coeff)!=self.c_mc.get_number_of_atoms():
            raise ValueError('coeff has to be a single number or a list of length'
                             'corresponding to the number of atoms in the box')
        self.c_mc.set_landau_coeff(coeff, deg, index)

    def set_heisenberg_coeff(self, coeff, me, neigh, deg=1, index=0):
        """
            Args:
                coeff (float/list/ndarray): Heisenberg coefficient. If a single number is given,
                                            the same parameter is applied to all the pairs defined
                                            in me and neigh.
                me (list/ndarray): list of indices i (s. definition in the comment)
                neigh (list/ndarray): list of indices j (s. definition in the comment)
                deg (int): polynomial degree
                index (int): potential index for thermodynamic integration (0 or 1; choose 0 if
                             not thermodynamic integration)
            Comment:
                Heisenberg term is given by: -sum_ij coeff_ij*(m_i*m_j)^deg. Beware of the negative sign.
        """
        coeff = np.array([coeff]).flatten()
        me = np.array(me).flatten()
        neigh = np.array(neigh).flatten()
        if len(coeff)==1:
            coeff = np.tile(coeff, len(me))
        if len(coeff)!=len(me) or len(me)!=len(neigh):
            raise ValueError('Length of vectors not the same')
        self.c_mc.set_heisenberg_coeff(coeff, me, neigh, deg, index)

    def clear_heisenberg_coeff(self, index=0):
        """
            Args:
                index (int): potential index for thermodynamic integration (0 or 1; choose 0 if
                             not thermodynamic integration)

            This function erases all the Heisenberg coefficients defined before.
        """
        self.c_mc.clear_heisenberg_coeff(index)

    def clear_landau_coeff(self, index=0):
        """
            Args:
                index (int): potential index for thermodynamic integration (0 or 1; choose 0 if
                             not thermodynamic integration)

            This function erases all the Landau coefficients defined before.
        """
        self.c_mc.clear_landau_coeff(index)

    def clear_all_coeff(self):
        """
            This function erases all the coefficients defined before.
        """
        self.clear_heisenberg_coeff()
        self.clear_heisenberg_coeff(1)
        self.clear_landau_coeff()
        self.clear_landau_coeff(1)

    def get_magnetic_moments(self):
        """
            Returns:
                nx3 array of magnetic moments
        """
        m = self.c_mc.get_magnetic_moments()
        return np.array(m).reshape(-1, 3)

    def run(self, temperature, number_of_iterations=1, reset=True):
        """
            Args:
                temperature (float): Temperature in K
                number_of_iterations (int): Number of MC steps (internally multiplied
                                            by the number of atoms)
                reset (bool): Resets statistics (s. get_mean_energy() etc.)
        """
        if reset:
            self.c_mc.reset()
        self.c_mc.run(temperature, number_of_iterations)

    def get_acceptance_ratio(self, individual=False):
        """
            Args:
                individual (bool): If true, an array containing each acceptance ratio is returned
            Returns:
                Acceptance ratio since the last reset (s. reset())
        """
        if individual:
            return np.array(self.c_mc.get_acceptance_ratios())
        return self.c_mc.get_acceptance_ratio()

    def get_energy(self, index=0):
        """
            Returns:
                Energy value of current snapshot in eV
        """
        return self.c_mc.get_energy(index)

    def get_mean_energy(self, index=0):
        """
            Returns:
                Mean energy value since the last reset (s. reset())
        """
        return self.c_mc.get_mean_energy(index)

    def get_energy_variance(self, index=0):
        """
            Returns:
                Energy variance since the last reset (s. reset())
        """
        return self.c_mc.get_energy_variance(index)

    def set_lambda(self, val):
        """
            Args:
                val (float): Lambda value for thermodynamic integration
        """
        if val<0 or val>1:
            raise ValueError("Lambda value has to be between 0 and 1")
        self.c_mc.set_lambda(val)

    def get_steps_per_second(self):
        """
            Returns:
                Number of MC steps performed per second
        """
        return self.c_mc.get_steps_per_second()

    def set_magnitude(self, dm, dphi, dtheta):
        """
            Args:
                dm (float/list/ndarray): Magnitude variation strength. If a single value is set, the same
                                         value is used for all the atoms (applies to dphi and dtheta)
                dphi (float/list/ndarray): Phi variation strength
                dtheta (float/list/ndarray): Theta variation strength

            Comment:
                MC algorithm proposes a new magnetic moment by

                m_new = m_old+dm*random_number # it is ensured afterwards that m_new won't be negative
                theta_new = acos(cos(theta_old)+dtheta*random_number); # it is ensured internally that the
                                                                         term inside acos is between -1 and 1
                phi_new = phi_old+dphi*2*PI*random_number;

                where random_number is a random number between -1 and 1.
        """
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

    def activate_debug(self):
        """
            Activate debug mode (not so helpful though...)
        """
        self.c_mc.activate_debug()
