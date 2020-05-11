# distutils: language = c++

from MC cimport MC as MCcpp
from libcpp.vector cimport vector
import numpy as np

cdef class MC:
    """

        Magnetic Metropolis Monte Carlo module. It does not have to work with pyiron, but for convenience
        it is strongly recommended to do so. Here's a short example:

            from pyiron import Project
            from mc import MC

            structure = Project('.').create_structure(element='Fe',
                                  bravais_basis='bcc',
                                  lattice_constant=2.85)
            structure.set_repeat(10)
            J = 0.1 # eV
            first_shell_tensor = structure.get_shell_matrix(1)

            mc = MC(len(structure))
            mc.set_heisenberg_coeff(J*first_shell_tensor)

            mc.run(temperature=300, number_of_iterations=1000)

        The results can be analysed by attributes like `get_mean_energy()` or `get_magnetic_moments()`.

        The Hamiltonian H is given in the form:

        H = -0.5*sum_ij J_ij*(m_i*m_j)^deg + sum_i A_i*m_i^deg

        The first coefficients (J_ij) can be inserted with set_heisenberg_coeff and the magnitude dependent
        terms (A_i) can be set via set_landau_coeff. Beware of the signs for the Heisenberg coefficients
    """
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

    def set_heisenberg_coeff(self, coeff, i=None, j=None, deg=1, index=0):
        """
            Args:
                coeff (float/list/ndarray): Heisenberg coefficient. If a single number is given,
                                            the same parameter is applied to all the pairs defined
                                            in me and neigh. Instead of giving me and neigh, you can
                                            also give a n_atom x n_atom tensor.
                me (list/ndarray): list of indices i (s. definition in the comment)
                neigh (list/ndarray): list of indices j (s. definition in the comment)
                deg (int): polynomial degree
                index (int): potential index for thermodynamic integration (0 or 1; choose 0 if
                             not thermodynamic integration)
            Comment:
                Heisenberg term is given by: -sum_ij coeff_ij*(m_i*m_j)^deg. Beware of the negative sign.
        """
        if i is None and j is None:
            n = self.c_mc.get_number_of_atoms()
            if np.array(coeff).shape!=(n, n):
                raise ValueError('If i and j are not specified, coeff has to be a 2d tensor with the length equal'
                                 'to the number of atoms in each direction.')
            i,j = np.where(coeff!=0)
            coeff = coeff[coeff!=0]
        coeff = np.array([coeff]).flatten()
        i = np.array(i).flatten()
        j = np.array(j).flatten()
        if len(coeff)==1:
            coeff = np.tile(coeff, len(i))
        if len(coeff)!=len(i) or len(i)!=len(j):
            raise ValueError('Length of vectors not the same')
        self.c_mc.set_heisenberg_coeff(coeff, i, j, deg, index)

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

    def get_magnetic_gradients(self):
        """
            Returns:
                nx3 array of magnetic gradients
        """
        m = self.c_mc.get_magnetic_gradients()
        return np.array(m).reshape(-1, 3)

    def get_magnetic_moments(self):
        """
            Returns:
                nx3 array of magnetic moments
        """
        m = self.c_mc.get_magnetic_moments()
        return np.array(m).reshape(-1, 3)

    def set_magnetic_moments(self, magmoms):
        """
            Args:
                magmoms (ndarray/list): nx3 magnetic moments to set
        """
        self.c_mc.set_magnetic_moments(np.array(magmoms).flatten())

    def get_average_magnetization(self, vector=False):
        """
            Returns:
                average magnetization value or vector
        """
        m = np.mean(self.get_magnetic_moments(), axis=0)
        if vector:
            return m
        return np.linalg.norm(m)

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

    def get_energy(self, index=0, per_atom=False):
        """
            Returns:
                Energy value of current snapshot in eV
        """
        if per_atom:
            return self.c_mc.get_energy(index)/self.c_mc.get_number_of_atoms()
        return self.c_mc.get_energy(index)

    def get_mean_energy(self, index=0, per_atom=False):
        """
            Returns:
                Mean energy value since the last reset (s. reset())
        """
        if per_atom:
            return self.c_mc.get_mean_energy(index)/self.c_mc.get_number_of_atoms()
        return self.c_mc.get_mean_energy(index)

    def get_energy_variance(self, index=0, per_atom=False):
        """
            Returns:
                Energy variance since the last reset (s. reset())
        """
        if per_atom:
            return self.c_mc.get_energy_variance(index)/self.c_mc.get_number_of_atoms()
        return self.c_mc.get_energy_variance(index)

    def revoke_qmc(self):
        self.set_eta(1)

    def set_eta(self, val):
        """
            Args:
                val (float): eta value for QMC
        """
        if val<0:
            raise ValueError("eta cannot be negative")
        self.c_mc.set_eta(val)

    def get_eta(self):
        """
            Returns:
                (float) current eta value
        """
        return self.c_mc.get_eta()

    def get_ground_state_energy(self, per_atom=False):
        """
            Args:
                per_atom (bool): return per atom energy instead of total energy

            Returns:
                energy value
        """
        if per_atom:
            return self.c_mc.get_ground_state_energy()/self.c_mc.get_number_of_atoms()
        return self.c_mc.get_ground_state_energy()

    def prepare_qmc(self, temperature, number_of_iterations=1, run_twice=True):
        """
            Args:
                temperature (float): Temperature in K
                number_of_iterations (int): Number of MC steps (internally multiplied
                                            by the number of atoms)
                run_twice (bool): run twice to make system converge in the first run
                                  and measure energy in the second run
        """
        if run_twice:
            self.set_eta(1)
            self.run(temperature, 1)
            self.set_eta(0)
            self.run(temperature, number_of_iterations)
        self.set_eta(0)
        self.run(temperature, number_of_iterations)
        self.set_eta((self.get_mean_energy(per_atom=True)
                      -self.get_ground_state_energy(per_atom=True))/(8.617e-5*temperature))

    def set_lambda(self, val):
        """
            Args:
                val (float): Lambda value for thermodynamic integration
        """
        if val<0 or val>1:
            raise ValueError("Lambda value has to be between 0 and 1")
        self.c_mc.set_lambda(val)

    def select_id(self, indices):
        """
            Args:
                indices (list): list of ID's that can be chosen
        """
        self.c_mc.select_ID(np.array(indices).tolist())

    def get_steps_per_second(self):
        """
            Returns:
                Number of MC steps performed per second
        """
        return self.c_mc.get_steps_per_second()

    def set_magnitude(self, dm, dphi, dtheta, flip=1):
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
        flip = np.array([flip]).flatten()
        if len(dm)==1:
            dm = np.array(self.c_mc.get_number_of_atoms()*dm.tolist())
        if len(dphi)==1:
            dphi = np.array(self.c_mc.get_number_of_atoms()*dphi.tolist())
        if len(dtheta)==1:
            dtheta = np.array(self.c_mc.get_number_of_atoms()*dtheta.tolist())
        if len(flip)==1:
            flip = np.array(self.c_mc.get_number_of_atoms()*flip.tolist())
        self.c_mc.set_magnitude(dm, dphi, dtheta, flip)

    def run_gradient_descent(self, max_iter=None, step_size=1, decrement=0.001, diff=1.0e-8):
        """
            args:
                max_iter (int): number of steps to perform
                step_size (float): step size
                decrement (float): decrement value
        """
        if max_iter is None:
            max_iter = self.c_mc.get_number_of_atoms()**2;
        return self.c_mc.run_gradient_descent(max_iter, step_size, decrement, diff)

    def run_debug(self):
        """
            run a few test calculations
        """
        self.c_mc.run_debug()

    def activate_debug(self):
        """
            Activate debug mode (not so helpful though...)
        """
        self.c_mc.activate_debug()
