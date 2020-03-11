# Magnetic Metropolis Monte Carlo following classical Boltzmann statistics

This code allows you to launch Metropolis Monte Carlo simulations via Heisenberg Landau models (with various polynomial degrees) from a jupyter notebook.

## How to compile

Download all files and run `python setup.py build_ext --inplace`.

## How to launch a simple calculation

Suppose we have a bcc Fe system created by [pyiron](http://github.com/pyiron/pyiron) (e.g. via `structure = Project('.').create_structure('Fe', 'bcc', 2.83).repeat(10)`) and a Heisenberg coefficient `J=0.1` (eV). Then the magnetic interactions can be calculated by:

```
from pyiron import Project
from mc import MC

structure = Project('.').create_structure('Fe', 'bcc', 2.85)
J = 0.1 # eV
first_shell_tensor = structure.get_shell_matrix(1)

mc = MC(len(structure))
mc.set_heisenberg_coeff(J*first_shell_tensor)

mc.run(temperature=300, number_of_iterations=1000)
```

The results can be analysed by attributes like `get_mean_energy()` or `get_magnetic_moments`.

For more info, take a look at `mc.pyx`.
