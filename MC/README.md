# Magnetic Metropolis Monte Carlo following classical Boltzmann statistics

This code allows you to launch Metropolis Monte Carlo simulations via Heisenberg Landau models (with various polynomial degrees) from a jupyter notebook.

## How to compile

Download all files and run `python setup.py build_ext --inplace`.

## How to launch a simple calculation

Suppose we have a bcc Fe system created by [pyiron](http://github.com/pyiron/pyiron) (e.g. via `structure = Project('.').create_structure('Fe', 'bcc', 2.83).repeat(10)`) and a Heisenberg coefficient `J=0.1` (eV). Then the magnetic interactions can be calculated by:

```
import numpy as np
from mc import MC

J = 0.1
neighbors = structure.get_neighbors(num_neighbors=8) # 8 NN atoms
neighbor_indices = neighbors.indices
my_indices = np.arange(len(structure))[:, np.newaxis]*np.ones_like(neighbor_indices)

mc = MC(len(structure))
mc.set_heisenberg_coeff(J, my_indices, neighbor_indices)

mc.run(temperature=300, number_of_iterations=1000)
```

The results can be analysed by attributes like `get_mean_energy()` or `get_magnetic_moments`.

For more info, take a look at `mc.pyx`.
