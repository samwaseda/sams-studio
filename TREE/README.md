# Tree structure

This code is meant for an efficient calculation of kinetic Monte Carlo method

## How to compile

Download all files and run `python setup.py build_ext --inplace`.

## How to use it

Suppose you have $N$ diffusing solutes (such as vacancies or interstitial atoms). Let's say each of them has 4 possible jumps. Then you can create the tree by:

```
from tree import Tree
import numpy as np

tree = Tree()
tree.append([[0.1, 0.2, 0.3, 0.4]
             [0.5, 0.6, 0.7, 0.8], # jump probabilities
	    [0, 1] # indices
	   )
```

Then you can put a random number between 0 and 1:

```
tree.choose_event(0.5)
tree.get_index() # This will return 1
tree.get_jump_id() # This will return 2

# Then either remove the leaf by:
tree.remove()

# Or update the values for example by:
tree.update([0.9, 1, 1.1, 1.3])

```
