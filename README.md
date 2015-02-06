**Quasiclass** -- simulation of proton-coupled electron transfer with a presence of bath.

The code randomizes a large number of initial trajectories guided by coherent state and temperature distributions. The trajectories are propagated using Hamilton's equations of motion. The average quantity in the end shows an overall dynamics with presence of quantum effects.

`source` directory contains `main.cpp`, which is the starting point.

The biggest challenges were nummerical stability at longer times and optimizing parameters to get correct quantum effects while minimizing the number of trajectories.

To increase redability and performance while reducing bugs, the code's style involves use of std/boost:: algorithms over for-loops, functors over large classes, tuples, and boost/std numerical libraries for ODE and distributions.

