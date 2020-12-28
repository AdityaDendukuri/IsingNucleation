
add_library(src transition_path_sampl.c simul.c monte_carlo.c ising_lattice.c util.c media.c)
set_property(TARGET src PROPERTY C_STANDARD 99)
