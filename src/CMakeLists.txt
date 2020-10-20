
add_library(src simul.c monte_carlo.c ising_lattice.c util.c media.c)
set_property(TARGET src PROPERTY C_STANDARD 99)
