

|--------------------------------------------------------------------|
|                        INTRODUCTION                                |
|--------------------------------------------------------------------|
|This is a Monte Carlo simulation program for the Ising spin model   |
|written by Aditya Dendukuri. Simulations include free energy        |
|construction [1], nucleus formation [2] with umbrella sampling [1]. |
|--------------------------------------------------------------------|

|--------------------------------------------------------------------|
|                   INSTRUCTIONS TO RUN                              |
|--------------------------------------------------------------------|
|1. Create the directories build/ and output/                        |
|          >> mkdir build output                                     |
|                                                                    |
|2. Go to build directory and run cmake as follows                   |
|          >> cd build                                               |
|          >> cmake -C ../ -S .                                      |
|                                                                    |
|3. This will generate the makefile. Compile the program by simply   |
|   typing:                                                          |
|          >> make                                                   |
|                                                                    |
|4. Run the executable                                               |
|          >> ./IsingNucl                                            |
|--------------------------------------------------------------------|

|--------------------------------------------------------------------|
|               SIMULATION SETUP AND ANALYSIS                        |
|--------------------------------------------------------------------|
|Pleaser refer to "main.c" for the steps to invoke the simulation.   |
|There are three experiments pre-coded in "src/simul.c".             |
|--------------------------------------------------------------------|
|1. run_sim()                                                        |
|    Classic Metropolis Monte Carlo 2D Ising Lattice Simulation.     |
|    [simulation video](https://youtu.be/HaPEz-NQ8I4)                |
|--------------------------------------------------------------------|
|2. validate_free_energy()                                           |
|   Validating simulation by running simulations at various          |
|   temperatures and comparing the minimum energy of each            |
|   temperature with the exact solution. Run the simulation          |
|   and free energy construction in the python jupyter notebook      |
|   "free_energy.ipynb"                                              |
|--------------------------------------------------------------------|
|3. run_sim_nucleation()                                             |
|    Simulation of nucleus formation in 2D Ising System simulated    |
|    via umbrella sampling (also known as non boltzmann sampling)    |    
|    with a overlapping windows from N = 15 to 450. Refer to         |
|    nucleation.ipynb for analysis of the simulation data and free   |
|    energy construction.                                            |
|    [simulaton video](https://youtu.be/6_lvSokWUsw)                 |
|--------------------------------------------------------------------|

|--------------------------------------------------------------------|
|                     IMPORTANT POINTS                               |
|--------------------------------------------------------------------|
|1. The lattice size, defined in main (int n) is the sidelength of   |
|   the lattice square (not number of spins). That means that        |
|   number of spins = n * n.                                         |
|--------------------------------------------------------------------|
|2. The resolution of the output video should be the same as the     |
|   lattice dimensions. The resolution of the output videl is set    |
|   int "make_video()" function in "src/media.c". In the ffmpeg      |
|   commanf simply change the desired size. for example if the       |
|   lattice is n=50, change the argument to "-s 50x50".              |
|--------------------------------------------------------------------|
|3. "int* simspace" is the memory buffer holding the lattice         |
|   snapshots which will be turned into a video. It is better to     |
|   skip some steps before adding a frame to the video to avoid      |
|   stack overflow issues.                                           |
|--------------------------------------------------------------------|


|--------------------------------------------------------------------|
|                       REFERENCES                                   |
|--------------------------------------------------------------------|
| [1] Chandler, David. Introduction to modern statistical mechanics. |
|     1987.                                                          |
|                                                                    |
| [2] Brendel, Kevin, G. T. Barkema, and Henk van Beijeren.          |
|     "Nucleation times in the two-dimensional Ising model." Physical|
|     Review E 71, no. 3 (2005): 031601.                             |
|--------------------------------------------------------------------|