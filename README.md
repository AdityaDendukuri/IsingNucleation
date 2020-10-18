


## INTRODUCTION                                
This is a Monte Carlo simulation program for the Ising spin model. Simulations include free energy construction [1], nucleus formation [2] with umbrella sampling [1]. This program is written in C with ffmpeg used for generation of the simulation video. Please make sure ffmpeg is installed in your system (can be found [here](https://ffmpeg.org/download.html)). Python is used to visualize the thermodynamics of the system via jupyter notebooks for the two specific experiments for validation, which can be found in analysis/ .  

## INSTRUCTIONS TO RUN 
1. Clone the repo
```
$ git clone https://github.com/AdityaDendukuri/IsingNucleation.git
```
2. Create the directories build/ and output/                        
``` 
$ cd IsingNucleation  
$ mkdir build output                                     
```                                    
                                                                    
3. Go to build directory and run cmake as follows (you only need to do this the first time)                   
```
$ cd build  
$ cmake -S ../ -B .                                      
```
          
4. This will generate the makefile in build/. Compile the program by simply   
   typing:       
```
$ make                                                   
```                                                       
5. Run the executable           
 ```
 $ ./IsingNucl                                            
 ```
 Simulation output and video will be saved in output/ .
## SIMULATION SETUP AND ANALYSIS                        
Pleaser refer to "main.c" for the steps to invoke the simulation. There are three experiments pre-coded in "src/simul.c".             

1. run_sim()                                                        
    Classic Metropolis Monte Carlo 2D Ising Lattice Simulation. [simulation video](https://youtu.be/HaPEz-NQ8I4)                
1. validate_free_energy()                                           
   Validating simulation by running simulations at various temperatures and comparing the minimum energy of each temperature with the exact solution. Free energy validation based on the output of the experiment is shown in the python jupyter notebook in analysis/ folder "free_energy.ipynb".                                              
1. run_sim_nucleation()                                             
   Simulation of nucleus formation in 2D Ising System simulated via umbrella sampling (also known as non boltzmann sampling) with a overlapping windows from N = 15 to 450. Refer    to nucleation.ipynb in analysis/ folder for analysis of the simulation data and free energy construction. [simulaton video](https://youtu.be/6_lvSokWUsw).     

## NOTES                             
1. The lattice size, defined in main (int n) is the sidelength of the lattice square (not number of spins). That means that number of spins = n * n.        

2. The resolution of the output video should be the same as the lattice dimensions. The resolution of the output video is set in the function "make_video()" in "src/media.c". The resolution should be defined explicitly in the string with the ffmpeg command. For example if the lattice is n=50, change the resolution as "ffmpeg <args> -s 50x50 <rest_of_the_command>".         

3. "int* simspace" is the memory buffer holding the lattice snapshots which will be turned into a video. It is better to skip some steps before adding a frame to the video to      avoid stack overflow issues.                                           

## REFERENCES                                   
 [1] Chandler, David. Introduction to modern statistical mechanics. 1987.                                                          
                                                                    
 [2] Brendel, Kevin, G. T. Barkema, and Henk van Beijeren. "Nucleation times in the two-dimensional Ising model." Physical Review E 71, no. 3 (2005): 031601.                             
