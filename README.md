This project is split into two parts. One part has the [simulation](https://github.com/ocatias/LCMC_Sim) (and jupyter notebooks) while the other part has the [visualisation](https://github.com/ocatias/LCMC_Vis) of the liquid crystal.

# A Monte Carlo Simulation of Biaxaial Liquid Crystals in Confinement

Liquid crystals are made out of many particles that can move and interact with each other. For the purpose of this simulation, these particles will be hard board-like particles (HPBs) this means that they basicly have the form of cuboids and cannot overlap. These particles are placed inside a system with hard walls; with the help of the Metropolis Monte Carlo method we then try to find equilibrium states and analyse the emerging patterns.
The simulation is written in C++ and uses the [Eigen Library](http://eigen.tuxfamily.org/index.php?title=Main_Page) to calculate eigenvectors and eigenvalues.   The visualisation is done in the Unity game engine (it also uses the  [Unity Standalone File Browser](https://github.com/gkngkc/UnityStandaloneFileBrowser)). Finally, some Jupyter Notebooks were used to analyse the results.

## Getting Started
Download the project and compile the simulation which is in the LCMC.cc file
```
g++ -std=c++1z LCMC.cc vector3d.h vector3d.cc box.h box.cc BoxCollision.h BoxCollision.cc RandomMove.h RandomMove.cc Statistics.h Statistics.cc -o lcmcsim.out
```

You can then start a simulation with

```
 ./lcmcsim.out output_filename particle_length particle_height density  confinement_geometry confinement_paramater1 confinement_parameter2 
```
* output_filename: The program will save the results in a file with this name 
* density: The density of the system (the total volume of the particles divided by the volume of th system), has to be a number between 0 and 1
* confinement_geometry: A number that corresponds to the confinement
* confinement_paramater: Parameters that define the size of your confinement

The possibilities for confinements and their parameters are:

0. **Cylinder** (radius, height)
1. **Cube** (length)
3. **Sphere** (radius)
4. **Half Sphere** (radius)

For example the following command starts a simulation inside a cylinder with radius 12, height 10 and a density of 0.20. The particles have a length of 2 and a height of 3
```
./lcmcsim.out test_cylinder 2 3 0.20 0 12 10
```

*Note:* you can change a lot of important parameters (like the number of timesteps) in the LCMC.cc file.

If you want to take a look at your simulation results, you can use [this](https://github.com/ocatias/LCMC_Vis) program.

### Examples
![](https://github.com/ocatias/LCMC_Sim/blob/master/Pictures/Sim_Screenshot1.png)

![](https://github.com/ocatias/LCMC_Sim/blob/master/Pictures/rod1d16_1.resized.jpg) ![](https://github.com/ocatias/LCMC_Sim/blob/master/Pictures/rod2d20_1.resized.jpg) ![](https://github.com/ocatias/LCMC_Sim/blob/master/Pictures/rod1d25_1.resized.jpg)

![](https://github.com/ocatias/LCMC_Sim/blob/master/Pictures/lyingN400_2.resized.jpg) ![](https://github.com/ocatias/LCMC_Sim/blob/master/Pictures/hs_side01.resized.jpg) ![](https://github.com/ocatias/LCMC_Sim/blob/master/Pictures/sphere002.resized.jpg)

 ![](https://github.com/ocatias/LCMC_Sim/blob/master/Pictures/OrderParamater1.resized.png) ![](https://github.com/ocatias/LCMC_Sim/blob/master/Pictures/COG.resized.png)




