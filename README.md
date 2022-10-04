# lamb_mode_elasticity
Code with the article Lamb modes and Born approximation for small shape defects inversion in elastic plates

Main programs: 
  - visu_lamb_ZGV.m
Programm to visualize the propagation of a wavefield in an elastic waveguide in 2D. This program generates the simulation using the Lamb decomposition presented in the article and compare it the a FEM/PML generated wavefield. It is used to produce Figure 6. It is possible to run a fast version of this program by commenting L63-L76

- recloinres.m
Program to recover width defects in waveguides in 2D. This program generates the data using FEM/PML methods and use the reconstruction method presented in the article to provide an approximation of the width defect. It is used to produce Figure 9 and Table 1. It is possible to run a fast version of this program by commenting L52-L56

- LAMB3D.m
Program to visualize the propagation of a wavefield in an elastic waveguide in 3D. This program generates the simulation using the Lamb decomposition presented in the article. It is possible to visualize the wavefield on a surface/on a side of the plate by commenting the approxiate section of the code in the sub-program solveLamb3D.m. This program is used to produce the left side of Figures 7 and 8. 

- LAMB_3D.edp, LAMB_3D2.edp and LAMB_3D.m
Program to visualize the propagation of a wavefield in an elastic waveguide in 3D. The .edp files are runned using Freefem++ and generate the wavefield using FEM/PML methods. Then, the code LAMB_3D.m is used to visualize the data from the .edp files and produce right side of Figures 7 and 8. 

