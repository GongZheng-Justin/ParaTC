# ParaTC
&emsp;**ParaTC** is a high order finite difference solver for simulations of turbidity currents with high parallel efficiency:
* Periodic conditions are imposed in streamwise (x) and spanwise (z) directions.![](doc/SchematicDiagram.png)
* Navier-Stokes equations,coupled with an active passive scalar transport equation, are simulated.![](doc/GoverningEquation.png)
* MPI parallelization by means of pencil distributed decomposition. In order to improve the parallel efficiency, we propose a new 2D pencil-like parallel configuration with totally 6 different pencil arrangements.![](doc/SixPencils.png)
* A [parallel Thomas algorithm](https://github.com/MPMC-Lab/PaScaL_TDMA) is included to further reduce the communication overhead when solving tridiagonal equations.![](doc/ParallelThomas.png)
* Fourth-order spatial scheme is used for periodic directions, i.e., streamwise and spanwise directions. Second-order scheme is used in wall-normal direction.
* 
