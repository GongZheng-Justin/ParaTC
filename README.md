# ParaTC
&emsp;**ParaTC** is a high order finite difference solver for simulations of turbidity currents with high parallel efficiency:
* MPI parallelization by means of pencil distributed decomposition. In order to improve the parallel efficiency, we propose a new 2D pencil-like parallel configuration with totally 6 different pencil arrangements.![](doc/SixPencils.png)
* A [parallel Thomas algorithm](https://github.com/MPMC-Lab/PaScaL_TDMA) is included to further reduce the communication overhead when solving tridiagonal equations.![](doc/ParallelThomas.png)
* An optimal search method is performed in the initializing stage to find the fast Poisson solver scheme among four alternatives for specific mesh configuration. The runtime ratio between traditional pencil-like Poisson solver and present solver is about 1.5.
```
  An example of the optimal PPE method search. Corresponding grid number is 9216×140×1400.
  Auto-tuning mode for Poisson Solver......
    Choice-1, time= 4350.00035871624
    Choice-2, time= 4099.76840879989
    Choice-3, time= 3699.08386557852
    Choice-4, time= 4507.20634813479
  The best Poisson Solver choice is probably Choice-3
```
*  Periodic conditions are imposed in streamwise (x) and spanwise (z) directions.![](doc/SchematicDiagram.png)
* Navier-Stokes equations,coupled with an active passive scalar transport equation, are simulated.![](doc/GoverningEquation.png)
* Fourth-order spatial scheme is used for periodic directions, i.e., streamwise and spanwise directions. Second-order scheme is used in wall-normal direction.
*  An approximate linear strong scaling performance is achieved, and the weak scaling performance is also improved.![](doc/Scaling.png)
