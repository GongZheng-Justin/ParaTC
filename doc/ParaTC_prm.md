# input parameter for `ParaTC` 

&emsp;Input file for `ParaTC` can be divided into 8 part. Every part contains one **namelist**, i.e. `BasicParam`, `MeshSection`, `MeshOptions`, `SpectraOptions`, `Spectra2DOpt`, `IO_Options`, `SaveScalarOption`, and `ScalarFlowOptions`.

&emsp;**NOTE:**
* The parameter sequence in a **namelist** is not important.
* A **namelist** is started with `&NamelistName` (e.g. `&BasicParam`), and ended with `/End`.
* You can also add blank lines or comments (after sign `!`) in the input file, just like what you do in a free format Fortran code.
