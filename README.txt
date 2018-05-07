NPRE 555 CP3 “Discrete Ordinate Solution Using Radau's Quadrature Set”
By Muhammad Abdelghany
PhD Student at the NPRE department at the UIUC
***********************************************

This code is made using MATLAB and is Composed of 3 folders. The first one is for Solving the discrete
ordinate using Radau's Quadrature Set and then comparing it with Gauss-Legendre's and Chebyshev’s 
quadrature set. This is used to solve a one-dimension slab geometry, with vacuum boundary conditions,
making use of Marshak’s boundary conditions. The second one, is to study the spatial convergence for a
given order of Radau's Quadrature. Finally, the third one, is to study the accuracy convergance of this 
approach.

The Details of Each File:
*************************

1. “cp3_Radau_Q1” in the Q1 folder
----------------------------------

In this file, you can calculate the partial fluxes, rightward angular flux, leftward angular flux,
and scalar flux with implementing the Radau's Quadrature Set, using 4 quadrature and 120 space meshes.

2. “cp3_Gauss_Q1” in the Q1 folder
----------------------------------

This file is exactly doing the same as the previous file but with implementing the Gauss-Legendre's 
Quadrature Set.

3. “cp3_Compare_Gauss_Q1”, “cp3_Compare_Radau_Q1”, and “cp3_Compare_Chebyshev_Q1” in the Q1 folder
--------------------------------------------------------------------------------------------------

These are three files used simultaneously to generate the graphs for comparing the results of implementing
Radau's Quadrature, and the  results obatined by applying the Gauss-Legendre's and Chebyshev’s approaches.

3. “cp3_Radau_Q2_2.m” in the Q2.2 folder
----------------------------------------

In this file, we do a spatial convergence study by changing the number of space meshes in a systematic 
way and calculate the difference in the scalar flux between the different iterations till it match a certain
convergence criterion. This study is done using 4 Radau's quadrature. A graph for the error estimated vs the
number of spatial meshes is generated.

2. “cp3_Radau_Q2_3.m” in the Q2.3 folder
----------------------------------------

In this file, we carry out an accuracy convergance study. By fixing the number of space meshes to be 
sufficiently a big number like (81), then we study the accuracy convergence of using the Radau’s quadrature
with increasing the considered order and calculate the flux change between these iterations. The values
obtained are plotted verses the order of the used quadrature.
