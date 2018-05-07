NPRE 555 CP2 “P1 and P3 Numerical Solutions for Neutron Transport in a 1D infinite slab”
By Muhammad Abdelghany
PhD Student at the NPRE department at the UIUC
***********************************************

This code is made using MATLAB and is Composed of 2 folders, one is for P1 solution and
the other is for P2 solution. Each folder have 2 files, one for the solution using Mark’s
Boundary conditions at the vacuum boundary and the other is for the solution using Marshak’s
Boundary Conditions.

The Details of Each File:
*************************

1. “P1Mark.m” in the P1 folder
------------------------------

In this file, you can calculate the P1 approximation of the scalar flux in the slab assuming
Mark’s boundary conditions at the vacuum boundaries, where the scalar flux is assumed to 
vanish there. The code starts with specifying the width of the slab (L) and the cross-sections,
then the number of mesh points (k) which you can changed according to your preference of 
accuracy. The code then builds the (S) and (A) matrices and in the first step it takes care
of the Mark’s BCs by specifying the 1st and kth rows of the (A) Matrix.
Next it finds the inverse of the (A) matrix and calculate the flux matrix (phi). 
At the end, a plot of flux (phi) vs the width of the slab (x) is generated.

2. “P1Marshak.m” in the P1 folder
---------------------------------

This file is exactly similar to the previous one “P1Mark.m”, and starts with specifying the
 number of space points (k) and the width of the slab (L). The only difference is that, 
while building the (S) and (A) matrices, in the first step it takes care of the Marshak’s 
BCs by specifying the 1st and kth rows of the (A) Matrix. Marshak’s boundary conditions 
assumes that the  incoming neutron currents are zero at the vacuum boundaries. 

3. “P3Mark.m” in the P3 folder
------------------------------

In this file, you can calculate the P3 approximation of the scalar flux in the slab assuming
Mark’s boundary conditions at the vacuum boundaries, where the scalar flux is assumed to 
vanish there. The code starts with specifying the width of the slab (L) and the 
cross-sections, then the number of mesh points (k) which you can changed according to your
preference of accuracy. The code is developed to implement an iterative scheme of the (F1)
and (F0) and to find the scalar flux after matching a certain conversion criterion (mx_dif).
The code starts with specifying an arbitrary initial guess for (F0) and (F1) corresponding
to the 1st iteration, then it calculates the values of (F0) and (F1) at each mesh point for
each iteration. Then, at each iteration it assumes the initial values of (F1) of the two 
boundary points that are not possible to be calculated from the (F1) equations to be equal
to the corresponding values of the nearest neighboring points.
Next, for each iteration the code implement Mark’s BCs for the two boundary points for (F0)
interms of (F1). Then it calculates the scalar flux (phi) for each mesh point and calculates
 the difference in (phi) between each two successive iteration for each mesh point and 
 ompare the maximum of this difference with the pre-specified convergence criterion (mx_dif).
After the convergence criterion is satisfied, the code displays the number of iterations (n)
required and then generate a plot of the flux (phi) vs the width of the slab (x).

2. “P3Marshak.m” in the P3 folder
---------------------------------

This file is exactly similar to the previous one “P3Mark.m”, and starts with specifying the
number of space points (k) and the width of the slab (L). The only difference is that, it 
implement the equations required by applying the Marshak’s boundary conditions during each
iteration for the two boundary points for (F0) in terms of (F1). At the end, the code 
displays the number of iterations (n) required to satisfy the conversion criterion and then
generate a plot of the flux (phi) vs the width of the slab (x).

