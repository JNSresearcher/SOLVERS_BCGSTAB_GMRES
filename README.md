The BiConjugate Gradient STABilization (BiCGSTAB) and Generalized Minimum RESidual (GMRES) methods for calculating linear sparse systems  
of equations are presented, and the application of these methods to the calculation of three-dimensional fields is shown.

 &emsp;Fortran codes are based on published algorithms:  
 
- the algorithm BCGSTAB is presented in Xianyi Zeng's lectures :  Algorithm 2.3 Lecture Note 7 
<https://utminers.utep.edu/xzeng/2017spring_math5330/MATH_5330_Computational_Methods_of_Linear_Algebra_files/ln07.pdf>   
- the algorithm GMRES is presented in Wikipedia (Matlab/Octave version): 
<https://en.wikipedia.org/wiki/Generalized_minimal_residual_method>   


### Repository Structure
&emsp; The _/SOLVERS_ directory contains :   
 
 * _/src_  
    - file _SOLVERS_base.f90_  contains an implementation of the BCGSTAB or GMRES methods for solving a linear sparse system: Ax = B. 
    - file _SOLVERS.f90_ oriented on implementation of BCGSTAB or GMRES for calculating 3D fields. Explanations are given below.
    - file _utilites.f90_ contains subroutines for calculating the product of a sparse matrix by a vector,  
converting a matrix from full format to sparse CSR format, inverting a matrix and a subroutine for outputting  
calculation results to csv and vtk files. Paraview may be required to visualization  vtk (free download: <https://www.paraview.org/>).
 * _/example_ - contains _example_SOLVERS_base.f90_  A simple example of solve a sparse system Ax=B using GMRES or BCGSTAB methods.  
A small matrix is entered in full format, which is automatically converted to CSR format.  
 * _/test_ - contains _test_SOLVERS.f90_ this is a test for solving systems of equations for calculating  
three-dimensional vector fields using the BCGSTAB or GMRES methods. Here the field equations are presented  
in the format: $L(X) = B,$ where $L()$ is a linear operator corresponding to the finite-difference approximation of the PDE.  
This version uses a finite-difference approximation of the Laplace equation $L(X) = Laplacian(X)$ and  
a finite-difference approximation of the equation $L(X) = grad(div(X)) - Laplacian(X)$ (this is equivalent to $L(X)=Curl(Curl(X))$  
without using a matrix.  Alternatively, for the Laplacian, on can use a sparse matrix,  
which is automatically generated to solve the equation $Ax=B$.  
&emsp; All files contain comments.
