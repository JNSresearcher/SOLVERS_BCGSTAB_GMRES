! Biconjugate Gradient Stabilized Method (BCGSTAB) and Generalized Minimal RESidual method (GMRES)  
! to solve a linear spaese System: Ax = B  

! The Fortran code is based on published algorithms:
! ---the algorithm BCGSTAB is presented in Xianyi Zeng's lectures :  Algorithm 2.3 Lecture Note 7 
! https://utminers.utep.edu/xzeng/2017spring_math5330/MATH_5330_Computational_Methods_of_Linear_Algebra_files/ln07.pdf
! ---the algorithm GMRES is presented in Wikipedia (Matlab/Octave version): 
! https://en.wikipedia.org/wiki/Generalized_minimal_residual_method

program example_SOLVERS_base
implicit none

!  INPUTS
!-----------------------------------------------------------------------------
integer Ndim                        ! Dimension of the system of equations
real(8), allocatable :: A(:,:)      ! Nonsingular matrix of general form
real(8), allocatable :: B(:)        ! Right-hand side vector (the field source)
real(8), allocatable :: Btst(:)     ! Copy of right-hand side vector for testing

REAL(8) :: TOLERANCE  ! convergence criterion.
INTEGER :: ITMAX      ! maximum number of iterations

! itol=1: iterations stop when the 2-norm (or abs for GMRES) of the residual divided by the 2-norm of the right hand side is less than TOLERANCE
! itol=2: iterations stop when the 2-norm (or abs for GMRES) of the residual is less than TOLERANCE
INTEGER :: itol       ! iterations stop condition

character(:), allocatable :: solv ! character string corresponding to the name of methods: BCGSTAB or GMRES

!  On input X is the initial guess for solution vector
!  On output X is the final approximate solution:
real(8),  allocatable :: X(:)

!  OUTPUTS
!-----------------------------------------------------------------------------
! Matrix A in the format CSR (Compressed Sparse Row) (see https://en.wikipedia.org/wiki/Sparse_matrix)
! valA() contains non-zero values of matrix A
! jcol() contains the indexes of columns for each row
! irow(j) points to the beginning of column j in the lists valA() and jcol()
real(8), allocatable :: valA(:)
integer, allocatable :: jcol(:), irow(:)

INTEGER nnz    ! number of non-zero elements of a sparse matrix

! iter - the number of iterations obtained to achieve convergence, 
! or ITMAX if the convergence criterion cannot be achieved in ITMAX iterations:
INTEGER iter 

integer :: k0,k1,k2     ! variables to measure calculation time
real Tcalc

integer ::i   ! internal working variable

! ENTERING INPUT DATA
!-----------------------------------------------------------------------------
Ndim=4                             ! Dimension of the system of equations
itmax =10                          ! Maximum number of iterations
tolerance = 1.0d-6                 ! Convergence criterion
itol=2                             ! iterations stop condition

! uncomment/comment what you need:
! solv='BCGSTAB'
solv='GMRES'

! Allocation of memory and input of matrix A
allocate(A(Ndim,Ndim), B(Ndim), Btst(Ndim), X(Ndim), source=0.d0 )
A(1,:) =  [ 4.,  0., -1.,  0.]
A(2,:) =  [-2.,  4.,  0., -1.]
A(3,:) =  [ 0., -3.,  4.,  0.]
A(4,:) =  [ 0.,  0.,  0. , 4.]

! call random_number(A)    ! full matrix for other sizes

Print*,'Matrix A in full format:'
do i=1,Ndim
    print '( *(1x,f10.3) )',  A(i,:)
enddo

! Entering the vector of the right side of the system of equations
B = 1.00d0

! counting the number of non-zero elements
nnz = count(A /= 0.d0)

Print*,' '
print '(a,i9,a,i9, a,g12.5,a /)', 'Ndim=', Ndim, ' non zero elem= ',nnz, ' Density of matrix:', 100.0* real(nnz)/real(Ndim*Ndim),'%'

! memory allocation for CSR storage format arrays
allocate(valA(nnz), source=0.d0)
allocate(jcol(nnz), irow(Ndim+1) , source=0)

! converting a matrix  A from full format to CSR format
call m2crs(A, Ndim, valA,  irow, jcol )
deallocate(A)

print '(a)', 'vector of the right side of the system of equations:'
print '( a, *(a,e10.3) )', 'B =', (' ',B(i), i=1,Ndim )
Print*
Print*,'Matrix A in CSR format:'
write(*, '(  a, *(i11) )')     ' index i =', (i, i=1,Ndim+1)
write(*, '(  a, *(i11) )')     ' irow(i) =', (irow(i),i=1,Ndim+1) 
write(*, '(  a )')             ' '
write(*, '(  a, *(i11) )')     ' index k =', (i, i=1,nnz)   
write(*, '(  a, *(i11) )')     ' jcol(k) =', (jcol(i),i=1,nnz) 
write(*, '(  a, *(a,f10.3) )') ' valA(k) =', (' ',valA(i),i=1,nnz)  


! NUMERICAL CALCULATION of sparse linear system 
!-----------------------------------------------------------------------------
call system_clock(k1,k0); 
    call SOLVERS_base (valA, irow, jcol, Ndim, X, B, solv, tolerance, itol, itmax, iter)
call system_clock(k2)

Tcalc = real(k2-k1)/real(k0)

Print*
print '(a)', 'RESULT of solve '//solv
print '( a, *(a,e10.3) )', 'X =', (' ',X(i), i=1,Ndim )
Print*
print '(a,i7, a, e10.3,$)', 'iterations=', iter, '   calculation time=',Tcalc

! residual check and printing of  the calculation time
!-----------------------------------------------------------------------------
call sprsAx(valA, irow, jcol, x, Btst, Ndim)
Btst = B - Btst

print '( a,e12.5)' ,   '      ||residual||=',norm2(Btst)

end program example_SOLVERS_base

