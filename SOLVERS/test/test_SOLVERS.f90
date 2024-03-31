program test_SOLVERS

! Biconjugate Gradient Stabilized Method (BCGSTAB) and Generalized Minimal RESidual method (GMRES)  
! to solve a linear System: L(X) = B  is presented,
! where L() - linear operator corresponding to the finite-difference approximation PDE.
! Íere use a finite-difference approximation of Laplace's equation L(X) = Laplacian(X) for a three-dimensional vector field,
! and the finite-difference approximation of the equation L(X) = grad(div(X)) - Laplacian(X) 
! (this is equivalent L(X)=Curl(Curl(X)) without using a matrix. 
! For the finite difference approximation of the Laplacian, a variant is given using a sparse matrix.

! The Fortran code is based on published algorithms:
! ---the algorithm BCGSTAB is presented in Xianyi Zeng's lectures :  Algorithm 2.3 Lecture Note 7 
! https://utminers.utep.edu/xzeng/2017spring_math5330/MATH_5330_Computational_Methods_of_Linear_Algebra_files/ln07.pdf
! ---the algorithm GMRES is presented in Wikipedia (Matlab/Octave version): 
! https://en.wikipedia.org/wiki/Generalized_minimal_residual_method

implicit none 

!  INPUTS:
!-----------------------------------------------------------------------------
integer Nx, Ny, Nz              ! number of nodes in the computational domain along the X, Y, Z axes
real(8) dx, dy, dz              ! grid spacing along the X, Y, Z axes (in meters)

! PDE - character string corresponding to the name of the finite-difference approximation PDE:
character(:), allocatable :: PDE  ! 

! solv - character string corresponding to the name of methods:
character(:), allocatable :: solv 

real(8),  allocatable :: B(:)    ! Right-hand side vector (the field source)
real(8),  allocatable :: Btst(:) ! Copy of right-hand side vector for testing
real(8) Amx, Amy                 ! Quantities for calculating the field source

! itol=1: iterations stop when the 2-norm (or abs for GMRES) of the residual divided by the 2-norm of the right hand side is less than TOLERANCE
! itol=2: iterations stop when the 2-norm (or abs for GMRES) of the residual is less than TOLERANCE
INTEGER :: itol       ! iterations stop condition
REAL(8) :: TOLERANCE  ! convergence criterion
INTEGER :: ITMAX      ! maximum number of iterations

!  On input  X is the initial guess for solution vector
!  On output X is the final approximate solution:
real(8),  allocatable :: X(:)

!  OUTPUTS:
!-----------------------------------------------------------------------------
! iter - the number of iterations obtained to achieve convergence, 
! or ITMAX if the convergence criterion cannot be achieved in ITMAX iterations:
INTEGER iter 

INTEGER nnz                       ! number of non-zero elements of a sparse matrix

! names of files with csv and vtk extensions in which the solution results will be written:
character(:), allocatable :: files

! internal working variables:
integer :: i,j,k,L,m,n,kdz,Ndim
integer :: kim,kjm,kkm,kip,kjp,kkp
real(8) :: sdx, sdy, sdz, fi

integer k0,k1,k2                   ! variables to measure calculation time
real Tcalc

! ENTERING INPUT DATA
!-----------------------------------------------------------------------------
Nx = 50; Ny = 50; Nz = 40          ! number of nodes along the X, Y, Z axes
dx  = 1.d0; dy = 1.d0; dz = 1.d0   ! ggrid spacing along the X, Y, Z axes (in meters)

! uncomment/comment what you need
solv='BCGSTAB'
! solv='GMRES'

! PDE approximation type
! uncomment/comment what you need
PDE='Laplacian'                    ! Use only "Laplacian"  or "CurlCurl"  or "Msparse". Case sensitive
! PDE='CurlCurl'                        ! Note:  Curl(Curl(A)) = grad(div(A)) - Laplacian(A)
! PDE='Msparse'

itmax =  600                       ! Maximum number of iterations
tolerance = 1.0d-5                 ! Convergence criterion
itol=1                             ! iterations stop condition
files='X3d'                        ! names of files with results 

! START of calculation
!-----------------------------------------------------------------------------
kdz = Nx*Ny
Ndim = Nx*Ny*Nz  
allocate( B(3*Ndim), Btst(3*Ndim), X(3*Ndim), source=0.d0 )

! SETTING Dirichlet boundary conditions
!-----------------------------------------------------------------------------
if (PDE == 'Laplacian' .or. PDE=='CurlCurl' ) then
    call boundary(x)
elseif (PDE == 'Msparse' ) then
    call boundary(B)
ELSE
   print*, PDE, ': use only Laplacian or CurlCurl or Msparse. Case sensitive'
   stop 'error in name PDE'
endif

! CALCULATION of the vector of the right side of the system of linear equations : in this version, this is 
! the distribution of the field source vector along the ellipse with semi-axes Amx Amy, 
! located in the XY plane in the center of the calculation area
!-----------------------------------------------------------------------------
Amx = real(Nx/4); Amy = real(Ny/4); 
do n = 1,360
    fi = real(n)
    ! ellipse coordinates 
    i = Nx/2 + int(Amx * cosd(fi))
    j = Ny/2 + int(Amy * sind(fi))
    k =  i + Nx*(j-1) + kdz*(Nz/2-1)
    ! X and Y components of the source vector directed tangentially to the ellipse.
    ! Amplitude is 10
    B(k) = -10.d0*sind(fi);  B(k+Ndim) = +10.d0*cosd(fi); 
enddo

! NUMERICAL CALCULATION
!-----------------------------------------------------------------------------
print '(a,i3,a,i3,a,i3,a)', 'START test calculation 3d field, with options "'//PDE//'" and "'//solv//'" on grid {', Nx,' x ',Ny,' x ',Nz, ' )'

call system_clock(k1,k0);
    call SOLVERS  (Nx,Ny,Nz,dx,dy,dz, 3*Ndim, PDE, X, B, solv, tolerance, itol, itmax, iter)
call system_clock(k2)
print '(a, e10.3 )',  '     Calculation time=', real(k2-k1)/real(k0) 

! Saving results to files
!-----------------------------------------------------------------------------
print*
print '(a)', 'RESULT of solution with options "'//PDE//'" and "'//solv//'" saved in files:'
call writeCsvVtk (Nx,Ny,Nz,dx,dy,dz, X, B, files)

contains

subroutine boundary(V) !  Dirichlet boundary conditions
REAL(8)  :: V(3*Ndim)
    n=0
    do k = 1,Nz;  do j = 1,Ny;   do i = 1,Nx
        n=n+1
        if (i==1 .or. i==Nx  ) then
            V(n)=0.d0;  V(n+Ndim)=0.d0;  V(n+2*Ndim)=0.d0;  !boundary for X
        elseif (j==1 .or. j==Ny ) then
            V(n)=0.d0;  V(n+Ndim)=0.d0;  V(n+2*Ndim)=0.d0;  !boundary for Y
        elseif (k==1 .or. k==Nz ) then
            V(n)=-0.d0; V(n+Ndim)=0.d0;  V(n+2*Ndim)=0.d0;   !boundary for Z
        endif
     enddo; enddo;enddo
end subroutine boundary


end program test_SOLVERS



