subroutine SOLVERS  (Nx, Ny,Nz,dx,dy,dz, n, PDE, X, B, solv, tolerance, itol, itmax, iter)

! Biconjugate Gradient Stabilized Method (BCGSTAB) and Generalized Minimal RESidual method (GMRES)  
! to solve a linear System: L(X) = B  is presented,
! where L() - linear operator corresponding to the finite-difference approximation PDE.
! Íere use a finite-difference approximation of Laplace's equation L(X) = Laplacian(X) for a three-dimensional vector field,
! and the finite-difference approximation of the equation L(X) = grad(div(X)) - Laplacian(X) 
! (this is equivalent L(X)=Curl(Curl(X)) without using a matrix. 
! For the finite difference approximation of the Laplacian, a variant is given using a sparse matrix.
! 
! The Fortran code is based on published algorithms:
! ---the algorithm BCGSTAB is presented in Xianyi Zeng's lectures :  Algorithm 2.3 Lecture Note 7 
! https://utminers.utep.edu/xzeng/2017spring_math5330/MATH_5330_Computational_Methods_of_Linear_Algebra_files/ln07.pdf
! ---the algorithm GMRES is presented in Wikipedia (Matlab/Octave version, some comments saved): 
! https://en.wikipedia.org/wiki/Generalized_minimal_residual_method

! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )

implicit none

! character string corresponding to the name of methods: BCGSTAB or GMRES
character(*), intent(IN)   :: solv 

! PDE - character string corresponding to the name of the finite-difference approximation PDE
! Use only "Laplacian"  or "CurlCurl", or "Msparse" . Case sensitive
character(*), intent(IN)   :: PDE 

INTEGER, intent(IN) :: n          ! Dimension of the system of equations
REAL(8), intent(IN) :: B(n)       ! Right-hand side vector
integer, intent(IN) :: Nx, Ny, Nz ! number of nodes in the computational domain along the X, Y, Z axes
real(8), intent(IN) :: dx,dy,dz   ! grid spacing along the X, Y, Z axes (in meters)

! itol=1: iteration stops when the 2-norm of the residual divided by the 2-norm of the right hand side is less than TOLERANCE
! itol=2: iteration stops when the 2-norm of the residual is less than TOLERANCE
INTEGER, intent(IN) :: itol
REAL(8), intent(IN) :: TOLERANCE ! Convergence criterion.

INTEGER, intent(IN) :: ITMAX      ! maximum number of iterations

! iter - the number of iterations obtained to achieve convergence, 
! or ITMAX if the convergence criterion cannot be achieved in ITMAX iterations:
INTEGER, intent(out):: iter  

!  On input X is the initial guess for solution vector
!  On output X is the final approximate solution
REAL(8), intent(INout):: X(n) 

! internal variables
real(8), allocatable :: valA(:)
integer, allocatable :: jcol(:), irow(:)
integer nnz

INTEGER i,j,k,l, m, Ndim, kdz, nim,nip,njm,njp,nkp,nkm
real(8) sdx,sdy,sdz, fi, w, R(n), Btst(n)

! START of calculation
Ndim= Nx*Ny*Nz
kdz = Nx*Ny
sdx = 1.d0/(dx*dx); sdy = 1.d0/(dy*dy); sdz = 1.d0/(dz*dz)
w= 2.0d0*(sdx +sdy + sdz)

if ( PDE == 'Msparse' )  call gen_sparse_matrix

if (solv == 'BCGSTAB') then

BCGSTAB:BLOCK

    ! internal variables for BCGSTAB:
    REAL(8):: alpha, beta, omega,  rr0, Bnorm , & 
            R0(n), P(n), AP(n), S(n), AS(n)
            
    R = approximation(PDE, X) !this function calculates the finite-difference approximation of PDE

    Bnorm=norm2(B)
    if (Bnorm == 0.d0 ) then
        print*,'Warning: Right-hand side vector = 0'
        return
    endif
    R = B - R
    R0 = R  
    P = R
    
    iter=0
    do while (iter <= itmax)
        iter = iter + 1
        AP = approximation(PDE, P)
        rr0 = dot_product(R,R0)   !1.0
        alpha = rr0/dot_product(AP,R0) 
        S = R - alpha*AP
        AS = approximation(PDE, S)
        omega = dot_product(AS,S)/dot_product(AS,AS)
        X = X + alpha*P + omega*S
        R = S - omega*AS 
        if (itol == 1) then
            if ( norm2(R)/Bnorm < TOLERANCE) exit  
        else
            if ( norm2(R) < TOLERANCE) exit  
        endif
        beta = (alpha/omega)*dot_product(R,R0)/rr0
        P = R + beta*(P - omega*AP)
    enddo
    
END BLOCK BCGSTAB

elseif (solv == 'GMRES') then

GMRES: BLOCK
	! internal working variables for GMRES:
	real(8), allocatable :: H(:,:),   Q(:,:)  
	real(8), allocatable :: sn(:) , cs(:)  ,beta_r(:)  , e1(:) , q0(:), y(:)
	real(8) ::   b_norm, r_norm, temp  ! error, 
    
	allocate(sn(ITMAX), cs(ITMAX), beta_r(ITMAX+1),  y(ITMAX), source=0.d0)
	allocate(H(ITMAX+1,ITMAX), Q(n,ITMAX+1), q0(n), e1(ITMAX+1), source=0.d0)

    R = approximation(PDE, X) !this function calculates the finite-difference approximation of PDE
	R = b - R
    
	b_norm = norm2(b)
	! error  = norm2(r) / b_norm
	e1(1) = 1.0d0
	r_norm = norm2(r)
	Q(:,1) = r / r_norm
    
	! Note: this is not the beta_r scalar in section "The method" above but 
	! the beta_r scalar multiplied by e1 
	beta_r = r_norm * e1
    
    iter=1
    do while (iter < ITMAX)

        q0 = approximation(PDE, Q(:,iter) )  ! q0 - Krylov Vector

		! Modified Gram-Schmidt, keeping the Hessenberg matrix 
		do i=1,iter                                    
			H(i,iter) = dot_product(q0,Q(:, i))
			q0 = q0 - H(i,iter) * Q(:, i)
		end do 
		H(iter+1,iter) = norm2(q0)
		q0 = q0 / H(iter+1,iter)
		Q(:,iter+1) = q0

		! eliminate the last element in H ith row and update the rotation matrix 
		! apply for ith column 
		do i=1,iter-1 
			temp        =  cs(i) * H(i,iter) + sn(i) * H(i+1,iter)
			H(i+1,iter) = -sn(i) * H(i,iter) + cs(i) * H(i+1,iter)
			H(i,iter) = temp
		end do 

		! update the next sin cos values for rotation 
		temp = sqrt(H(iter,iter)*H(iter,iter) + H(iter+1,iter)*H(iter+1,iter))
		cs(iter) = H(iter,iter)   / temp          ! see http://www.netlib.org/eispack/comqr.f
		sn(iter) = H(iter+1,iter) / temp

		! eliminate H(i + 1, i) 
		H(iter,iter) = cs(iter) * H(iter,iter) + sn(iter) * H(iter + 1,iter)
		H(iter + 1,iter) = 0.0
	
		! update the residual vector 
		beta_r(iter + 1) = -sn(iter) * beta_r(iter)
		beta_r(iter)     =  cs(iter) * beta_r(iter)
	
		! error = abs(beta_r(iter + 1)) / b_norm
        if (itol == 1) then
            if ( abs(beta_r(iter + 1)) / b_norm < TOLERANCE) exit  
        else
            if ( abs(beta_r(iter + 1)) < TOLERANCE) exit  
        endif
        
		! if (error <= TOLERANCE) exit
        iter = iter + 1

	end do 
	! if TOLERANCE is not reached, k = m at this point (and not m+1) 
	! calculate the result 
	call invers(H(1:iter, 1:iter), iter)
	y = matmul(H(1:iter, 1:iter),beta_r(1:iter) )
	X = X + matmul( Q(:, 1:iter) , y)
        
end BLOCK GMRES

endif


Btst = approximation(PDE, X)
Btst = B - Btst

print '(a,i7, a, e10.3, $)', 'Solution complete: iterations=', iter, '    ||residual||=',norm2(btst)

contains

function approximation( PDE, V)
real(8)  approximation( n), V( n )
character(*) PDE

    if   (PDE == 'Laplacian') then
        approximation = -Laplacian( V )
        
    elseif(PDE == 'CurlCurl') then             !   Curl(Curl(A)) = grad(div(A)) - Laplacian(A)
        approximation = grad(div(V)) - Laplacian(V) 
        
    elseif(PDE == 'Msparse') then  
        approximation = sprsAx(V)
        
    endif
end function approximation


function Laplacian (V)
real(8)  V(n) 
real(8)  Laplacian (n)

do k = 2,Nz-1;  do j = 2,Ny-1;  do i = 2,Nx-1

! Laplacian = ddV/ddx + ddV/ddy + ddV/ddz
    m   = i   + Nx*(j-1)  + kdz*(k-1)
    nim = m-1;  njm = m-Nx;  nkm = m-kdz;
    nip = m+1;  njp = m+Nx;  nkp = m+kdz; 
!Laplacian_x
    Laplacian(m) =  +( V(nip) + V(nim) )*sdx   &
                    +( V(njp) + V(njm) )*sdy   &
                    +( V(nkp) + V(nkm) )*sdz   &
                    -w*V(m)  
!Laplacian_y
    Laplacian(m+Ndim) = +( V(nip+Ndim) + V(nim+Ndim) )*sdx   &
                        +( V(njp+Ndim) + V(njm+Ndim) )*sdy   &
                        +( V(nkp+Ndim) + V(nkm+Ndim) )*sdz   &
                        -w*V(m+Ndim)  
!Laplacian_z
    Laplacian(m+2*Ndim) = +( V(nip+2*Ndim) + V(nim+2*Ndim) )*sdx   &
                          +( V(njp+2*Ndim) + V(njm+2*Ndim) )*sdy   &
                          +( V(nkp+2*Ndim) + V(nkm+2*Ndim) )*sdz  &
                          -w*V(m+2*Ndim)  
enddo;  enddo;  enddo

end function Laplacian


function Div( V )
   real(8)  V(n) 
   real(8)  Div(Ndim)
! v(x,y,z) = div(V) = dVx/dx + dVy/dy + dVz/dz

do k = 2,Nz-1;  do j = 2,Ny-1;  do i = 2,Nx-1

    m   = i   + Nx*(j-1)  + kdz*(k-1)
    nim = m-1;  njm = m-Nx;  nkm = m-kdz;
    nip = m+1;  njp = m+Nx;  nkp = m+kdz; 

    Div(m) =  0.5d0*( (V(nip)        - V(nim))       / dx    &     ! dVx/dx
                     +(V(njp+Ndim)   - V(njm+Ndim))  / dy    &     ! dVy/dy
                     +(V(nkp+2*Ndim) - V(nkm+2*Ndim))/ dz )        ! dVz/dz
                     
end do;  end do;  end do

end function Div

function grad( v )
real(8)  v(Ndim) 
real(8)  grad (n)
! V= grad(v)
! Vx= da/dx;   Vy = da/dy;   Vz= da/dz

do k = 2,Nz-1;  do j = 2,Ny-1;  do i = 2,Nx-1

    m   = i   + Nx*(j-1)  + kdz*(k-1)
    nim = m-1;  njm = m-Nx;  nkm = m-kdz;
    nip = m+1;  njp = m+Nx;  nkp = m+kdz; 

    grad(m)        = 0.5d0*(v(nip) - v(nim)) / dx   ! Vx = dv/dx
    grad(m+Ndim)   = 0.5d0*(v(njp) - v(njm)) / dy   ! Vy = dv/dy 
    grad(m+2*Ndim) = 0.5d0*(v(nkp) - v(nkm)) / dz   ! Vz = dv/dz

end do;  end do;  end do

end function grad 
 
function Curl( V ) ! this function is not yet used
real(8)  V(n) 
real(8)  Curl (n)

do k = 2,Nz-1;  do j = 2,Ny-1;  do i = 2,Nx-1

    m   = i   + Nx*(j-1)  + kdz*(k-1)
    nim = m-1;  njm = m-Nx;  nkm = m-kdz;
    nip = m+1;  njp = m+Nx;  nkp = m+kdz; 

    Curl(m)         = 0.5d0*(V(2*Ndim+njp)- V(2*Ndim+njm)) /dy      &
                    - 0.5d0*( V(Ndim+nkp) - V(Ndim+nkm))/dz             ! Curl_x = dVz/dy - dVy/dz 

    Curl(m+Ndim)   = 0.5d0*(V(nkp) - V(nkm)) /dz                    &
                   - 0.5d0*( V(2*Ndim+nip) - V(2*Ndim+nim)) /dx         ! Curl_y = dVx/dz - dVz/dx 
            
    Curl(m+2*Ndim) = 0.5d0*(V(Ndim+nip) - V(Ndim+nim)) /dx          &
                   - 0.5d0*( V(njp) - V(njm)) /dy                       ! Curl_z = dVy/dx - dVx/dy 
                   
end do;  end do;  end do

end function Curl


function sprsAx  (V) 
! Matrix-vector multiplication in CSR format using the dot product form
    REAL(8)  :: V(n)
    REAL(8) :: sprsAx(n)
    INTEGER i1,i2

    do  i=1,n
        i1=irow(i); i2=irow(i+1)-1
        sprsAx(i) =  dot_product(valA(i1:i2), V(jcol(i1:i2) ) )
    enddo 

END function sprsAx

subroutine gen_sparse_matrix
integer nn
! GENERATION of compressed sparse matrix discretizing the 3d Laplace operator
! d^2 / dx^2 + d^2 / dy^2 + d^2 / dz^2  is created from a 7-point
!  stencil on an Nx by Ny by Nz grid.
!-----------------------------------------------------------------------------

! COUNTING the number of non-zero elements nnz
nnz = 0
do k = 1, nz;  do j = 1, ny;  do i = 1, nx
    if (i==1 .or.j==1 .or.k==1 .or.i==nx .or.j==ny .or.k==nz  ) then
        nnz = nnz + 1
    else
        nnz = nnz + 7
    endif
end do; end do; end do;
nnz = 3*nnz

! memory allocation for CSR storage format arrays
allocate (valA(nnz), source=0.d0)
allocate (jcol(nnz), irow(n+1), source=0)

! FILLING ARRAYS for CSR storage format 
sdx = 1.d0/(dx*dx); sdy = 1.d0/(dy*dy); sdz = 1.d0/(dz*dz)
fi  = 2.d0*(sdx + sdy + sdz )
L = 0;  nn = 0
irow(1) = 1
do m=1,3
    do k = 1, nz;  do j = 1, ny;  do i = 1, nx
        nn = nn + 1
        irow(nn)=L+1
        if (i==1 .or.j==1 .or.k==1 .or.i==nx .or.j==ny .or.k==nz  ) then
            L = L + 1
            jcol(L) = nn
            valA(L) = 1.d0
        else
            nim = nn-1;  njm = nn-Nx;  nkm = nn-kdz;
            nip = nn+1;  njp = nn+Nx;  nkp = nn+kdz; 
            jcol(L+1:L+7) = [ nkm,  njm,  nim, nn,  nip,  njp,  nkp]
            valA(L+1:L+7) = [-sdz, -sdy, -sdx, fi, -sdx, -sdy, -sdz]
            L = L + 7
        endif
    end do; end do; end do;
enddo
irow(nn+1)=L+1

print '(a,i3,a,i3,a,i3,     a,i9, a,g12.5,a/)', 'sparse matrix generation completed on grid (', Nx,' x ',Ny,' x ',Nz, &
                 ' ), Non zero elem= ',nnz, ' Density of matrix:', 100.0* real(nnz)/real(n)/real(n),'%'


end subroutine gen_sparse_matrix



end subroutine SOLVERS

