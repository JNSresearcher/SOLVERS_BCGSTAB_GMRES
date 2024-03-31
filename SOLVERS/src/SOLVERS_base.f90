SUBROUTINE SOLVERS_base (valA, irow, jcol, n, X, B, solv, tolerance, itol, itmax, iter)

! Biconjugate Gradient Stabilized Method (BCGSTAB) and Generalized Minimal RESidual method (GMRES)  
! to solve a linear spaese System: Ax = B  is presented

! The Fortran code is based on published algorithms:
! ---the algorithm BCGSTAB is presented in Xianyi Zeng's lectures :  Algorithm 2.3 Lecture Note 7 
! https://utminers.utep.edu/xzeng/2017spring_math5330/MATH_5330_Computational_Methods_of_Linear_Algebra_files/ln07.pdf
! ---the algorithm GMRES is presented in Wikipedia (Matlab/Octave version): 
! https://en.wikipedia.org/wiki/Generalized_minimal_residual_method

! Fortran code created by J.Sochor   ( https://github.com/JNSresearcher )

implicit none 

! character string corresponding to the name of methods: BCGSTAB or GMRES
character(*), intent(IN)   :: solv 

! Matrix A in Compressed Sparse Row (CSR) format:
INTEGER, intent(IN) :: jcol(*), irow(n+1)
REAL(8), intent(IN) :: valA(*)

INTEGER, intent(IN) :: n          ! Dimension of the system of equations
REAL(8), intent(IN) :: B(n)       ! Right-hand side vector

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
REAL(8):: R(n)
integer :: i 

if (solv == 'BCGSTAB') then

BCGSTAB:BLOCK
    ! internal variables for BCGSTAB:
    REAL(8):: alpha, beta, omega,  rr0, Bnorm , & 
            R0(n), P(n), AP(n), S(n), AS(n)

    iter=0 
    call sprsAx(valA, irow, jcol, X, R, n)
    
    Bnorm = norm2(b)    
    if (Bnorm == 0.d0 ) then
        print*,'Warning: Right-hand side vector = 0'
        return
    endif
     
    R = B - R; R0 = R
    P = R
    do while (iter <= ITMAX)
        iter = iter + 1
        call sprsAx(valA, irow, jcol, P, AP, n)
        rr0 = dot_product(R,R0)
        alpha = rr0/dot_product(AP,R0) 
        S = R - alpha*AP
        call sprsAx(valA, irow, jcol, S, AS, n)
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
	real(8) ::   b_norm, r_norm, temp 
	
	allocate(sn(ITMAX), cs(ITMAX), beta_r(ITMAX+1),  y(ITMAX), source=0.d0)
	allocate(H(ITMAX+1,ITMAX), Q(n,ITMAX+1), q0(n), e1(ITMAX+1), source=0.d0)

	! use x as the initial vector 
	call sprsAx ( valA, irow, jcol, x, R, n )
	R = b - R

	b_norm = norm2(b)
	e1(1) = 1.0d0
	r_norm = norm2(r)
	Q(:,1) = r / r_norm

	! Note: this is not the beta_r scalar in section "The method" above but 
	! the beta_r scalar multiplied by e1 
	beta_r = r_norm * e1

	do iter=1,ITMAX

		call sprsAx(valA, irow, jcol, Q(:,iter), q0, n)  ! q0 - Krylov Vector

		! Modified Gram-Schmidt, keeping the Hessenberg matrix 
		do i=1,iter                                    
			H(i,iter) = dot_product(q0,Q(:, i))
			q0 = q0 - H(i,iter) * Q(:, i)
		end do 
		H(iter + 1,iter) = norm2(q0)
		q0 = q0 / H(iter + 1,iter)
		Q(:,iter+1) = q0

		! eliminate the last element in H ith row and update the rotation matrix 
		! apply for ith column 
		do i=1,iter-1 
			temp     =  cs(i) * H(i,iter) + sn(i) * H(i+1,iter)
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
	
        if (itol == 1) then
            if ( abs(beta_r(iter + 1)) / b_norm < TOLERANCE) exit  
        else
            if ( abs(beta_r(iter + 1)) < TOLERANCE) exit  
        endif
	
	end do 
	! if TOLERANCE is not reached, k = m at this point (and not m+1) 
	! calculate the result 
	call invers(H(1:iter, 1:iter),iter)
	y=matmul(H(1:iter, 1:iter),beta_r(1:iter))
	X = X + matmul( Q(:, 1:iter) , y)
    
end BLOCK GMRES
	
else
   print*, solv, ': use only BCGSTAB or GMRES. Case sensitive'
   stop 'error in name solv'
endif

end SUBROUTINE SOLVERS_base


