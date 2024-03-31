SUBROUTINE m2crs(A, n, valA,  irow, jcol)
! Convert the full matrix to CSR format

implicit none
integer, intent(IN)  :: n 
real(8), intent(IN)  :: A(n,n)
real(8), intent(OUT) :: valA(*)
integer, intent(OUT) :: irow(n+1), jcol(*)

integer i,j,k

    k=0
    irow(1)=1
    do i=1,n 
        do j=1,n 
            if( abs(a(i,j)) > 0.d0) then
                    k = k + 1
                    valA(k) = a(i,j)
                    jcol(k) = j
            endif
        enddo 
        irow(i+1) = k + 1
    enddo 
END SUBROUTINE m2crs

SUBROUTINE sprsAx(valA, irow, jcol,x,b,n) 
! Matrix-vector multiplication in CSR format using the dot product form

    implicit none
    INTEGER, intent(IN)  :: n, jcol(*),irow(n+1)
    REAL(8), intent(IN)  :: valA(*),X(n)
    REAL(8), intent(out) :: B(n)
    INTEGER i,i1,i2

    do  i=1,n
        i1=irow(i); i2=irow(i+1)-1
        b(i) =  dot_product(valA(i1:i2), x(jcol(i1:i2) ) )
    enddo 

END SUBROUTINE sprsAx


SUBROUTINE invers(a,n)

IMPLICIT NONE
INTEGER         ,INTENT(in)    :: n
DOUBLE PRECISION,INTENT(inout) :: a(n,n)

INTEGER :: k,i,j
real(8) :: p,sp
DO k=1,n
    p=1.d0/a(k,k)
    DO j=1,N
        a(k,j)=a(k,j)*P
    ENDDO
    a(k,k)=p
    DO i=1,n
        IF(i == k) CYCLE
        sp=-a(i,k)
        DO j=1,n
            a(i,j)=a(i,j)+sp*a(k,j)
        ENDDO
        a(i,k)=sp*p
    ENDDO
ENDDO

END SUBROUTINE invers


subroutine writeCsvVtk (Nx, Ny, Nz, dx,dy,dz, X, B, files)
implicit none

integer,      intent(IN):: Nx,Ny,Nz
real(8),      intent(IN):: X(*) , B(*),dx,dy,dz
character(*), intent(IN) ::files
    
character (len = 1)  ci
character (len = 30) buf1
integer:: i,j,k,n,kdz,Ndim, nim,njm,nkm,nip,njp,nkp,ios,n_new
real(8):: s,sx,sy,sz
    
    kdz=Nx*Ny
    Ndim = Nx * Ny * Nz

! 3D field slice at z/2 level in XY plane - only for csv
    open(newunit=n_new,file=trim(files)//'.csv' )

    write(n_new, '(a)'  ) 'x,      y,     Xm,     Xx,    Xy,    Xz'
    k=Nz/2
    do j=1,Ny
        sy = real(j,8)*dy - 0.5d0*dy 
        do i=1,Nx
            sx = real(i,8)*dx - 0.5d0*dx 
            n  = i + Nx*(j-1) + kdz*(k-1)
            s= sqrt(X(n)**2 + X(n+Ndim)**2 + X(n+2*Ndim)**2)
            
            write(n_new,'( es12.4e3,  a,  es12.4e3,  a,   es12.4e3,    a,   es12.4e3,     a,    es12.4e3,    a,    es12.4e3)') &  
                        sx,     ',  ',  sy,     ',       ',  s ,   ',  ', X(n), ',  ', X(n+Ndim), ',  ', X(n+2*Ndim)
        enddo
    enddo
    close(n_new)

    print*, 'file:' //trim(files)//'.csv'! //' saved'

    ios=0; 
    open(newunit=n_new,file=files//'.vtk',form="UNFORMATTED", access="STREAM", convert="big_endian", iostat=ios) 
    if (ios/=0) then
        print *,'Error! Could not open the outfile:'//files//'.vtk';   stop
    end if

    write (n_new) "# vtk DataFile Version 3.0"//new_line(ci)//"out data result"//new_line(ci)//"BINARY"//new_line(ci)
    write(buf1,'(i8," ",i8," ",i8)') Nx, Ny, Nz 
    write(n_new)"DATASET STRUCTURED_GRID"//new_line(ci)//"DIMENSIONS "//trim(adjustl(buf1))//new_line(ci)
    write(buf1,'(i8)')   Nx*Ny*Nz 
    write(n_new) "POINTS "//trim(adjustl(buf1))//" float"//new_line(ci)

    do k=1,Nz
        sz = real(k,8)*dz - 0.5d0*dz 
        do j=1,Ny
            sy = real(j,8)*dy - 0.5d0*dy 
            do i=1,Nx
                sx = real(i,8)*dx - 0.5d0*dx 
                write(n_new) real(sx,4), real(sy,4), real(sz,4)
            enddo
        enddo
    enddo
    write(n_new)  new_line(ci)

    write(n_new) "POINT_DATA "//trim(adjustl(buf1))//new_line(ci)

!     X is the final approximate solution
    write(n_new) "VECTORS "//"Vector_X"//" float"//new_line(ci)
    n=0
    do k=1,Nz
    do j=1,Ny
    do i=1,Nx
        n=n+1
        sx=X(n)
        sy=X(n+Ndim)
        sz=X(n+2*Ndim)
        write(n_new)  real(sx,4) ,  real(sy,4), real(sz,4)
    enddo
    enddo
    enddo
    write(n_new)  new_line(ci)

!  
    write(n_new) "VECTORS "//"Vector_Curl(X)"//" float"//new_line(ci)
    do k = 1,Nz
    do j = 1,Ny
    do i = 1,Nx
        if (k==1 .or. j==1 .or. i==1 .or.  k==Nz .or. j==Ny .or. i==Nx ) then 
            sx=0.d0; sy=0.d0; sz=0.d0
        else
            n   = i   + Nx*(j-1)  + kdz*(k-1)
            nim = n-1;  njm = n-Nx;  nkm = n-kdz;
            nip = n+1;  njp = n+Nx;  nkp = n+kdz; 

            sx   = 0.5d0*(X(2*Ndim+njp)- X(2*Ndim+njm))/dy      &
                      - 0.5d0*( X(Ndim+nkp) - X(Ndim+nkm))/dz            ! Curl_x = dXz/dy - dXy/dz 

            sy   = 0.5d0*(X(nkp) - X(nkm))/dz                     &
                      - 0.5d0*( X(2*Ndim+nip) - X(2*Ndim+nim))/dx        ! Curl_y = dXx/dz - dXz/dx 
            
            sz = 0.5d0*(X(Ndim+nip) - X(Ndim+nim))/dx           &
                      - 0.5d0*( X(njp) - X(njm))/dy                      ! Curl_z = dXy/dx - dXx/dy 
        endif
         write(n_new)  real(sx,4) ,  real(sy,4), real(sz,4)
    end do
    end do
    end do
    write(n_new)  new_line(ci)

!  Right-hand side vector (RHS)
    write(n_new) "VECTORS "//"Vector_RHS"//" float"//new_line(ci)
    n=0
    do k=1,Nz
    do j=1,Ny
    do i=1,Nx
        n=n+1
        sx=B(n)
        sy=B(n+Ndim)
        sz=B(n+2*Ndim)
        write(n_new)  real(sx,4) ,  real(sy,4), real(sz,4)
    enddo
    enddo
    enddo
    write(n_new)  new_line(ci)

    close(n_new)

    print*, 'file:' //trim(files)//'.vtk' !//' saved'

end subroutine writeCsvVtk