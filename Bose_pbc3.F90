!The same as Bose_pbc2.F90
!Exact Solution derived from Wang's paper
program main
implicit none
real*8 pi,Jint,temp,t0,ampR,ampI,theta,fidel,t1,t2
real*8 dt,h,temp1
integer ntime,hgap,N,j1,k1,i,m,ii
real*8,dimension(:),allocatable :: E
!
write(*,*) 'tgap = ?'
read(*,*) ntime
!
write(*,*) 'hgap = ?'
read(*,*) hgap
!
call cpu_time(t1)
Jint = 1.d0
pi = 4.d0*atan(1.d0)
dt = 4000.d0/ntime
j1 = 1
!
do N = 3, 50
    allocate(E(N))
    !
    if (mod(N,2) .eq. 1) then
        k1 = (N + 1)/2
    else
        k1 = N/2 +1
    end if
    !
    temp1 = 0.d0
    do ii = 1, hgap + 1
        h = 2.d0/hgap*(ii - 1)
        do i = 1, N
            E(i) = 2.d0*h - 2.d0*Jint*cos(2.d0*pi*i/N)
        end do
        !==============================
        !Determine the optimal solution
        !==============================
        temp = 0.d0
        do m = 1, ntime + 1
            t0 = dt*(m - 1)
            ampR = 0.d0
            ampI = 0.d0
            do i = 1, N
                theta = 2.d0*pi*i/N
                ampR = ampR + 1.d0/N*cos(E(i)*t0 - theta*(j1 - k1))
                ampI = ampI + 1.d0/N*sin(E(i)*t0 - theta*(j1 - k1))
            end do
            fidel = ampR/3.d0 + (ampR**2 + ampI**2)/6.d0 +0.5d0
            if(fidel .gt. temp) then
                temp = fidel
            end if
        end do
        if (temp .gt. temp1) then
            temp1 = temp
        end if
    end do
    fidel = temp1
    !open(1,access='append',file='dat_pbc3.txt')
    if (N .eq. 3) then
        write(*,*) '          N','         fidelity'
        !write(1,*) '          N','         fidelity'
    end if
    write(*,*) N,fidel
    !write(1,*) N,fidel
    !close(1)
    deallocate(E)
end do
!
call cpu_time(t2)
write(*,*) 'Time = ', t2-t1
end program