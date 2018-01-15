program spin_chain
implicit none
real*8 pi,B,Jint,tend,a,IDCTR,IDCTI,fidel,temp,dt
real*8 t0,phaseR,phaseI,concu,theta
real*8,dimension(:),allocatable :: E,coeff
integer i,m,ntime,j,N
!
write(*,*) 'Interval t0 = ?'
read(*,*) ntime
!
Jint = 1.d0
B = 0.d0
tend = 4000/Jint
dt = tend/ntime
pi = 4.d0*atan(1.d0)
!========================
!Values of some variables
!========================
do N = 2, 80
    allocate(E(N))
    allocate(coeff(N))
    do m = 1, N
        theta = pi*(m-1)/N
        if (m .eq. 1) then
            a = sqrt(1.d0/N)
        else
            a = sqrt(2.d0/N)
        end if
        E(m) = 2.d0*B + 2.d0*Jint*(1.d0 - cos(theta))
        coeff(m) = a**2*cos(0.5d0*theta)*cos(theta*(2*N-1)*0.5d0)
    end do
    !------------------------------
    !Determine the optimal solution
    !------------------------------
    temp = 0.d0
    do i = 1, ntime+1
        t0 = dt*(i-1)
        IDCTR = 0.d0
        IDCTI = 0.d0
        do m = 1, N
            phaseR = cos(E(m)*t0)
            phaseI = sin(E(m)*t0)
            IDCTR = IDCTR + coeff(m)*phaseR
            IDCTI = IDCTI + coeff(m)*phaseI
        end do
        concu = sqrt(IDCTR**2 + IDCTI**2)
        if (concu .gt. temp) then
            temp = concu
        end if
    end do
    concu = temp
    fidel = concu/3.d0 + concu**2/6.d0 + 0.5d0
    !open(1,access='append',file='data_chain.txt')
    if (N .eq. 2) then
        write(*,*) '          N', '         concu','                    fidelity'
    !    write(1,*) '          N', '         concu','                    fidelity'
    end if
    write(*,*) N, concu,fidel
    !write(1,*) N, concu,fidel
    !write(1,*) 2.d0/3.d0
    !close(1)
    deallocate(E)
    deallocate(coeff)
end do
end program