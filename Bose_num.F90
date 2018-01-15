program main    
implicit none
logical EV
real*8 magnet,Jint,tend,IDCTR,IDCTI,fidel,concu,temp,dt
real*8 phaseR,phaseI,t0
integer IROT,i,m,ntime,j,N
real*8,dimension(:),allocatable :: E,B,Z,coeff
real*8,dimension(:,:),allocatable :: A,V
!
write(*,*) 'Interval t0 = ?'
read(*,*), ntime
!
IROT = 1
EV = .true.
Jint = 1.d0
magnet = 0.d0
tend = 4000/Jint
dt = tend/ntime
!
do N = 2, 80
    allocate(E(N))
    allocate(coeff(N))
	allocate(B(N))
	allocate(Z(N))
	allocate(A(N,N))
	allocate(V(N,N))
	B = 0.d0
	Z = 0.d0
	A = 0.d0
	V = 0.d0
	do m = 1, N-1
		A(m,m+1) = 1.d0
		A(m+1,m) = 1.d0
	end do
	!A(1,N) = 1.d0
	!A(N,1) = 1.d0
	!
	call JACOBI(N,N,EV,A,E,V,IROT,B,Z)
	!
    do m = 1, N
		E(m) = -2.d0*magnet - Jint*(E(m) - 2.d0)
        coeff(m) = V(1,m)*V(N,m)
    end do
    !------------------------------
    !Determine the optimal solution
    !------------------------------
    temp = 0.d0
    do i = 1, ntime + 1
		t0 = dt*(i - 1)
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
    !open(1,access='append',file='bosenum.txt')
    if (N .eq. 2) then
        write(*,*) '          N','     fidelity'
    !    write(1,*) '          N','     fidelity'
    end if
    write(*,*) N, fidel
    !write(1,*) N, fidel
    !close(1)
    deallocate(E)
    deallocate(coeff)
	deallocate(B)
	deallocate(Z)
	deallocate(A)
	deallocate(V)
end do
end program 


	subroutine JACOBI(N1,N,EV,A,D,V,IROT,B,Z)
	implicit none
	integer N1,N,IROT
	integer IP,IQ,NM1,I,IPP1,IPM1,IQM1,J,IQP1
	real*8 SM,TRESH,G,H,THETA,T,C,S
	real*8 A(N1,N1),D(N1),V(N1,N1),B(N1),Z(N1)
	logical EV
	!dimension A(N1,N1),D(N1),V(N1,N1),B(N1),Z(N1)
	if (.not. EV) goto 10
	do IP = 1, N
		do IQ = 1, N
	    	if (IP - IQ) 50,60,50
60	    	V(IP,IP) = 1.d0
	    	goto 70
50	    	V(IP,IQ) = 0.d0
70  	end do
	end do
10	do IP = 1, N
	  	D(IP) = A(IP,IP)
	  	B(IP) = D(IP)
	  	Z(IP) = 0.d0
	end do
	IROT = 0
	do I = 1, 50
	  	SM = 0.d0
	  	NM1 = N - 1
	  	do IP = 1, NM1
	  	  	IPP1 = IP + 1
	  	  	do IQ = IPP1, N
	  	  	  	SM = SM + abs(A(IP,IQ))
			end do
		end do
	  	if (SM) 110,120,110
110	  	if (I - 4) 130,140,140
130	  	TRESH = 0.2d0*SM/(N*N)
	  	goto 150
140	  	TRESH = 0.d0
150	  	do IP = 1, NM1
	  		IPP1 = IP + 1
	  		do IQ = IPP1, N
	  		  	G = 100.d0*abs(A(IP,IQ))
	  		  	if (I .gt. 4 .and. abs(D(IP)) + G .eq. abs(D(IP)) &
      		  	         .and. abs(D(IQ)) + G .eq. abs(D(IQ))) goto 200
	  		  	if (abs(A(IP,IQ)) .le. TRESH) goto 160
	  		  	H = D(IQ) - D(IP)
	  		  	if (abs(H) + G .eq. abs(H)) goto 240
	  		  	THETA = 0.5*H/A(IP,IQ)
	  		  	T = 1.d0/(abs(THETA) + sqrt(1.d0 + THETA*THETA))
	  		  	if (THETA .lt. 0.d0) T = -T
	  		  	goto 250
240	    	 	T = A(IP,IQ)/H
250	    	 	C = 1.d0/sqrt(1.d0 + T*T)
	  		  	S = T*C
	  		  	H = T*A(IP,IQ)
	  		  	Z(IP) = Z(IP) - H
	  		  	Z(IQ) = Z(IQ) + H
	    	  	D(IP) = D(IP) - H
	    	  	D(IQ) = D(IQ) + H
	    	  	A(IP,IQ) = 0.d0
	    	  	IPM1 = IP - 1
	    	  	if (IPM1) 260,260,270
270	    	  	do J = 1, IPM1
	    	  	  	G = A(J,IP)
	    	  	  	H = A(J,IQ)
	    	  	  	A(J,IP) = C*G - S*H
	    	  	  	A(J,IQ) = S*G + C*H
	    	  	end do
260	    	  	IQM1 = IQ - 1
	    	  	if (IQM1 - IPP1) 300,290,290
290	    	  	do J = IPP1, IQM1
	    	    	G = A(IP,J)
	    	    	H = A(J,IQ)
	    	    	A(IP,J) = C*G - S*H
        	    	A(J,IQ) = S*G + C*H
        	    end do
300	      		IQP1 = IQ + 1
	      		if (N - IQP1) 330,320,320
320	      		do J = IQP1, N
	      		  	G = A(IP,J)
	      		  	H = A(IQ,J)
	      		  	A(IP,J) = C*G - S*H
	      		  	A(IQ,J) = S*G + C*H
	      		end do
330	      		if (.not. EV) goto 350
	      		do J = 1, N
	      		  	G = V(J,IP)
	      		  	H = V(J,IQ)
	      		  	V(J,IP) = C*G - S*H
	      		  	V(J,IQ) = S*G + C*H
          		end do
350	      		IROT = IROT + 1
	      		goto 160
200	      		A(IP,IQ) = 0.d0
160	    	end do
		end do
		do IP = 1, N
		  	B(IP) = B(IP) + Z(IP)
		  	D(IP) = B(IP)
		  	Z(IP) = 0.d0
		end do
	end do
120	return
	end subroutine JACOBI