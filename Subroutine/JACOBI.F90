	SUBROUTINE JACOBI(N1,N,EV,A,D,V,IROT,B,Z)
	LOGICAL EV
	DIMENSION A(N1,N1),D(N1),V(N1,N1),B(N1),Z(N1)
	IF(.NOT.EV) GOTO 10
	DO IP=1,N
	  DO IQ=1,N
	    IF(IP-IQ) 50,60,50
60	    V(IP,IP)=1.
	    GOTO 70
50	    V(IP,IQ)=0.
70      END DO
	END DO
10	DO IP=1,N
	  D(IP)=A(IP,IP)
	  B(IP)=D(IP)
	  Z(IP)=0.
	END DO
	IROT=0
	DO I=1,50
	  SM=0.
	  NM1=N-1
	  DO IP=1,NM1
	    IPP1=IP+1
	    DO IQ=IPP1,N
	      SM=SM+ABS(A(IP,IQ))
	    END DO
	  END DO
	  IF(SM) 110,120,110
110	  IF(I-4) 130,140,140
130	  TRESH=0.2*SM/(FLOAT(N)*FLOAT(N))
	  GOTO 150
140	  TRESH=0.
150	  DO IP=1,NM1
	    IPP1=IP+1
	    DO IQ=IPP1,N
	      G=100.*ABS(A(IP,IQ))
	      IF(I.GT.4.AND.ABS(D(IP))+G.EQ.ABS(D(IP))
     &	           .AND.ABS(D(IQ))+G.EQ.ABS(D(IQ))) GOTO 200
		  IF(ABS(A(IP,IQ)).LE.TRESH) GOTO 160
	      H=D(IQ)-D(IP)
	      IF(ABS(H)+G.EQ.ABS(H)) GOTO 240
	      THETA=0.5*H/A(IP,IQ)
	      T=1./(ABS(THETA)+SQRT(1.+THETA*THETA))
	      IF(THETA.LT.0.) T=-T
	      GOTO 250
240	      T=A(IP,IQ)/H
250	      C=1./SQRT(1.+T*T)
	      S=T*C
	      H=T*A(IP,IQ)
	      Z(IP)=Z(IP)-H
	      Z(IQ)=Z(IQ)+H
	      D(IP)=D(IP)-H
	      D(IQ)=D(IQ)+H
	      A(IP,IQ)=0.
	      IPM1=IP-1
	      IF(IPM1) 260,260,270
270	      DO J=1,IPM1
	        G=A(J,IP)
	        H=A(J,IQ)
	        A(J,IP)=C*G-S*H
	        A(J,IQ)=S*G+C*H
	      END DO
260	      IQM1=IQ-1
	      IF(IQM1-IPP1) 300,290,290
290	      DO J=IPP1,IQM1
	        G=A(IP,J)
	        H=A(J,IQ)
	        A(IP,J)=C*G-S*H
              A(J,IQ)=S*G+C*H
            END DO
300	      IQP1=IQ+1
	      IF(N-IQP1) 330,320,320
320	      DO J=IQP1,N
	        G=A(IP,J)
	        H=A(IQ,J)
	        A(IP,J)=C*G-S*H
	        A(IQ,J)=S*G+C*H
	      END DO
330	      IF(.NOT.EV) GOTO 350
	      DO J=1,N
	        G=V(J,IP)
	        H=V(J,IQ)
	        V(J,IP)=C*G-S*H
	        V(J,IQ)=S*G+C*H
            END DO
350	      IROT=IROT+1
	      GOTO 160
200	      A(IP,IQ)=0.
160	    END DO
	  END DO
	  DO IP=1,N
	    B(IP)=B(IP)+Z(IP)
	    D(IP)=B(IP)
	    Z(IP)=0.
	  END DO
	END DO
120	RETURN
	END
