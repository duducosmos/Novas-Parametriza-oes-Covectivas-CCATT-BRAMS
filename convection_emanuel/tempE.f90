C
C   ***  CALCULATE TENDENCIES OF POTENTIAL TEMPERATURE AND MIXING RATIO  ***
C   ***               AT LEVELS ABOVE THE LOWEST LEVEL                   ***
C
C   ***  FIRST FIND THE NET SATURATED UPDRAFT AND DOWNDRAFT MASS FLUXES  ***
C   ***                      THROUGH EACH LEVEL                          ***
C
        DO 500 I=2,INB
        DPINV=0.01/(PH(I)-PH(I+1))
        CPINV=1.0/CPN(I)
        AMP1=0.0
        AD=0.0
        IF(I.GE.NK)THEN
         DO 440 K=I+1,INB+1
  440    AMP1=AMP1+M(K)
        END IF
        DO 450 K=1,I
        DO 450 J=I+1,INB+1
         AMP1=AMP1+MENT(K,J)
  450   CONTINUE
        IF((2.*G*DPINV*AMP1).GE.DELTI)IFLAG=4
        DO 470 K=1,I-1
        DO 470 J=I,INB
         AD=AD+MENT(J,K)
  470   CONTINUE
        FT(I)=FT(I)+G*DPINV*(AMP1*(T(I+1)-T(I)+(GZ(I+1)-GZ(I))*
     1   CPINV)-AD*(T(I)-T(I-1)+(GZ(I)-GZ(I-1))*CPINV))
     2   -SIGD*LVCP(I)*EVAP(I)
        FT(I)=FT(I)+G*DPINV*MENT(I,I)*(HP(I)-H(I)+
     1    T(I)*(CPV-CPD)*(Q(I)-QENT(I,I)))*CPINV
        FT(I)=FT(I)+SIGD*WT(I+1)*(CL-CPD)*WATER(I+1)*
     1    (T(I+1)-T(I))*DPINV*CPINV
        FQ(I)=FQ(I)+G*DPINV*(AMP1*(Q(I+1)-Q(I))-
     1    AD*(Q(I)-Q(I-1)))
        FU(I)=FU(I)+G*DPINV*(AMP1*(U(I+1)-U(I))-
     1    AD*(U(I)-U(I-1)))
        FV(I)=FV(I)+G*DPINV*(AMP1*(V(I+1)-V(I))-
     1    AD*(V(I)-V(I-1)))
        DO K=1,NTRA
         FTRA(I,K)=FTRA(I,K)+G*DPINV*(AMP1*(TRA(I+1,K)-
     1    TRA(I,K))-AD*(TRA(I,K)-TRA(I-1,K)))
        END DO
        DO 480 K=1,I-1
         AWAT=ELIJ(K,I)-(1.-EP(I))*CLW(I)
         AWAT=MAX(AWAT,0.0)
         FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-AWAT-Q(I))
         FU(I)=FU(I)+G*DPINV*MENT(K,I)*(UENT(K,I)-U(I))
         FV(I)=FV(I)+G*DPINV*MENT(K,I)*(VENT(K,I)-V(I))
         DO J=1,NTRA
          FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)-
     1     TRA(I,J))
         END DO
  480   CONTINUE
        DO 490 K=I,INB
         FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-Q(I))
         FU(I)=FU(I)+G*DPINV*MENT(K,I)*(UENT(K,I)-U(I))
         FV(I)=FV(I)+G*DPINV*MENT(K,I)*(VENT(K,I)-V(I))
         DO J=1,NTRA
          FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)-
     1     TRA(I,J))
         END DO
  490   CONTINUE
        FQ(I)=FQ(I)+SIGD*EVAP(I)+G*(MP(I+1)*
     1    (QP(I+1)-Q(I))-MP(I)*(QP(I)-Q(I-1)))*DPINV
        FU(I)=FU(I)+G*(MP(I+1)*(UP(I+1)-U(I))-MP(I)*
     1    (UP(I)-U(I-1)))*DPINV
        FV(I)=FV(I)+G*(MP(I+1)*(VP(I+1)-V(I))-MP(I)*
     1    (VP(I)-V(I-1)))*DPINV
        DO J=1,NTRA
         FTRA(I,J)=FTRA(I,J)+G*DPINV*(MP(I+1)*(TRAP(I+1,J)-TRA(I,J))-
     1    MP(I)*(TRAP(I,J)-TRA(I-1,J)))
        END DO
  500   CONTINUE
C
