C
C   ***  CHECK WHETHER EP(INB)=0, IF SO, SKIP PRECIPITATING    ***
C   ***             DOWNDRAFT CALCULATION                      ***
C
        IF(EP(INB).LT.0.0001)GOTO 405
C
C   ***  INTEGRATE LIQUID WATER EQUATION TO FIND CONDENSED WATER   ***
C   ***                AND CONDENSED WATER FLUX                    ***
C
        JTT=2
C
C    ***                    BEGIN DOWNDRAFT LOOP                    ***
C
        DO 400 I=INB,1,-1
C
C    ***              CALCULATE DETRAINED PRECIPITATION             ***
C
        WDTRAIN=G*EP(I)*M(I)*CLW(I)
        IF(I.GT.1)THEN
         DO 320 J=1,I-1
         AWAT=ELIJ(J,I)-(1.-EP(I))*CLW(I)
         AWAT=MAX(0.0,AWAT)
  320    WDTRAIN=WDTRAIN+G*AWAT*MENT(J,I)
        END IF
C
C    ***    FIND RAIN WATER AND EVAPORATION USING PROVISIONAL   ***
C    ***              ESTIMATES OF QP(I)AND QP(I-1)             ***
C     
c
c  ***  Value of terminal velocity and coefficient of evaporation for snow   ***
c 
        COEFF=COEFFS
        WT(I)=OMTSNOW
c      
c  ***  Value of terminal velocity and coefficient of evaporation for rain   ***
c
        IF(T(I).GT.273.0)THEN
         COEFF=COEFFR
         WT(I)=OMTRAIN
        END IF
        QSM=0.5*(Q(I)+QP(I+1))
        AFAC=COEFF*PH(I)*(QS(I)-QSM)/(1.0E4+2.0E3*PH(I)*QS(I))
        AFAC=MAX(AFAC,0.0)
        SIGT=SIGP(I)
        SIGT=MAX(0.0,SIGT)
        SIGT=MIN(1.0,SIGT)
        B6=100.*(PH(I)-PH(I+1))*SIGT*AFAC/WT(I)
        C6=(WATER(I+1)*WT(I+1)+WDTRAIN/SIGD)/WT(I)
        REVAP=0.5*(-B6+SQRT(B6*B6+4.*C6))
        EVAP(I)=SIGT*AFAC*REVAP
        WATER(I)=REVAP*REVAP
C
C    ***  CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER     ***
C    ***              HYDROSTATIC APPROXIMATION                 ***
C   
        IF(I.EQ.1)GOTO 360
        DHDP=(H(I)-H(I-1))/(P(I-1)-P(I))
        DHDP=MAX(DHDP,10.0)
        MP(I)=100.*GINV*LV(I)*SIGD*EVAP(I)/DHDP
        MP(I)=MAX(MP(I),0.0)
C
C   ***   ADD SMALL AMOUNT OF INERTIA TO DOWNDRAFT              ***
C
        FAC=20.0/(PH(I-1)-PH(I))
        MP(I)=(FAC*MP(I+1)+MP(I))/(1.+FAC)
C   
C    ***      FORCE MP TO DECREASE LINEARLY TO ZERO                 ***
C    ***      BETWEEN ABOUT 950 MB AND THE SURFACE                  ***
C
          IF(P(I).GT.(0.949*P(1)))THEN
           JTT=MAX(JTT,I)
           MP(I)=MP(JTT)*(P(1)-P(I))/(P(1)-P(JTT))
          END IF              
  360   CONTINUE
C
C    ***       FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT     ***
C
        IF(I.EQ.INB)GOTO 400
        IF(I.EQ.1)THEN
         QSTM=QS(1)
        ELSE
         QSTM=QS(I-1)
        END IF
        IF(MP(I).GT.MP(I+1))THEN
          RAT=MP(I+1)/MP(I)
          QP(I)=QP(I+1)*RAT+Q(I)*(1.0-RAT)+100.*GINV*
     1       SIGD*(PH(I)-PH(I+1))*(EVAP(I)/MP(I))
          UP(I)=UP(I+1)*RAT+U(I)*(1.-RAT)
          VP(I)=VP(I+1)*RAT+V(I)*(1.-RAT)
          DO J=1,NTRA
           TRAP(I,J)=TRAP(I+1,J)*RAT+TRAP(I,J)*(1.-RAT)
          END DO
         ELSE
          IF(MP(I+1).GT.0.0)THEN
            QP(I)=(GZ(I+1)-GZ(I)+QP(I+1)*(LV(I+1)+T(I+1)*(
     1        CL-CPD))+CPD*(T(I+1)-T(I)))/(LV(I)+T(I)*(CL-CPD))
            UP(I)=UP(I+1)
            VP(I)=VP(I+1)
            DO J=1,NTRA
             TRAP(I,J)=TRAP(I+1,J)
            END DO
          END IF
        END IF
        QP(I)=MIN(QP(I),QSTM)
        QP(I)=MAX(QP(I),0.0)
  400   CONTINUE
  
  C
C   ***  CALCULATE SURFACE PRECIPITATION IN MM/DAY     ***
C
        PRECIP=PRECIP+WT(1)*SIGD*WATER(1)*3600.*24000./(ROWL*G)
C
  405   CONTINUE
  
  

C
C   ***  CALCULATE DOWNDRAFT VELOCITY SCALE AND SURFACE TEMPERATURE AND  ***
c   ***                    WATER VAPOR FLUCTUATIONS                      ***
C
      WD=BETA*ABS(MP(ICB))*0.01*RD*T(ICB)/(SIGD*P(ICB))
      QPRIME=0.5*(QP(1)-Q(1))
      TPRIME=LV0*QPRIME/CPD
C
C   ***  CALCULATE TENDENCIES OF LOWEST LEVEL POTENTIAL TEMPERATURE  ***
C   ***                      AND MIXING RATIO                        ***
C
        DPINV=0.01/(PH(1)-PH(2))
        AM=0.0
        IF(NK.EQ.1)THEN
         DO 410 K=2,INB
  410    AM=AM+M(K)
        END IF
        IF((2.*G*DPINV*AM).GE.DELTI)IFLAG=4
        FT(1)=FT(1)+G*DPINV*AM*(T(2)-T(1)+(GZ(2)-GZ(1))/CPN(1))
        FT(1)=FT(1)-LVCP(1)*SIGD*EVAP(1)
        FT(1)=FT(1)+SIGD*WT(2)*(CL-CPD)*WATER(2)*(T(2)-
     1   T(1))*DPINV/CPN(1)
        FQ(1)=FQ(1)+G*MP(2)*(QP(2)-Q(1))*
     1    DPINV+SIGD*EVAP(1)
        FQ(1)=FQ(1)+G*AM*(Q(2)-Q(1))*DPINV
        FU(1)=FU(1)+G*DPINV*(MP(2)*(UP(2)-U(1))+AM*(U(2)-U(1)))
        FV(1)=FV(1)+G*DPINV*(MP(2)*(VP(2)-V(1))+AM*(V(2)-V(1)))
        DO J=1,NTRA
         FTRA(1,J)=FTRA(1,J)+G*DPINV*(MP(2)*(TRAP(2,J)-TRA(1,J))+
     1    AM*(TRA(2,J)-TRA(1,J)))
        END DO
        AMDE=0.0
        DO 415 J=2,INB
         FQ(1)=FQ(1)+G*DPINV*MENT(J,1)*(QENT(J,1)-Q(1))
         FU(1)=FU(1)+G*DPINV*MENT(J,1)*(UENT(J,1)-U(1))
         FV(1)=FV(1)+G*DPINV*MENT(J,1)*(VENT(J,1)-V(1))
         DO K=1,NTRA
          FTRA(1,K)=FTRA(1,K)+G*DPINV*MENT(J,1)*(TRAENT(J,1,K)-
     1     TRA(1,K))
         END DO
  415   CONTINUE
