!   ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***
!   ***                    VIRTUAL TEMPERATURE                    ***
!
        do_sixty_four: DO I=ICB+1,NL
             TVP(I)=TVP(I)-TP(I)*Q(NK)
        end do_sixty_four
        
        TVP(NL+1)=TVP(NL)-(GZ(NL+1)-GZ(NL))/CPD
!
!   ***        NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS       ***
!
        !DO 70 I=1,NL+1
         HP(1:NL+1)=H(1:NL+1);         NENT(1:NL+1)=0;          WATER(1:NL+1)=0.0;         EVAP(1:NL+1)=0.0
         WT(1:NL+1)=OMTSNOW;         MP(1:NL+1)=0.0;         M(1:NL+1)=0.0;         LVCP(1:NL+1)=LV(1:NL+1)/CPN(1:NL+1)
         !DO 70 J=1,NL+1
         QENT(1:NL+1,1:NL+1)=Q(1:NL+1);         ELIJ(1:NL+1,1:NL+1)=0.0;         MENT(1:NL+1,1:NL+1)=0.0
         SIJ(1:NL+1,1:NL+1)=0.0;         UENT(1:NL+1,1:NL+1)=U(1:NL+1);         VENT(1:NL+1,1:NL+1)=V(1:NL+1)
         !DO 70 K=1,NTRA
         TRAENT(1:NL+1,1:NL+1,1:NTRA)=TRA(1:NL+1,NTRA)
   !70   CONTINUE
   
   
   
        QP(1)=Q(1); UP(1)=U(1);  VP(1)=V(1);
        
        !DO 71 I=1,NTRA
        TRAP(1,1:NTRA)=TRA(1,1:NTRA)
   !71	CONTINUE
        do_seventy_tow: DO I=2,NL+1
            QP(I)=Q(I-1)
            UP(I)=U(I-1)
            VP(I)=V(I-1)
            !DO  J=1,NTRA
            TRAP(I,1:NTRA)=TRA(I-1,1:NTRA)
        enddo do_seventy_tow
!
!  ***  FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S      ***
!  ***          HIGHEST LEVEL OF NEUTRAL BUOYANCY                 ***
!  ***     AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)           ***
!
        CAPE=0.0
        CAPEM=0.0
        INB=ICB+1
        INB1=INB
	    BYP=0.0
        do_eighty_tow: DO  I=ICB+1,NL-1
             BY=(TVP(I)-TV(I))*(PH(I)-PH(I+1))/P(I)
             CAPE=CAPE+BY
             IF(BY.GE.0.0)INB1=I+1
             IF(CAPE.GT.0.0)THEN
                 INB=I+1
                 BYP=(TVP(I+1)-TV(I+1))*(PH(I+1)-PH(I+2))/P(I+1)
                 CAPEM=CAPE
             END IF
        end do do_eighty_tow
        INB=MAX(INB,INB1)
        CAPE=CAPEM+BYP
        DEFRAC=CAPEM-CAPE
        DEFRAC=MAX(DEFRAC,0.001)
        FRAC=-CAPE/DEFRAC
        FRAC=MIN(FRAC,1.0)
        FRAC=MAX(FRAC,0.0)
!
!   ***   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL   ***
!
        !DO 95 I=ICB,INB
         HP(ICB:INB)=H(NK)+(LV(ICB:INB)+(CPD-CPV)*T(ICB:INB))*EP(ICB:INB)*CLW(ICB:INB)
   !95   CONTINUE                  
!
!   ***  CALCULATE CLOUD BASE MASS FLUX AND RATES OF MIXING, M(I),  ***
!   ***                   AT EACH MODEL LEVEL                       ***
!
        DBOSUM=0.0
!   
!   ***     INTERPOLATE DIFFERENCE BETWEEN LIFTED PARCEL AND      ***
!   ***  ENVIRONMENTAL TEMPERATURES TO LIFTED CONDENSATION LEVEL  ***
!	
        TVPPLCL=TVP(ICB-1)-RD*TVP(ICB-1)*(P(ICB-1)-PLCL)/(CPN(ICB-1)*P(ICB-1))
        TVAPLCL=TV(ICB)+(TVP(ICB)-TVP(ICB+1))*(PLCL-P(ICB))/(P(ICB)-P(ICB+1))
        DTPBL=0.0
        
        do_ninety_six: DO  I=NK,ICB-1
            DTPBL=DTPBL+(TVP(I)-TV(I))*(PH(I)-PH(I+1))
        end do do_ninety_six
        
        DTPBL=DTPBL/(PH(NK)-PH(ICB))
        DTMIN=TVPPLCL-TVAPLCL+DTMAX+DTPBL
        DTMA=DTMIN
!
!   ***  ADJUST CLOUD BASE MASS FLUX   ***
!
      CBMFOLD=CBMF
	  DELT0=300.0
      DAMPS=DAMP*DELT/DELT0 
      CBMF=(1.-DAMPS)*CBMF+0.1*ALPHA*DTMA 
      CBMF=MAX(CBMF,0.0)
!
!   *** If cloud base mass flux is zero, skip rest of calculation  ***
!
      IF(CBMF.EQ. 0.0 .AND. CBMFOLD.EQ.0.0)THEN
          RETURN
      END IF
!
!   ***   CALCULATE RATES OF MIXING,  M(I)   ***
!
      M(ICB)=0.0
      do_a_hundred_three: DO I=ICB+1,INB
          K=MIN(I,INB1)
          DBO=ABS(TV(K)-TVP(K))+ENTP*0.02*(PH(K)-PH(K+1))
          DBOSUM=DBOSUM+DBO
          M(I)=CBMF*DBO
      do_a_hundred_three
      
      !DO 110 I=ICB+1,INB
      M(ICB+1:INB)=M(ICB+1:INB)/DBOSUM  
  !110 CONTINUE     
!
!   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
!   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
!   ***                        FRACTION (SIJ)                          ***
!
         do_a_hundred_seventy: DO 170 I=ICB+1,INB
             QTI=Q(NK)-EP(I)*CLW(I)
         	 do_a_hundred_sixty: DO  J=ICB,INB
                 BF2=1.+LV(J)*LV(J)*QS(J)/(RV*T(J)*T(J)*CPD)
                 ANUM=H(J)-HP(I)+(CPV-CPD)*T(J)*(QTI-Q(J))
                 DENOM=H(I)-HP(I)+(CPD-CPV)*(Q(I)-QTI)*T(J)
                 DEI=DENOM
                 IF(ABS(DEI).LT.0.01)DEI=0.01
			         SIJ(I,J)=ANUM/DEI
			         SIJ(I,I)=1.0
			         ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)
			         ALTEM=ALTEM/BF2
			         CWAT=CLW(J)*(1.-EP(J))
			         STEMP=SIJ(I,J)
			         IF((STEMP.LT.0.0 .OR. STEMP .GT. 1.0 .OR. ALTEM.GT.CWAT).AND.J.GT.I)THEN
			             ANUM=ANUM-LV(J)*(QTI-QS(J)-CWAT*BF2)
			             DENOM=DENOM+LV(J)*(Q(I)-QTI)
			             IF(ABS(DENOM).LT.0.01)DENOM=0.01
			             SIJ(I,J)=ANUM/DENOM
			             ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)
			             ALTEM=ALTEM-(BF2-1.)*CWAT
			         END IF
			         IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
			             QENT(I,J)=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI
			             UENT(I,J)=SIJ(I,J)*U(I)+(1.-SIJ(I,J))*U(NK)
			             VENT(I,J)=SIJ(I,J)*V(I)+(1.-SIJ(I,J))*V(NK)
			         !DO K=1,NTRA
			             TRAENT(I,J,1:NTRA)=SIJ(I,J)*TRA(I,1:NTRA)+(1.-SIJ(I,J))*TRA(NK,1:NTRA)
			         !END DO
			             ELIJ(I,J)=ALTEM
			             ELIJ(I,J)=MAX(0.0,ELIJ(I,J))
			             MENT(I,J)=M(I)/(1.-SIJ(I,J))
			             NENT(I)=NENT(I)+1
			         END IF
			         SIJ(I,J)=MAX(0.0,SIJ(I,J))
			         SIJ(I,J)=MIN(1.0,SIJ(I,J))
             end do do_a_hundred_sixty
!
!   ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS  ***
!   ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES  ***
!
             IF(NENT(I).EQ.0)THEN
                 MENT(I,I)=M(I)
                 QENT(I,I)=Q(NK)-EP(I)*CLW(I)
                 UENT(I,I)=U(NK)
                 VENT(I,I)=V(NK)
          !DO J=1,NTRA
                 TRAENT(I,I,1:NTRA)=TRA(NK,1:NTRA)
          !END DO
                 ELIJ(I,I)=CLW(I)
                 SIJ(I,I)=1.0
             END IF 
         end do do_a_hundred_seventy
         
         SIJ(INB,INB)=1.0
!
!   ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***
!   ***              PROBABILITIES OF MIXING                     ***
!
        do_tow_hundred: DO 200 I=ICB+1,INB
            IF(NENT(I).NE.0)THEN
                QP1=Q(NK)-EP(I)*CLW(I)
                ANUM=H(I)-HP(I)-LV(I)*(QP1-QS(I))
                DENOM=H(I)-HP(I)+LV(I)*(Q(I)-QP1)
                IF(ABS(DENOM).LT.0.01)DENOM=0.01
                SCRIT=ANUM/DENOM
                ALT=QP1-QS(I)+SCRIT*(Q(I)-QP1)
                IF(ALT.LT.0.0)SCRIT=1.0
                SCRIT=MAX(SCRIT,0.0)
                ASIJ=0.0
                SMIN=1.0
                DO 175 J=ICB,INB
                IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
                    IF(J.GT.I)THEN
                        SMID=MIN(SIJ(I,J),SCRIT)
                        SJMAX=SMID
                        SJMIN=SMID
                        IF(SMID.LT.SMIN.AND.SIJ(I,J+1).LT.SMID)THEN
                            SMIN=SMID
                            SJMAX=MIN(SIJ(I,J+1),SIJ(I,J),SCRIT)
                            SJMIN=MAX(SIJ(I,J-1),SIJ(I,J))
                            SJMIN=MIN(SJMIN,SCRIT)
                        END IF
                    ELSE
                        SJMAX=MAX(SIJ(I,J+1),SCRIT)
                        SMID=MAX(SIJ(I,J),SCRIT)
                        SJMIN=0.0
                        IF(J.GT.1)SJMIN=SIJ(I,J-1)
                        SJMIN=MAX(SJMIN,SCRIT)
                    END IF
                    DELP=ABS(SJMAX-SMID)
                               
                    !   ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***
                    !   ***                    VIRTUAL TEMPERATURE                    ***
                    !
                    !DO 64 I=ICB+1,NL
                    TVP(ICB+1:NL)=TVP(ICB+1:NL)-TP(ICB+1:NL)*Q(NK)
                    !64   CONTINUE
                    TVP(NL+1)=TVP(NL)-(GZ(NL+1)-GZ(NL))/CPD
                    !
                    !   ***        NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS       ***
                    !
                    !DO 70 I=1,NL+1
                    HP(1:NL+1)=H(1:NL+1)
                    NENT(1:NL+1)=0
                    WATER(1:NL+1)=0.0
                    EVAP(1:NL+1)=0.0
                    WT(1:NL+1)=OMTSNOW
                    MP(1:NL+1)=0.0
                    M(1:NL+1)=0.0
                    LVCP(1:NL+1)=LV(1:NL+1)/CPN(1:NL+1)
                    ELIJ(1:NL+1,1:NL+1)=0.0
                    MENT(1:NL+1,1:NL+1)=0.0
                    SIJ(1:NL+1,1:NL+1)=0.0
                    DO  J=1,NL+1
                         QENT(1:NL+1,J)=Q(J)
                         ELIJ(1:NL+1,J)=0.0
                         MENT(1:NL+1,J)=0.0
                         SIJ(1:NL+1,J)=0.0
                         UENT(1:NL+1,J)=U(J)
                         VENT(1:NL+1,J)=V(J)
                         TRAENT(1:NL+1,J,1:NTRA)=TRA(J,1:NTRA)
                    enddo
                  
                  
                  
                  
                    QP(1)=Q(1)
                    UP(1)=U(1)
                    VP(1)=V(1)
                    !DO 71 I=1,NTRA
                    TRAP(1,1:NTRA)=TRA(1,1:NTRA)
                    !71	CONTINUE
                    !DO 72 I=2,NL+1
                    QP(2:NL+1)=Q(1:NL)
                    UP(2:NL+1)=U(1:NL)
                    VP(2:NL+1)=V(1:NL)
                    !DO 72 J=1,NTRA
                    TRAP(2:NL+1,1:NTRA)=TRA(1:NL,1:NTRA)
                    !72	   CONTINUE
!
!  ***  FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S      ***
!  ***          HIGHEST LEVEL OF NEUTRAL BUOYANCY                 ***
!  ***     AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)           ***
!
          
                    CAPE=0.0
                    CAPEM=0.0
                    INB=ICB+1
                    INB1=INB
	                BYP=0.0
                    eight_tow: DO I=ICB+1,NL-1
                        BY=(TVP(I)-TV(I))*(PH(I)-PH(I+1))/P(I)
                        CAPE=CAPE+BY
                        IF(BY.GE.0.0)INB1=I+1
                        IF(CAPE.GT.0.0)THEN
                            INB=I+1
                            BYP=(TVP(I+1)-TV(I+1))*(PH(I+1)-PH(I+2))/P(I+1)
                            CAPEM=CAPE
                        END IF
                    enddo eight_tow
                    INB=MAX(INB,INB1)
                    CAPE=CAPEM+BYP
                    DEFRAC=CAPEM-CAPE
                    DEFRAC=MAX(DEFRAC,0.001)
                    FRAC=-CAPE/DEFRAC
                    FRAC=MIN(FRAC,1.0)
                    FRAC=MAX(FRAC,0.0)
!
!   ***   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL   ***
!
                    do_ninety_five: DO 95 I=ICB,INB
                        HP(I)=H(NK)+(LV(I)+(CPD-CPV)*T(I))*EP(I)*CLW(I)
                    end do do_ninety_five                  
!
!   ***  CALCULATE CLOUD BASE MASS FLUX AND RATES OF MIXING, M(I),  ***
!   ***                   AT EACH MODEL LEVEL                       ***
!
                    DBOSUM=0.0
!   
!   ***     INTERPOLATE DIFFERENCE BETWEEN LIFTED PARCEL AND      ***
!   ***  ENVIRONMENTAL TEMPERATURES TO LIFTED CONDENSATION LEVEL  ***
!	
                    TVPPLCL=TVP(ICB-1)-RD*TVP(ICB-1)*(P(ICB-1)-PLCL)/(CPN(ICB-1)*P(ICB-1))
                    TVAPLCL=TV(ICB)+(TVP(ICB)-TVP(ICB+1))*(PLCL-P(ICB))/(P(ICB)-P(ICB+1))
                    DTPBL=0.0
                    DO 96 I=NK,ICB-1
                        DTPBL=DTPBL+(TVP(I)-TV(I))*(PH(I)-PH(I+1))
                    96   CONTINUE
                    DTPBL=DTPBL/(PH(NK)-PH(ICB))
                    DTMIN=TVPPLCL-TVAPLCL+DTMAX+DTPBL
                    DTMA=DTMIN
!
!   ***  ADJUST CLOUD BASE MASS FLUX   ***
!
                    CBMFOLD=CBMF
                	DELT0=300.0
                    DAMPS=DAMP*DELT/DELT0 
                    CBMF=(1.-DAMPS)*CBMF+0.1*ALPHA*DTMA 
                    CBMF=MAX(CBMF,0.0)
!
!   *** If cloud base mass flux is zero, skip rest of calculation  ***
!
                    IF(CBMF.EQ.0.0.AND.CBMFOLD.EQ.0.0)THEN
                        RETURN
                    END IF
!
!   ***   CALCULATE RATES OF MIXING,  M(I)   ***
!
                    M(ICB)=0.0
                    do_a_hundred_three: DO I=ICB+1,INB
                        K=MIN(I,INB1)
                        DBO=ABS(TV(K)-TVP(K))+ENTP*0.02*(PH(K)-PH(K+1))
                        DBOSUM=DBOSUM+DBO
                        M(I)=CBMF*DBO
                    enddo do_a_hundred_three
                    do_a_hundred_ten:DO 110 I=ICB+1,INB
                        M(I)=M(I)/DBOSUM  
                    enddo do_a_hundred_ten     
C
C   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
C   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
C   ***                        FRACTION (SIJ)                          ***
C
        DO 170 I=ICB+1,INB
         QTI=Q(NK)-EP(I)*CLW(I)
         DO 160 J=ICB,INB
          BF2=1.+LV(J)*LV(J)*QS(J)/(RV*T(J)*T(J)*CPD)
          ANUM=H(J)-HP(I)+(CPV-CPD)*T(J)*(QTI-Q(J))
          DENOM=H(I)-HP(I)+(CPD-CPV)*(Q(I)-QTI)*T(J)
          DEI=DENOM
          IF(ABS(DEI).LT.0.01)DEI=0.01
          SIJ(I,J)=ANUM/DEI
          SIJ(I,I)=1.0
          ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)
          ALTEM=ALTEM/BF2
          CWAT=CLW(J)*(1.-EP(J))
          STEMP=SIJ(I,J)
          IF((STEMP.LT.0.0.OR.STEMP.GT.1.0.OR.
     1      ALTEM.GT.CWAT).AND.J.GT.I)THEN
           ANUM=ANUM-LV(J)*(QTI-QS(J)-CWAT*BF2)
           DENOM=DENOM+LV(J)*(Q(I)-QTI)
           IF(ABS(DENOM).LT.0.01)DENOM=0.01
           SIJ(I,J)=ANUM/DENOM
           ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)
           ALTEM=ALTEM-(BF2-1.)*CWAT
          END IF
          IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
           QENT(I,J)=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI
           UENT(I,J)=SIJ(I,J)*U(I)+(1.-SIJ(I,J))*U(NK)
           VENT(I,J)=SIJ(I,J)*V(I)+(1.-SIJ(I,J))*V(NK)
           DO K=1,NTRA
            TRAENT(I,J,K)=SIJ(I,J)*TRA(I,K)+(1.-SIJ(I,J))*
     1       TRA(NK,K)
           END DO
           ELIJ(I,J)=ALTEM
           ELIJ(I,J)=MAX(0.0,ELIJ(I,J))
           MENT(I,J)=M(I)/(1.-SIJ(I,J))
           NENT(I)=NENT(I)+1
          END IF
          SIJ(I,J)=MAX(0.0,SIJ(I,J))
          SIJ(I,J)=MIN(1.0,SIJ(I,J))
  160    CONTINUE
C
C   ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS  ***
C   ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES  ***
C
         IF(NENT(I).EQ.0)THEN
          MENT(I,I)=M(I)
          QENT(I,I)=Q(NK)-EP(I)*CLW(I)
          UENT(I,I)=U(NK)
          VENT(I,I)=V(NK)
          DO J=1,NTRA
           TRAENT(I,I,J)=TRA(NK,J)
          END DO
          ELIJ(I,I)=CLW(I)
          SIJ(I,I)=1.0
         END IF 
  170   CONTINUE
        SIJ(INB,INB)=1.0
C
C   ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***
C   ***              PROBABILITIES OF MIXING                     ***
C
        DO 200 I=ICB+1,INB
        IF(NENT(I).NE.0)THEN
         QP1=Q(NK)-EP(I)*CLW(I)
         ANUM=H(I)-HP(I)-LV(I)*(QP1-QS(I))
         DENOM=H(I)-HP(I)+LV(I)*(Q(I)-QP1)
         IF(ABS(DENOM).LT.0.01)DENOM=0.01
         SCRIT=ANUM/DENOM
         ALT=QP1-QS(I)+SCRIT*(Q(I)-QP1)
         IF(ALT.LT.0.0)SCRIT=1.0
	   SCRIT=MAX(SCRIT,0.0)
         ASIJ=0.0
         SMIN=1.0
         DO 175 J=ICB,INB
          IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
           IF(J.GT.I)THEN
            SMID=MIN(SIJ(I,J),SCRIT)
            SJMAX=SMID
            SJMIN=SMID
            IF(SMID.LT.SMIN.AND.SIJ(I,J+1).LT.SMID)THEN
             SMIN=SMID
             SJMAX=MIN(SIJ(I,J+1),SIJ(I,J),SCRIT)
             SJMIN=MAX(SIJ(I,J-1),SIJ(I,J))
             SJMIN=MIN(SJMIN,SCRIT)
            END IF
           ELSE
            SJMAX=MAX(SIJ(I,J+1),SCRIT)
            SMID=MAX(SIJ(I,J),SCRIT)
            SJMIN=0.0
            IF(J.GT.1)SJMIN=SIJ(I,J-1)
            SJMIN=MAX(SJMIN,SCRIT)
           END IF
           DELP=ABS(SJMAX-SMID)
           DELM=ABS(SJMIN-SMID)
           ASIJ=ASIJ+(DELP+DELM)*(PH(J)-PH(J+1))
           MENT(I,J)=MENT(I,J)*(DELP+DELM)*(PH(J)-PH(J+1))
          END IF
  175    CONTINUE
         ASIJ=MAX(1.0E-21,ASIJ)
         ASIJ=1.0/ASIJ
         DO 180 J=ICB,INB
          MENT(I,J)=MENT(I,J)*ASIJ
  180    CONTINUE
         BSUM=0.0
         DO 190 J=ICB,INB
          BSUM=BSUM+MENT(I,J)
  190    CONTINUE
         IF(BSUM.LT.1.0E-18)THEN
          NENT(I)=0
          MENT(I,I)=M(I)
          QENT(I,I)=Q(NK)-EP(I)*CLW(I)
          UENT(I,I)=U(NK)
          VENT(I,I)=V(NK)
          DO J=1,NTRA
           TRAENT(I,I,J)=TRA(NK,J)
          END DO
          ELIJ(I,I)=CLW(I)
          SIJ(I,I)=1.0
         END IF
        END IF
        end do do_tow_hundred
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
C   *** Adjust tendencies at top of convection layer to reflect  ***
C   ***       actual position of the level zero CAPE             ***
C
        FQOLD=FQ(INB)
        FQ(INB)=FQ(INB)*(1.-FRAC)
        FQ(INB-1)=FQ(INB-1)+FRAC*FQOLD*((PH(INB)-PH(INB+1))/
     1   (PH(INB-1)-PH(INB)))*LV(INB)/LV(INB-1)
        FTOLD=FT(INB)
        FT(INB)=FT(INB)*(1.-FRAC)
        FT(INB-1)=FT(INB-1)+FRAC*FTOLD*((PH(INB)-PH(INB+1))/
     1   (PH(INB-1)-PH(INB)))*CPN(INB)/CPN(INB-1)
        FUOLD=FU(INB)
        FU(INB)=FU(INB)*(1.-FRAC)
        FU(INB-1)=FU(INB-1)+FRAC*FUOLD*((PH(INB)-PH(INB+1))/
     1   (PH(INB-1)-PH(INB)))
        FVOLD=FV(INB)
        FV(INB)=FV(INB)*(1.-FRAC)
        FV(INB-1)=FV(INB-1)+FRAC*FVOLD*((PH(INB)-PH(INB+1))/
     1   (PH(INB-1)-PH(INB)))
        DO K=1,NTRA
         FTRAOLD=FTRA(INB,K)
         FTRA(INB,K)=FTRA(INB,K)*(1.-FRAC)
         FTRA(INB-1,K)=FTRA(INB-1,K)+FRAC*FTRAOLD*(PH(INB)-PH(INB+1))/
     1    (PH(INB-1)-PH(INB))
        END DO
C
C   ***   Very slightly adjust tendencies to force exact   ***
C   ***     enthalpy, momentum and tracer conservation     ***
C
        ENTS=0.0
        UAV=0.0
        VAV=0.0
        DO 680 I=1,INB
         ENTS=ENTS+(CPN(I)*FT(I)+LV(I)*FQ(I))*(PH(I)-PH(I+1))	
         UAV=UAV+FU(I)*(PH(I)-PH(I+1))
         VAV=VAV+FV(I)*(PH(I)-PH(I+1))
  680	CONTINUE
        ENTS=ENTS/(PH(1)-PH(INB+1))
        UAV=UAV/(PH(1)-PH(INB+1))
        VAV=VAV/(PH(1)-PH(INB+1))
        DO 640 I=1,INB
         FT(I)=FT(I)-ENTS/CPN(I)
         FU(I)=(1.-CU)*(FU(I)-UAV)
         FV(I)=(1.-CU)*(FV(I)-VAV)
  640	CONTINUE
        DO 700 K=1,NTRA
         TRAAV=0.0
         DO 690 I=1,INB
          TRAAV=TRAAV+FTRA(I,K)*(PH(I)-PH(I+1))
  690    CONTINUE
         TRAAV=TRAAV/(PH(1)-PH(INB+1))
         DO 695 I=1,INB
          FTRA(I,K)=FTRA(I,K)-TRAAV
  695    CONTINUE
  700	CONTINUE
