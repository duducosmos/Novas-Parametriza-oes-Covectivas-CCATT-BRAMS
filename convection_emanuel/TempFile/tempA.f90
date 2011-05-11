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
          
        enddo do_a_hundred_three
      
      !DO 110 I=ICB+1,INB
        M(ICB+1:INB)=M(ICB+1:INB)/DBOSUM  
  !110 CONTINUE     
!
!   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
!   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
!   ***                        FRACTION (SIJ)                          ***
!
        do_a_hundred_seventy: DO  I=ICB+1,INB
         
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
