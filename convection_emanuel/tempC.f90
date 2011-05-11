
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
C
C   ***           RETURN           ***
C
