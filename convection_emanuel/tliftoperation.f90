module tliftoperation
    private CPD,CPV,CL,RV,RD,LV0
C ---------------------------------------------------------------------------
C
        SUBROUTINE TLIFT(P,T,Q,QS,GZ,ICB,NK,TVP,TPK,CLW,ND,NL,KK)
        REAL GZ(ND),TPK(ND),CLW(ND),P(ND)
        REAL T(ND),Q(ND),QS(ND),TVP(ND),LV0
C
C   ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***
C
        CPD=1005.7
        CPV=1870.0
        CL=2500.0
        RV=461.5
        RD=287.04
        LV0=2.501E6
C
        CPVMCL=CL-CPV
        EPS=RD/RV
        EPSI=1./EPS
C
C   ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***
C
        AH0=(CPD*(1.-Q(NK))+CL*Q(NK))*T(NK)+Q(NK)*(LV0-CPVMCL*(
     1   T(NK)-273.15))+GZ(NK)
        CPP=CPD*(1.-Q(NK))+Q(NK)*CPV
        CPINV=1./CPP
C
        IF(KK.EQ.1)THEN
C
C   ***   CALCULATE LIFTED PARCEL QUANTITIES BELOW CLOUD BASE   ***
C
        DO 50 I=1,ICB-1
         CLW(I)=0.0
   50   CONTINUE
        DO 100 I=NK,ICB-1
         TPK(I)=T(NK)-(GZ(I)-GZ(NK))*CPINV
         TVP(I)=TPK(I)*(1.+Q(NK)*EPSI)
  100   CONTINUE
        END IF
C
C    ***  FIND LIFTED PARCEL QUANTITIES ABOVE CLOUD BASE    ***
C
        NST=ICB
        NSB=ICB
        IF(KK.EQ.2)THEN  
         NST=NL
         NSB=ICB+1
        END IF
        DO 300 I=NSB,NST
         TG=T(I)
         QG=QS(I)
         ALV=LV0-CPVMCL*(T(I)-273.15)
         DO 200 J=1,2
          S=CPD+ALV*ALV*QG/(RV*T(I)*T(I))
          S=1./S
          AHG=CPD*TG+(CL-CPD)*Q(NK)*T(I)+ALV*QG+GZ(I)
          TG=TG+S*(AH0-AHG)
          TG=MAX(TG,35.0)
          TC=TG-273.15
          DENOM=243.5+TC
          IF(TC.GE.0.0)THEN  
           ES=6.112*EXP(17.67*TC/DENOM)
          ELSE  
           ES=EXP(23.33086-6111.72784/TG+0.15215*LOG(TG))
          END IF  
          QG=EPS*ES/(P(I)-ES*(1.-EPS))
  200    CONTINUE
         TPK(I)=(AH0-(CL-CPD)*Q(NK)*T(I)-GZ(I)-ALV*QG)/CPD
         CLW(I)=Q(NK)-QG
         CLW(I)=MAX(0.0,CLW(I))
         RG=QG/(1.-Q(NK))
         TVP(I)=TPK(I)*(1.+RG*EPSI)
  300   CONTINUE
        RETURN
        END
