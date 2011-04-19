                  
module convect_emanuel
    implicit none

! -----------------------------------------------------------------------
!
!   ***                     Specify Switches                         ***
!
!   ***   IPBL: Set to zero to bypass dry adiabatic adjustment       ***
!   ***    Any other value results in dry adiabatic adjustment       ***
!   ***     (Zero value recommended for use in models with           ***
!   ***                   boundary layer schemes)                    ***
!
!   ***   MINORIG: Lowest level from which convection may originate  ***
!   ***     (Should be first model level at which T is defined       ***
!   ***      for models using bulk PBL schemes; otherwise, it should ***
!   ***      be the first model level at which T is defined above    ***
!   ***                      the surface layer)                      ***

    integer, parameter :: IPBL=0,MINORIG=1
!
!------------------------------------------------------------------------------
!
!   ***                    SPECIFY PARAMETERS                        ***
!
!   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
!   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
!   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
!   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
!   ***               BETWEEN 0 C AND TLCRIT)                        ***
!   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
!   ***                       FORMULATION                            ***
!   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
!   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
!   ***                        OF CLOUD                              ***
!   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
!   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
!   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!   ***                          OF RAIN                             ***
!   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!   ***                          OF SNOW                             ***
!   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
!   ***                         TRANSPORT                            ***
!   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
!   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
!   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
!   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
!   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
!   ***                   (DAMP MUST BE LESS THAN 1)                 ***
!

    real(kind=8), parameter :: ELCRIT=0.0011,TLCRIT=-55.0,ENTP=1.5,SIGD=0.05,SIGS=0.12, &
                               OMTRAIN=50.0,OMTSNOW=5.5,COEFFR=1.0,COEFFS=0.8,CU=0.7,   &
                               BETA=10.0,DTMAX=0.9,ALPHA=0.2,DAMP=0.1

!   ***        ASSIGN VALUES OF THERMODYNAMIC CONSTANTS,        ***
!   ***            GRAVITY, AND LIQUID WATER DENSITY.           ***
!   ***             THESE SHOULD BE CONSISTENT WITH             ***
!   ***              THOSE USED IN CALLING PROGRAM              ***
!   ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***

    real(kind=8), parameter :: CPD=1005.7,CPV=1870.0,CL=2500.0,RV=461.5,RD=287.04,LV0=2.501E6,&
                               G=9.8,ROWL=1000.0 
    real(kind=8), parameter :: CPVMCL=CL-CPV, EPS=RD/RV, EPSI=1.0/EPS, GINV=1.0/G 
    private AdiabaticAdjustment

    contains

!-----------------------------------------------------------------------------
!    *** On input:      ***
!
!     T:   Array of absolute temperature (K) of dimension ND, with first
!           index corresponding to lowest model level. Note that this array
!           will be altered by the subroutine if dry convective adjustment
!           occurs and if IPBL is not equal to 0.
!
!     Q:   Array of specific humidity (gm/gm) of dimension ND, with first
!            index corresponding to lowest model level. Must be defined
!            at same grid levels as T. Note that this array will be altered
!            if dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     QS:  Array of saturation specific humidity of dimension ND, with first
!            index corresponding to lowest model level. Must be defined
!            at same grid levels as T. Note that this array will be altered
!            if dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     U:   Array of zonal wind velocity (m/s) of dimension ND, witth first
!            index corresponding with the lowest model level. Defined at
!            same levels as T. Note that this array will be altered if
!            dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     V:   Same as U but for meridional velocity.
!
!     TRA: Array of passive tracer mixing ratio, of dimensions (ND,NTRA),
!            where NTRA is the number of different tracers. If no
!            convective tracer transport is needed, define a dummy
!            input array of dimension (ND,1). Tracers are defined at
!            same vertical levels as T. Note that this array will be altered
!            if dry convective adjustment occurs and if IPBL is not equal to 0.
!
!     P:   Array of pressure (mb) of dimension ND, with first
!            index corresponding to lowest model level. Must be defined
!            at same grid levels as T.
!
!     PH:  Array of pressure (mb) of dimension ND+1, with first index
!            corresponding to lowest level. These pressures are defined at
!            levels intermediate between those of P, T, Q and QS. The first
!            value of PH should be greater than (i.e. at a lower level than)
!            the first value of the array P.
!
!     ND:  The dimension of the arrays T,Q,QS,P,PH,FT and FQ
!
!     NL:  The maximum number of levels to which convection can
!            penetrate, plus 1.
!            NL MUST be less than or equal to ND-1.
!
!     NTRA:The number of different tracers. If no tracer transport
!            is needed, set this equal to 1. (On most compilers, setting
!            NTRA to 0 will bypass tracer calculation, saving some CPU.)  
!
!     DELT: The model time step (sec) between calls to CONVECT
!
!----------------------------------------------------------------------------
!    ***   On Output:         ***
!
!     IFLAG: An output integer whose value denotes the following:
!
!                VALUE                        INTERPRETATION
!                -----                        --------------
!                  0               No moist convection; atmosphere is not
!                                 unstable, or surface temperature is less
!                                 than 250 K or surface specific humidity
!                                 is non-positive.
!
!                 1               Moist convection occurs.
!
!                 2               No moist convection: lifted condensation
!                                 level is above the 200 mb level.
!
!                 3               No moist convection: cloud base is higher
!                                 then the level NL-1.
!
!                 4               Moist convection occurs, but a CFL condition
!                                 on the subsidence warming is violated. This
!                                 does not cause the scheme to terminate.
!
!    FT:   Array of temperature tendency (K/s) of dimension ND, defined at same
!            grid levels as T, Q, QS and P.

!    FQ:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
!            defined at same grid levels as T, Q, QS and P.
!
!    FU:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
!            defined at same grid levels as T.
!
!    FV:   Same as FU, but for forcing of meridional velocity.
!
!    FTRA: Array of forcing of tracer content, in tracer mixing ratio per
!            second, defined at same levels as T. Dimensioned (ND,NTRA).
!
!     PRECIP: Scalar convective precipitation rate (mm/day).
!
!     WD:    A convective downdraft velocity scale. For use in surface
!             flux parameterizations. See convect.ps file for details.
!
!     TPRIME: A convective downdraft temperature perturbation scale (K).
!              For use in surface flux parameterizations. See convect.ps
!              file for details.
!
!    QPRIME: A convective downdraft specific humidity
!              perturbation scale (gm/gm).
!              For use in surface flux parameterizations. See convect.ps
!              file for details.
!
!     CBMF:   The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
!              BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
!              ITS NEXT CALL. That is, the value of CBMF must be "remembered"
!              by the calling program between calls to CONVECT.


        subroutine convection(T,   Q,    QS,     U,    V,      TRA,    P,    PH,&
                              ND,  NL,   NTRA,   DELT, IFLAG,  FT,     FQ,   FU,&
                              FV,  FTRA, PRECIP, WD,   TPRIME, QPRIME, CBMF)

            

            
            

            !    ***  THE PARAMETER NA SHOULD IN GENERAL BE GREATER THAN   ***
            !    ***                OR EQUAL TO  ND + 1                    ***

            integer, parameter :: NA=70
            
            !
            real(kind=8), dimension(NA) :: NENT,UP,VP,M,MP,TVP,TV,WATER,QP,EP,TH,WT,EVAP,CLW, &
                                      SIGP,TP,TOLD,CPN,LV,LVCP,H,HP,GZ,HM,TRATM
            !
            integer :: ND,NL,NTRA,IFLAG
            integer :: i,j,k,l
            real(kind=8), dimension(ND) ::  T,Q,QS,U,V,P,PH,FT,FQ,FU,FV
            real(kind=8), dimension(ND,NTRA) :: TRA,FTRA
            real(kind=8), dimension(NA,NA) :: UENT,VENT,MENT,QENT,ELIJ,SIJ
            real(kind=8), dimension(NA,NTRA) :: TRAP
            real(kind=8), dimension(NA,NA,NTRA) :: TRAENT
            real(kind=8) ::DELT,PRECIP,WD,TPRIME,QPRIME,DELTI,RDCP,CBMF
            

            DELTI=1.0/DELT


           




            !           ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***

            FT=0.0;FQ=0.0;FU=0.0;FV=0.0;FTRA=0.0 ! Equivale a do 5 de convect43c.f

            do_seven: DO  I=1,NL+1
                          RDCP=(RD*(1.-Q(I))+Q(I)*RV)/(CPD*(1.-Q(I))+Q(I)*CPV)
                          TH(I)=T(I)*(1000.0/P(I))**RDCP !Temperatura Potencial do ar umido
                      enddo do_seven

            PRECIP=0.0; WD=0.0; TPRIME=0.0; QPRIME=0.0; IFLAG=0

            if(IPBL /= 0 ) then
               call AdiabaticAdjustment(NL,NA,ND,NTRA,TH,Q,U,V,PH,TRA,TRATM,DELT)
            endif

        end subroutine convection
!
!!
!! PERFORM DRY ADIABATIC ADJUSTMENT
!
        subroutine AdiabaticAdjustment(NL,NA,ND,NTRA,TH,Q,U,V,PH,TRA,TRATM,DELT)
            logical :: lcomp
            integer :: JC,JN,NL,NTRA,ND,NA
            integer :: I,J,K
            real(kind=8) :: soma,THBAR,AHM,RM,UM,VM,DPHINV,A2,x,RDCP,TNEW,ALV,ALVNEW,DELT, &
                            PRECIP,QNEW,TC
            real(kind=8), dimension(NTRA) :: TRANTM
            real(kind=8), dimension(ND) ::  T,Q,QS,U,V,P,PH,FT,FQ,FU,FV
            real(kind=8), dimension(ND,NTRA) :: TRA,FTRA
            real(kind=8), dimension(NA) :: NENT,UP,VP,M,MP,TVP,TV,WATER,QP,EP,TH,WT,EVAP,CLW, &
                                           SIGP,TP,TOLD,CPN,LV,LVCP,H,HP,GZ,HM,TRATM
  
 

            lcomp=.true.




            JC=0
            do_thrity: DO  I=NL-1,1,-1

                           JN=0
                           soma=TH(I)*(1.+Q(I)*EPSI-Q(I))

                           do_ten: DO J=I+1,NL
                                      soma=soma+TH(J)*(1.+Q(J)*EPSI-Q(J))
                                      THBAR=soma/FLOAT(J+1-I)
                                      IF((TH(J)*(1.+Q(J)*EPSI-Q(J))).LT.THBAR)JN=J
                                   enddo do_ten    
                     
                           IF(I.EQ.1)JN=MAX(JN,2)
                           IF(JN.EQ.0) cycle do_thrity

                           do_twelve: do while(lcomp .eqv. .true.)

                               AHM=0.0; RM=0.0; UM=0.0; VM=0.0

                               TRATM=0.0

                           !DO K=1,NTRA
                           !   TRATM(K)=0.0
                           !END DO

                               do_fifteen: DO  J=I,JN
                                            AHM=AHM+(CPD*(1.-Q(J))+Q(J)*CPV)*T(J)*(PH(J)-PH(J+1)) !(c_{pd)(1-q^{j})+c_{pv}q^{j})T^{j}\delta p^{j}
                                            RM=RM+Q(J)*(PH(J)-PH(J+1))
                                            UM=UM+U(J)*(PH(J)-PH(J+1))
                                            VM=VM+V(J)*(PH(J)-PH(J+1))

                                            DO K=1,NTRA
                                               TRATM(K)=TRATM(K)+TRA(J,K)*(PH(J)-PH(J+1))
                                            END DO
   
                                         end do do_fifteen
   
                               DPHINV=1./(PH(I)-PH(JN+1))
                               RM=RM*DPHINV !Media de Q no nivel J
                               UM=UM*DPHINV !Media de U no nivel J
                               VM=VM*DPHINV !Media de V no nivel J

                               !DO K=1,NTRA
                               !   TRATM(K)=TRATM(K)*DPHINV
                               !END DO
                               TRATM =TRATM*DPHINV !Media de TRA no nível J
 
                               A2=0.0
 
                               do_twenty: DO J=I,JN
                                  Q(J)=RM
                                  U(J)=UM                  
                                  V(J)=VM 

                                  !DO K=1,NTRA
                                  TRA(J,:)=TRATM
                                  !END DO

                                  RDCP=(RD*(1.-Q(J))+Q(J)*RV)/(CPD*(1.-Q(J))+Q(J)*CPV)  !R/c_{p}
                                  X=(0.001*P(J))**RDCP !(p^{j}/p_{0})^{R/c_{p}}
                                  TOLD(J)=T(J) !Valor Original da Temperatura
                                  T(J)=X
                                  A2=A2+(CPD*(1.-Q(J))+Q(J)*CPV)*X*(PH(J)-PH(J+1)) !(c_{pd}(1-q_{m})+c_{pv}q_{m})(p^{j}/p_{0})^{R/c_{p}}\delta p^{j}

                               enddo do_twenty
   
                               do_twenty_five: DO J=I,JN
                                  TH(J)=AHM/A2 !\teta^{j}=((c_{pd)(1-q^{j})+c_{pv}q^{j})T^{j}\delta p^{j})/((c_{pd}(1-q_{m})+c_{pv}q_{m})(p^{j}/p_{0})^{R/c_{p}}\delta p^{j})
                                  T(J)=T(J)*TH(J) !T^{j}=\teta^{i}(p^{j}/p_{0})^{R/c_{p}}
                                  TC=TOLD(J)-273.15
                                  ALV=LV0-CPVMCL*TC !Calor Latente L_{v} = L_{v0}+(c_{pv}-c_{l})(T-273.15)
                                  QS(J)=QS(J)+QS(J)*(1.+QS(J)*(EPSI-1.))*ALV*(T(J)- &
                                        TOLD(J))/(RV*TOLD(J)*TOLD(J))
                               enddo do_twenty_five

                       IF((TH(JN+1)*(1.+Q(JN+1)*EPSI-Q(JN+1))) < (TH(JN)*(1.+Q(JN)*EPSI-Q(JN))))THEN
                           JN=JN+1
                       else
                           lcomp=.false.
                       END IF
                enddo do_twelve


                IF(I.EQ.1)JC=JN 

            enddo do_thrity

            !   ***   Remove any supersaturation that results from adjustment ***

            IF(JC.GT.1)THEN
                do_thirty_eight: DO J=1,JC
                    IF(QS(J).LT.Q(J))THEN 
                        ALV=LV0-CPVMCL*(T(J)-273.15)  
                        TNEW=T(J)+ALV*(Q(J)-QS(J))/(CPD*(1.-Q(J))+CL*Q(J)+QS(J)*(CPV-CL+ALV*ALV/(RV*T(J)*T(J))))
                        ALVNEW=LV0-CPVMCL*(TNEW-273.15)
                        QNEW=(ALV*Q(J)-(TNEW-T(J))*(CPD*(1.-Q(J))+CL*Q(J)))/ALVNEW
                        PRECIP=PRECIP+24.*3600.*1.0E5*(PH(J)-PH(J+1))*(Q(J)-QNEW)/(G*DELT*ROWL)
                        T(J)=TNEW
                        Q(J)=QNEW
                        QS(J)=QNEW
                    END IF     
                enddo do_thirty_eight  
            END IF

        end subroutine AdiabaticAdjustment
end module  convect_emanuel


