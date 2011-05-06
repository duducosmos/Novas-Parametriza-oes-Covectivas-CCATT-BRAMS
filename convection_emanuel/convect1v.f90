! Convective parametrization based in K. A. Emanuel (1991,1999) scheme
!    Copyright (C) 2011  Eduardo S. Pereira
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

                  
module copes !>Convective parametrization based in K. A. Emanuel (1991,1999) scheme
    implicit none
    private    
!>-----------------------------------------------------------------------------
!>    *** On input:      ***
!>
!>     T:   Array of absolute temperature (K) of dimension ND, with first
!>           index corresponding to lowest model level. Note that this array
!>           will be altered by the subroutine if dry convective adjustment
!>           occurs and if IPBL is not equal to 0.
!>
!>     Q:   Array of specific humidity (gm/gm) of dimension ND, with first
!>            index corresponding to lowest model level. Must be defined
!>            at same grid levels as T. Note that this array will be altered
!>            if dry convective adjustment occurs and if IPBL is not equal to 0.
!>
!>     QS:  Array of saturation specific humidity of dimension ND, with first
!>            index corresponding to lowest model level. Must be defined
!>            at same grid levels as T. Note that this array will be altered
!>            if dry convective adjustment occurs and if IPBL is not equal to 0.
!>
!>     U:   Array of zonal wind velocity (m/s) of dimension ND, witth first
!>            index corresponding with the lowest model level. Defined at
!>            same levels as T. Note that this array will be altered if
!>            dry convective adjustment occurs and if IPBL is not equal to 0.
!>
!>     V:   Same as U but for meridional velocity.
!>
!>     TRA: Array of passive tracer mixing ratio, of dimensions (ND,NTRA),
!>            where NTRA is the number of different tracers. If no
!>            convective tracer transport is needed, define a dummy
!>            input array of dimension (ND,1). Tracers are defined at
!>            same vertical levels as T. Note that this array will be altered
!>            if dry convective adjustment occurs and if IPBL is not equal to 0.
!>
!>     P:   Array of pressure (mb) of dimension ND, with first
!>            index corresponding to lowest model level. Must be defined
!>            at same grid levels as T.
!>
!>     PH:  Array of pressure (mb) of dimension ND+1, with first index
!>            corresponding to lowest level. These pressures are defined at
!>            levels intermediate between those of P, T, Q and QS. The first
!>            value of PH should be greater than (i.e. at a lower level than)
!>            the first value of the array P.
!>
!>     ND:  The dimension of the arrays T,Q,QS,P,PH,FT and FQ
!>
!>     NL:  The maximum number of levels to which convection can
!>            penetrate, plus 1.
!>            NL MUST be less than or equal to ND-1.
!>
!>     NTRA:The number of different tracers. If no tracer transport
!>            is needed, set this equal to 1. (On most compilers, setting
!>            NTRA to 0 will bypass tracer calculation, saving some CPU.)  
!>
!>     DELT: The model time step (sec) between calls to CONVECT
!>

    real,allocatable ,dimension(:) ::  T,Q,QS,U,V,P,PH
    real,allocatable, dimension(:,:) :: TRA
    integer :: ND !>It will be obtained by using size function.
    integer :: NTRA !>It will be obtained by using size function.
    integer :: NL
!>
!>*****************************************************************************************************
!>

!>----------------------------------------------------------------------------
!>    ***   On Output:         ***
!>
!>     IFLAG: An output integer whose value denotes the following:
!>
!>                VALUE                        INTERPRETATION
!>                -----                        --------------
!>                  0               No moist convection; atmosphere is not
!>                                 unstable, or surface temperature is less
!>                                 than 250 K or surface specific humidity
!>                                 is non-positive.
!>
!>                 1               Moist convection occurs.
!>
!>                 2               No moist convection: lifted condensation
!>                                 level is above the 200 mb level.
!>
!>                 3               No moist convection: cloud base is higher
!>                                 then the level NL-1.
!>
!>                 4               Moist convection occurs, but a CFL condition
!>                                 on the subsidence warming is violated. This
!>                                 does not cause the scheme to terminate.

    integer :: IFLAG

!>
!>************************************************************************************************************
!>

!>
!>    FT:   Array of temperature tendency (K/s) of dimension ND, defined at same
!>            grid levels as T, Q, QS and P.

!>    FQ:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
!>            defined at same grid levels as T, Q, QS and P.
!>
!>    FU:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
!>            defined at same grid levels as T.
!>
!>    FV:   Same as FU, but for forcing of meridional velocity.

    real,allocatable, dimension(:) ::  FT,FQ,FU,FV !>Dimension ND
!>
!>************************************************************************************************************
!>
!>
!>    FTRA: Array of forcing of tracer content, in tracer mixing ratio per
!>            second, defined at same levels as T. Dimensioned (ND,NTRA).

!>
!>     PRECIP: Scalar convective precipitation rate (mm/day).
!>
!>     WD:    A convective downdraft velocity scale. For use in surface
!>             flux parameterizations. See convect.ps file for details.
!>
!>     TPRIME: A convective downdraft temperature perturbation scale (K).
!>              For use in surface flux parameterizations. See convect.ps
!>              file for details.
!>
!>    QPRIME: A convective downdraft specific humidity
!>              perturbation scale (gm/gm).
!>              For use in surface flux parameterizations. See convect.ps
!>              file for details.
!>
!>     CBMF:   The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
!>              BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
!>              ITS NEXT CALL. That is, the value of CBMF must be "remembered"
!>              by the calling program between calls to CONVECT.

!>    ***  THE PARAMETER NA SHOULD IN GENERAL BE GREATER THAN   ***
            !>    ***                OR EQUAL TO  ND + 1                    ***

    integer :: NA !> NA >= ND+1            
    !>

    real, allocatable,dimension(:,:) :: FTRA !>ND, NTRA
    real ::PRECIP,WD,TPRIME,QPRIME,CBMF


    !>Variables that are common in the subroutines
    !>DELT is the time steep that copem_convection is called by principal program
    real DELT,DELTI,RDCP
    real, allocatable,dimension(:,:) :: UENT,VENT,MENT,QENT,ELIJ,SIJ  !>NA:NA
    real, allocatable,dimension(:,:) :: TRAP  !>Na:NTRA
    real, allocatable,dimension(:,:,:) :: TRAENT !>NA:NA:NTRA


    real, allocatable,dimension(:) :: NENT,UP,VP,M,MP,TVP,TV,WATER,QP,EP,TH,WT,EVAP,CLW, &
                                  SIGP,TP,TOLD,CPN,LV,LVCP,H,HP,GZ,HM  !>All this matriz have dimension equal NA
    real,allocatable,dimension(:) :: TRATM,TRAE

    !>


!>
!>************************************************************************************************************************************************
!>

!> -----------------------------------------------------------------------
!>
!>   ***                     Specify Switches                         ***
!>
!>   ***   IPBL: Set to zero to bypass dry adiabatic adjustment       ***
!>   ***    Any other value results in dry adiabatic adjustment       ***
!>   ***     (Zero value recommended for use in models with           ***
!>   ***                   boundary layer schemes)                    ***
!>
!>   ***   MINORIG: Lowest level from which convection may originate  ***
!>   ***     (Should be first model level at which T is defined       ***
!>   ***      for models using bulk PBL schemes; otherwise, it should ***
!>   ***      be the first model level at which T is defined above    ***
!>   ***                      the surface layer)                      ***

    integer, parameter :: IPBL=0,MINORIG=1
!>
!>*********************************************************************************************************************************************
!>

!>
!>------------------------------------------------------------------------------
!>
!>   ***                    SPECIFY PARAMETERS                        ***
!>
!>   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
!>   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
!>   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
!>   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
!>   ***               BETWEEN 0 C AND TLCRIT)                        ***
!>   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
!>   ***                       FORMULATION                            ***
!>   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
!>   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
!>   ***                        OF CLOUD                              ***
!>   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
!>   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
!>   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!>   ***                          OF RAIN                             ***
!>   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
!>   ***                          OF SNOW                             ***
!>   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
!>   ***                         TRANSPORT                            ***
!>   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
!>   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
!>   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
!>   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
!>   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
!>   ***       real DELT,DELTI,RDCP,             (DAMP MUST BE LESS THAN 1)                 ***
!>

    real, parameter :: ELCRIT=0.0011,TLCRIT=-55.0,ENTP=1.5,SIGD=0.05,SIGS=0.12, &
                               OMTRAIN=50.0,OMTSNOW=5.5,COEFFR=1.0,COEFFS=0.8,CU=0.7,   &
                               BETA=10.0,DTMAX=0.9,ALPHA=0.2,DAMP=0.1


!>
!>****************************************************************************************************************************************************************************************
!>

!>   ***        ASSIGN VALUES OF THERMODYNAMIC CONSTANTS,        ***
!>   ***            GRAVITY, AND LIQUID WATER DENSITY.           ***
!>   ***             THESE SHOULD BE CONSISTENT WITH             ***
!>   ***              THOSE USED IN CALLING PROGRAM              ***
!>   ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***

    real, parameter :: CPD=1005.7,CPV=1870.0,CL=2500.0,RV=461.5,RD=287.04,LV0=2.501E6,&
                               G=9.8,ROWL=1000.0 
    real, parameter :: CPVMCL=CL-CPV, EPS=RD/RV, EPSI=1.0/EPS, GINV=1.0/G
 
    logical :: T1_EntryMatriz !> This Variable is used by copem_set_matriz_at1, if the input parameters are O.K. the code continue and this varyable is set as .TRUE. 
    
    public copes_convection,copes_TLIFT



!>
!>****************************************************************************************************************************************************************************************
!>
    contains

!>
!>*********************************************************************************************************************************************************************************
!>!>
!>!>
!>*********************************************************************************************************************************************************************************
!>
        subroutine copes_convection(T1,   Q1,    QS1,     U1,    V1,      TRA1,    P1,    PH1,&
                              NL, DELT, IFLAG,  FT,     FQ,   FU,&
                              FV,  FTRA, PRECIP, WD,   TPRIME, QPRIME, CBMF)    
            implicit none 
            real, dimension(:) :: T1,Q1,QS1,U1,V1,P1,PH1
            real, dimension(:,:) :: TRA1,FTRA
   
            real, dimension(:) :: FT,FQ,FU,FV
            
            integer :: NL,IFLAG
            real :: DELT,PRECIP,WD,TPRIME,QPRIME,CBMF


            integer :: IHMIN,NK
            real :: AHMIN, AHMAX
            integer :: i,j,k,l
           

            call copes_set_matriz_at1(T1,Q1,QS1,U1,V1,TRA1,P1,PH1) !> The first automatic test of this parametrization, difine if  T1_EntryMatriz will be .TRUE. or .FALSE.


            if( T1_EntryMatriz .eqv. .TRUE.) then
                DELTI=1.0/DELT

                !>           ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***

                FT=0.0;FQ=0.0;FU=0.0;FV=0.0;FTRA=0.0 !> Equivale a do 5 de convect43c.f

                do_seven: DO  I=1,NL+1
                              RDCP=(RD*(1.-Q(I))+Q(I)*RV)/(CPD*(1.-Q(I))+Q(I)*CPV)
                              TH(I)=T(I)*(1000.0/P(I))**RDCP !>Temperatura Potencial do ar umido
                          enddo do_seven

                PRECIP=0.0; WD=0.0; TPRIME=0.0; QPRIME=0.0; IFLAG=0

                if(IPBL /= 0 ) then
                   call copes_AdiabaticAdjustment()
                endif

                call copes_Geo_Heat_SEnergy(NL,IHMIN,AHMIN) !CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY

!
                call copes_Level(NL,IHMIN,AHMAX,NK) !Find that model level below the level of minimum moist !

!>
!>
!>
!>
!>
                call copes_free_matriz!>
            else
                print*,'An erro in the input parameters'
                continue
            end if

            


        end subroutine  copes_convection
!>
!>********************************************************************************************************************************
!>!> *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY
!>!>
!>********************************************************************************************************************************
!>

        subroutine copes_Geo_Heat_SEnergy(NL,IHMIN,AHMIN)
            implicit none
            real :: AHMIN,TVX,TVY
            integer :: NL,IHMIN
            integer :: i,j,k,l

            GZ(1)=0.0
            CPN(1)=CPD*(1.-Q(1))+Q(1)*CPV
            H(1)=T(1)*CPN(1)
            LV(1)=LV0-CPVMCL*(T(1)-273.15)
            HM(1)=LV(1)*Q(1)
            TV(1)=T(1)*(1.+Q(1)*EPSI-Q(1))
            AHMIN=1.0E12
            IHMIN=NL
            do_forty: DO I=2,NL+1
                TVX=T(I)*(1.+Q(I)*EPSI-Q(I))
                TVY=T(I-1)*(1.+Q(I-1)*EPSI-Q(I-1))
                GZ(I)=GZ(I-1)+0.5*RD*(TVX+TVY)*(P(I-1)-P(I))/PH(I)
                CPN(I)=CPD*(1.-Q(I))+CPV*Q(I)
                H(I)=T(I)*CPN(I)+GZ(I)
                LV(I)=LV0-CPVMCL*(T(I)-273.15)
                HM(I)=(CPD*(1.-Q(I))+CL*Q(I))*(T(I)-T(1))+LV(I)*Q(I)+GZ(I)
                TV(I)=T(I)*(1.+Q(I)*EPSI-Q(I))
!>
!>  ***  Find level of minimum moist static energy    ***

                IF(I.GE.MINORIG.AND.HM(I).LT.AHMIN.AND.HM(I).LT.HM(I-1))THEN
                    AHMIN=HM(I)
                    IHMIN=I
                END IF

            end do do_forty

            IHMIN=MIN(IHMIN, NL-1)

        end subroutine copes_Geo_Heat_SEnergy

!**************************************************************************************************************************************
!
!>  ***     Find that model level below the level of minimum moist       ***
!>  ***  static energy that has the maximum value of moist static energy ***
!
!**************************************************************************************************************************************

        subroutine copes_Level(NL,IHMIN,AHMAX,NK)
            integer :: IHMIN,NK,ICB,NL
            real :: AHMAX,RH,CHI,PLCL
            integer :: i,j,k,l

!
!>  ***     Find that model level below the level of minimum moist       ***
!>  ***  static energy that has the maximum value of moist static energy ***
!
            AHMAX=0.0
            do_forty_tow: DO  I=MINORIG,IHMIN
                IF(HM(I).GT.AHMAX)THEN
                    NK=I
                    AHMAX=HM(I)
                END IF
            end do do_forty_tow
!
!  ***  CHECK WHETHER PARCEL LEVEL TEMPERATURE AND SPECIFIC HUMIDITY   ***
!  ***                          ARE REASONABLE                         ***
!  ***      Skip convection if HM increases monotonically upward       ***
!
            IF(T(NK).LT.250.0.OR.Q(NK).LE.0.0.OR.IHMIN.EQ.(NL-1))THEN
                IFLAG=0
                CBMF=0.0
                RETURN
            END IF
!
!   ***  CALCULATE LIFTED CONDENSATION LEVEL OF AIR AT PARCEL ORIGIN LEVEL ***
!   ***       (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)      ***
!
            RH=Q(NK)/QS(NK)
            CHI=T(NK)/(1669.0-122.0*RH-T(NK))
            PLCL=P(NK)*(RH**CHI)
            IF(PLCL.LT.200.0.OR.PLCL.GE.2000.0)THEN
                IFLAG=2
                CBMF=0.0
                RETURN
            END IF
!
!   ***  CALCULATE FIRST LEVEL ABOVE LCL (=ICB)  ***
!
            ICB=NL-1
            DO 50 I=NK+1,NL
               IF(P(I).LT.PLCL)THEN
                   ICB=MIN(ICB,I)
               END IF
            50   CONTINUE
            IF(ICB.GE.(NL-1))THEN
               IFLAG=3
               CBMF=0.0
               RETURN
            END IF
        end subroutine copes_Level




!***************************************************************************************************
!>   *** FIND TEMPERATURE UP THROUGH ICB AND TEST FOR INSTABILITY           ***
!***************************************************************************************************

        subroutine copes_TInstability(NL,NK)
            integer :: NL,NK,ICB
            integer :: i,j,k,l

!
!   *** SUBROUTINE TLIFT CALCULATES PART OF THE LIFTED PARCEL VIRTUAL      ***
!   ***  TEMPERATURE, THE ACTUAL TEMPERATURE AND THE ADIABATIC             ***
!   ***                   LIQUID WATER CONTENT                             ***
!
            CALL copes_TLIFT(ICB,NK,TVP,TP,CLW,NL,1)
            do_fifty_four: DO  I=NK,ICB
                TVP(I)=TVP(I)-TP(I)*Q(NK)
            end do do_fifty_four
!
!   ***  If there was no convection at last time step and parcel    ***
!   ***       is stable at ICB then skip rest of calculation        ***
!
            IF(CBMF.EQ.0.0.AND.TVP(ICB).LE.(TV(ICB)-DTMAX))THEN
                IFLAG=0
                RETURN
            END IF
!
!   ***  IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY ***
!
            IF(IFLAG.NE.4)IFLAG=1
!
!   ***  FIND THE REST OF THE LIFTED PARCEL TEMPERATURES          ***
!
            CALL copes_TLIFT(ICB,NK,TVP,TP,CLW,NL,2)

            do_fifty_seven: DO  I=1,NK
                EP(I)=0.0
                SIGP(I)=SIGS
            end do do_fifty_seven

            do_sixty: DO I=NK+1,NL
                TCA=TP(I)-273.15
                IF(TCA.GE.0.0)THEN
                 ELACRIT=ELCRIT
                ELSE
                 ELACRIT=ELCRIT*(1.0-TCA/TLCRIT)
                END IF
                ELACRIT=MAX(ELACRIT,0.0)
	          EPMAX=0.999
                EP(I)=EPMAX*(1.0-ELACRIT/MAX(CLW(I),1.0E-8))
                EP(I)=MAX(EP(I),0.0)
                EP(I)=MIN(EP(I),EPMAX)
                SIGP(I)=SIGS
            end do do_sixty

!
!>   ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***
!>  ***                    VIRTUAL TEMPERATURE                    ***
!
            !do_sity_four: DO  
            I=ICB+1,NL 
            TVP(I:NL)=TVP(I:NL)-TP(I:NL)*Q(NK)! Equivale ao do 64 de convect43c.f
            !end do do_sity_four

            TVP(NL+1)=TVP(NL)-(GZ(NL+1)-GZ(NL))/CPD


        end subroutine copes_TInstability



!>
!>************************************************************************************************************************************************
!>
!> PERFORM DRY ADIABATIC ADJUSTMENT
!>
!>************************************************************************************************************************************************
!>
        subroutine copes_AdiabaticAdjustment
            implicit none
            logical :: lcomp
            integer :: JC,JN,ierr
            integer :: i,j,k,l

            real :: soma,THBAR,AHM,RM,UM,VM,DPHINV,A2,x,RDCP,TNEW,ALV,ALVNEW, &
                            QNEW,TC
!>            real, dimension(NA) :: TH,TRATM,TOLD

  

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

                           !>DO K=1,NTRA
                           !>   TRATM(K)=0.0
                           !>END DO

                               do_fifteen: DO  J=I,JN
                                            AHM=AHM+(CPD*(1.-Q(J))+Q(J)*CPV)*T(J)*(PH(J)-PH(J+1)) !>(c_{pd)(1-q^{j})+c_{pv}q^{j})T^{j}\delta p^{j}
                                            RM=RM+Q(J)*(PH(J)-PH(J+1))
                                            UM=UM+U(J)*(PH(J)-PH(J+1))
                                            VM=VM+V(J)*(PH(J)-PH(J+1))

                                            DO K=1,NTRA
                                               TRATM(K)=TRATM(K)+TRA(J,K)*(PH(J)-PH(J+1))
                                            END DO
   
                                         end do do_fifteen
   
                               DPHINV=1./(PH(I)-PH(JN+1))
                               RM=RM*DPHINV !>Media de Q no nivel J
                               UM=UM*DPHINV !>Media de U no nivel J
                               VM=VM*DPHINV !>Media de V no nivel J

                               !>DO K=1,NTRA
                               !>   TRATM(K)=TRATM(K)*DPHINV
                               !>END DO
                               TRATM(1:NTRA) = TRATM(1:NTRA)*DPHINV !>Media de TRA no nível J
 
                               A2=0.0
 
                               do_twenty: DO J=I,JN
                                  Q(J)=RM
                                  U(J)=UM                  
                                  V(J)=VM 

                                  !>DO K=1,NTRA
                                  TRA(J,:)=TRATM
                                  !>END DO

                                  RDCP=(RD*(1.-Q(J))+Q(J)*RV)/(CPD*(1.-Q(J))+Q(J)*CPV)  !>R/c_{p}
                                  X=(0.001*P(J))**RDCP !>(p^{j}/p_{0})^{R/c_{p}}
                                  TOLD(J)=T(J) !>Valor Original da Temperatura
                                  T(J)=X
                                  A2=A2+(CPD*(1.-Q(J))+Q(J)*CPV)*X*(PH(J)-PH(J+1)) !>(c_{pd}(1-q_{m})+c_{pv}q_{m})(p^{j}/p_{0})^{R/c_{p}}\delta p^{j}

                               enddo do_twenty
   
                               do_twenty_five: DO J=I,JN
                                  TH(J)=AHM/A2 !>\teta^{j}=((c_{pd)(1-q^{j})+c_{pv}q^{j})T^{j}\delta p^{j})/((c_{pd}(1-q_{m})+c_{pv}q_{m})(p^{j}/p_{0})^{R/c_{p}}\delta p^{j})
                                  T(J)=T(J)*TH(J) !>T^{j}=\teta^{i}(p^{j}/p_{0})^{R/c_{p}}
                                  TC=TOLD(J)-273.15
                                  ALV=LV0-CPVMCL*TC !>Calor Latente L_{v} = L_{v0}+(c_{pv}-c_{l})(T-273.15)
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

            !>   ***   Remove any supersaturation that results from adjustment ***

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


        end subroutine copes_AdiabaticAdjustment

!>*****************************************************************************************************************************************
!>!>The Subroutine copem_set_matriz_at1 make a test of the dimension of the arrays  T1,Q1,QS1,P1,PH1,FT1,FQ1
!>!> And allocate the equivalente arrays T,Q,QS,P,PH,FT,FQ
!>!> These array must have same dimension ND. If its OK, all matriz are dimensioned.
!>*****************************************************************************************************************************************

        subroutine copes_set_matriz_at1(T1,Q1,QS1,U1,V1,TRA1,P1,PH1)
        integer :: N1,N2,N3,N4,N5,N6,N7,N8,N9
        real, dimension(:) :: T1,Q1,QS1,U1,V1,P1,PH1
        real, dimension(:,:) :: TRA1
        integer, dimension(2) :: NDNTRA
        integer :: ierr1,ierr2,ierr3,ierr4,ierr5,ierr6,ierr7,ierr8

            N1=size(T1); N2=size(Q1); N3=size(QS1); N4=size(U1); N5=size(V1); N6=size(P1); N7=size(PH1)
            NDNTRA = shape(TRA1)
            N8=NDNTRA(1)
            N9=NDNTRA(2)

            NTRA = N9

            
            if(N1 == N2 .and. N2 == N3 .and. N3 == N4 .and. N4 == N5 .and. N5 == N6 .and. N6 == N8) then
                ND= N6
                if(N7 == N6 +1 ) then
                    NA=ND+1
                    T1_EntryMatriz = .TRUE.
                    allocate(T(ND),stat=ierr1)
                    if(ierr1 /= 0) print*, 'Allocated error T'
                    allocate(Q(ND),stat=ierr2)
                    if(ierr2 /= 0) print*, 'Allocated error Q'
                    allocate(QS(ND),stat=ierr3)
                    if(ierr3 /= 0) print*, 'Allocated error QS'
                    allocate(U(ND),stat=ierr4)
                    if(ierr4 /= 0) print*, 'Allocated error U'
                    allocate(V(ND),stat=ierr5)
                    if(ierr5 /= 0) print*, 'Allocated error V'
                    allocate(P(ND),stat=ierr6)
                    if(ierr6 /= 0) print*, 'Allocated error P'
                    allocate(PH(NA),stat=ierr7)
                    if(ierr7 /= 0) print*, 'Allocated error U'
                    allocate(TRA(ND,NTRA),stat=ierr8)
                    if(ierr8 /= 0) print*, 'Allocated error TRA'

                    allocate(TH(NA))
                    allocate(FTRA(ND,NTRA))
                    allocate(FQ(ND)); allocate(FU(ND)); allocate(FV(ND))
                    allocate(FT(ND)); allocate(UENT(NA,NA)); allocate(VENT(NA,NA)); allocate(TRAENT(NA,NA,NTRA))
                    allocate(TRATM(NA)); allocate(UP(NA)); allocate(TRAP(NA,NTRA)); allocate(M(NA)); allocate(MP(NA))
                    allocate(MENT(NA,NA)); allocate(QENT(NA,NA)); allocate(ELIJ(NA,NA)); allocate(SIJ(NA,NA))
                    allocate(TVP(NA)); allocate(TV(NA)); allocate(WATER(NA)); allocate(QP(NA)); allocate(EP(NA))
                    allocate(WT(NA)); allocate(EVAP(NA)); allocate(CLW(NA)); allocate(SIGP(NA)); allocate(TP(NA))
                    allocate(TOLD(NA)); allocate(CPN(NA)); allocate(LV(NA)); allocate(LVCP(NA)); allocate(H(NA))
                    allocate(HP(NA)); allocate(GZ(NA)); allocate(HM(NA));allocate(NENT(NA))




                    if(ierr1 /= 0 .or. ierr2 /= 0 .or. ierr3 /= 0 .or. ierr4 /= 0 &
                       .or. ierr5 /= 0 .or. ierr6 /= 0 .or. ierr7 /= 0 .or. ierr8 /= 0 ) then

                        print*,"Allocated error, verify the code"
                        T1_EntryMatriz = .FALSE.
                        call copem_free_matriz

                    endif
                   
                else
                    print*,"The Dimension of PH array must be equal",ND,&
                           "+1, but it has dimension of ",N7
                    print*,"The  copem_convection will do nothing until this correction be made."
                    print*,"Verify your datas"
                    print*,"Good look for you"
                    T1_EntryMatriz = .FALSE.

                endif



            else
                print*,"The Dimension of arrays T,Q,QS,U,V and P, and the numbers of lines of TRA, must be equals,",&
                       "but they  have, respecively, the following dimensions: ",N1,N2,N3,N4,N5,N6,N8
                print*,"The  copem_convection will do nothing until this correction be made."
                print*,"Verify your datas."
                print*,"Good look for you."
                T1_EntryMatriz = .FALSE.

            endif

            

        end subroutine copes_set_matriz_at1
!>
!>************************************************************************************************************************************
!>
!> This subroutine is used to deallocate arrays T,Q,QS,P,PH,FT,FQ and matriz TRA
!>
!>************************************************************************************************************************************
!>
        subroutine copes_free_matriz
            deallocate(T); deallocate(Q); deallocate(QS); deallocate(U)
            deallocate(V); deallocate(P); deallocate(PH); deallocate(TRA)
            deallocate(TH); deallocate(FV); deallocate(FU); deallocate(FQ)
            deallocate(FT); deallocate(FTRA)

            deallocate(TH); deallocate(FTRA); deallocate(FQ); deallocate(FU); deallocate(FV)
            deallocate(FT); deallocate(UENT); deallocate(VENT); deallocate(TRAENT)
            deallocate(TRATM); deallocate(UP); deallocate(TRAP); deallocate(M); deallocate(MP)
            deallocate(MENT); deallocate(QENT); deallocate(ELIJ); deallocate(SIJ)
            deallocate(TVP); deallocate(TV); deallocate(WATER); deallocate(QP); deallocate(EP)
            deallocate(WT); deallocate(EVAP); deallocate(CLW); deallocate(SIGP); deallocate(TP)
            deallocate(TOLD); deallocate(CPN); deallocate(LV); deallocate(LVCP); deallocate(H)
            deallocate(HP); deallocate(GZ); deallocate(HM); deallocate(NENT)
           

        end subroutine copes_free_matriz



        SUBROUTINE copes_TLIFT(ICB,NK,TVP,TPK,CLW,NL,KK)
            integer :: NL, KK,NK,NST,NSB,ICB
            integer :: i,j,k,l


            real :: CPINV, CPP,TG,QG,ALV,S,TC,DENOM,ES,AH0,AHG,RG
            real, dimension(ND) :: CLW,TPK,TVP 
!   ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***
!
            AH0=(CPD*(1.-Q(NK))+CL*Q(NK))*T(NK)+Q(NK)*(LV0-CPVMCL*(T(NK)-273.15))+GZ(NK)
            CPP=CPD*(1.-Q(NK))+Q(NK)*CPV
            CPINV=1./CPP
!
            IF(KK.EQ.1)THEN
!
!   ***   CALCULATE LIFTED PARCEL QUANTITIES BELOW CLOUD BASE   ***
!
               do_fifty: DO  I=1,ICB-1
                   CLW(I)=0.0
               end do do_fifty
               do_a_hundred: DO I=NK,ICB-1
                   TPK(I)=T(NK)-(GZ(I)-GZ(NK))*CPINV
                   TVP(I)=TPK(I)*(1.+Q(NK)*EPSI)
               end do do_a_hundred
            END IF

!    ***  FIND LIFTED PARCEL QUANTITIES ABOVE CLOUD BASE    ***
!
            NST=ICB
            NSB=ICB
            IF(KK.EQ.2)THEN  
                NST=NL
                NSB=ICB+1
            END IF
            do_three_hunded: DO I=NSB,NST
                TG=T(I)
                QG=QS(I)
                ALV=LV0-CPVMCL*(T(I)-273.15)
                do_tow_hundred: DO J=1,2
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
                end do do_tow_hundred
                TPK(I)=(AH0-(CL-CPD)*Q(NK)*T(I)-GZ(I)-ALV*QG)/CPD
                CLW(I)=Q(NK)-QG
                CLW(I)=MAX(0.0,CLW(I))
                RG=QG/(1.-Q(NK))
                TVP(I)=TPK(I)*(1.+RG*EPSI)
            end do do_three_hunded
            RETURN
        end subroutine copes_TLIFT




end module  copes


