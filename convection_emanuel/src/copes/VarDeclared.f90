!> Convective parametrization based in K. A. Emanuel (1991,1999) scheme
!>   Copyright (C) 2011  Grupo de Modelagem da Atmosfera e Interfaces (GMAI)
!>                http://meioambiente.cptec.inpe.br/gmai/index.php
!>
!> @Author Eduardo S. Pereira
!>
!>    This program is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License as published by
!>    the Free Software Foundation, either version 3 of the License, or
!>    (at your option) any later version.
!>
!>    This program is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> This module contain the commom block constants and variables used by copes module

module VarDeclared
    implicit none
!>-----------------------------------------------------------------------------
!>    *** On input:      ***
!>
!>     T:   Array of absolute temperature (K) of dimension ND, with first
!>           index corresponding to lowest model level. Note that this array
!>           will be altered by the subroutine if dry convective adjustment
!>           occurs and if IPBL is not equal to 0.
!>copes
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
!>******************************************************************************
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
!>*******************************************************************************
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
!>********************************************************************************
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
!>*********************************************************************************************
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
!>*******************************************************************************************
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
!>   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***Constants
!>   ***       real DELT,DELTI,RDCP,             (DAMP MUST BE LESS THAN 1)                 ***
!>

    real, parameter :: ELCRIT=0.0011,TLCRIT=-55.0,ENTP=1.5,SIGD=0.05,SIGS=0.12, &
                               OMTRAIN=50.0,OMTSNOW=5.5,COEFFR=1.0,COEFFS=0.8,CU=0.7,   &
                               BETA=10.0,DTMAX=0.9,ALPHA=0.2,DAMP=0.1


!>
!>******************************************************************************************
!>

!>   ***        ASSIGN VALUES OF THERMODYNAMIC CONSTANTS,        ***
!>   ***            GRAVITY, AND LIQUID WATER DENSITY.           ***
!>   ***             THESE SHOULD BE CONSISTENT WITH             ***
!>   ***              THOSE USED IN CALLING PROGRAM              ***
!>   ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***

    real, parameter :: CPD=1005.7,CPV=1870.0,CL=2500.0,RV=461.5,RD=287.04,LV0=2.501E6,&
                               G=9.8,ROWL=1000.0 
    real, parameter :: CPVMCL=CL-CPV, EPS=RD/RV, EPSI=1.0/EPS, GINV=1.0/G


end module VarDeclared
