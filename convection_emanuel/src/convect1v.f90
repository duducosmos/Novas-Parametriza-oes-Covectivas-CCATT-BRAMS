!> Convective parametrization based in K. A. Emanuel (1991,1999) scheme
!>   Copyright (C) 2011  Eduardo S. Pereira
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

                  
module copes !>Convective parametrization based in K. A. Emanuel (1991,1999) scheme
    
     
    use VarDeclared
    use CopesStartEnd
    use TLIFT
    implicit none
    private   
 
    
    public copes_convection,copes_TLIFT


    contains

!
!*********************************************************************************************************************************************************************************
!*********************************************************************************************************************************************************************************
!*********************************************************************************************************************************************************************************
!

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
!> *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY
!>
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
            real :: ELACRIT,EPMAX,TCA

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
            I=ICB+1!,NL 
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
        
end module  copes


