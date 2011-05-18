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

!***************************************************************************************************
!>   *** FIND TEMPERATURE UP THROUGH ICB AND TEST FOR INSTABILITY           ***
!***************************************************************************************************
module INSTA
    use VarDeclared
    use LEVEL
    use TLIFT
    implicit none
    
    contains 
        subroutine copes_TInstability
            implicit none
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
end module INSTA
