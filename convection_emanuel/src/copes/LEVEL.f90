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

!**************************************************************************************************************************************
!
!>  ***     Find that model level below the level of minimum moist       ***
!>  ***  static energy that has the maximum value of moist static energy ***
!
!**************************************************************************************************************************************

module LEVEL

    use VarDeclared
    implicit none
    integer :: IHMIN,NK
    real :: AHMAX

    contains 

        subroutine copes_Level
            implicit none
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
            IF(T(NK).LT.250.0 .OR. Q(NK).LE.0.0 .OR. IHMIN.EQ.(NL-1))THEN
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
            
            IF(PLCL.LT.200.0 .OR. PLCL.GE.2000.0)THEN
                IFLAG=2
                CBMF=0.0
                RETURN
            END IF
!
!   ***  CALCULATE FIRST LEVEL ABOVE LCL (=ICB)  ***
!
            ICB=NL-1
            
            DO  I=NK+1,NL
               IF(P(I).LT.PLCL)THEN
                   ICB=MIN(ICB,I)
               END IF
            END DO
            
            IF(ICB.GE.(NL-1))THEN
               IFLAG=3
               CBMF=0.0
               RETURN
            END IF
            
        end subroutine copes_Level


end module LEVEL


