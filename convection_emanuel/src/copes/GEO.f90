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
!>

!>********************************************************************************************************************************
!> *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY
!>
!>********************************************************************************************************************************
!>

module GEO
    use VarDeclared
    implicit none
    real :: AHMIN
    integer :: IHMIN
    
    contains

       subroutine copes_Geo_Heat_SEnergy
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
        
end module GEO
