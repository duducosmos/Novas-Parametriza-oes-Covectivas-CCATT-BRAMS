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




module TLIFT

    use VarDeclared
    implicit none

    contains
    
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
        
end module TLIFT
