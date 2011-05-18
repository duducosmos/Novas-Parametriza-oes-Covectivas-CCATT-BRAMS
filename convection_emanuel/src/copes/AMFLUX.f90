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

module AMFLUX

    use VarDeclared
    implicit none

    contains
    
    subroutine copes_flux()
        integer :: I,J
        real :: QTI,BF2,ANUM,DENOM,DEI,ALTEM,CWAT,STEMP
        !
        !   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
        !   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
        !   ***                        FRACTION (SIJ)                          ***
        !
            do_a_hundred_seventy: DO  I=ICB+1,INB
         
                QTI=Q(NK)-EP(I)*CLW(I)
             
                do_a_hundred_sixty: DO  J=ICB,INB
             
                    BF2=1.+LV(J)*LV(J)*QS(J)/(RV*T(J)*T(J)*CPD)
                    ANUM=H(J)-HP(I)+(CPV-CPD)*T(J)*(QTI-Q(J))
                    DENOM=H(I)-HP(I)+(CPD-CPV)*(Q(I)-QTI)*T(J)
                    DEI=DENOM
                    IF(ABS(DEI).LT.0.01)DEI=0.01
                    SIJ(I,J)=ANUM/DEI
                    SIJ(I,I)=1.0
                    ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)
                    ALTEM=ALTEM/BF2
                    CWAT=CLW(J)*(1.-EP(J))
                    STEMP=SIJ(I,J)
                     
                    IF((STEMP.LT.0.0 .OR. STEMP .GT. 1.0 .OR. ALTEM.GT.CWAT).AND.J.GT.I)THEN
                                      
                       ANUM=ANUM-LV(J)*(QTI-QS(J)-CWAT*BF2)
                       DENOM=DENOM+LV(J)*(Q(I)-QTI)                    
                       IF(ABS(DENOM).LT.0.01)DENOM=0.01
                       SIJ(I,J)=ANUM/DENOM
                       ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)
                       ALTEM=ALTEM-(BF2-1.)*CWAT
                        
                    END IF
                        
                    IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
                        
                        QENT(I,J)=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI
                        UENT(I,J)=SIJ(I,J)*U(I)+(1.-SIJ(I,J))*U(NK)
                        VENT(I,J)=SIJ(I,J)*V(I)+(1.-SIJ(I,J))*V(NK)
                            
                         !DO K=1,NTRA
                        TRAENT(I,J,1:NTRA)=SIJ(I,J)*TRA(I,1:NTRA)+(1.-SIJ(I,J))*TRA(NK,1:NTRA)
                         !END DO
                         
                        ELIJ(I,J)=ALTEM
                        ELIJ(I,J)=MAX(0.0,ELIJ(I,J))
                        MENT(I,J)=M(I)/(1.-SIJ(I,J))
                        NENT(I)=NENT(I)+1
                            
                    END IF
                         
                    SIJ(I,J)=MAX(0.0,SIJ(I,J))
                    SIJ(I,J)=MIN(1.0,SIJ(I,J))
                         
                end do do_a_hundred_sixty
    !
    !   ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS  ***
    !   ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES  ***
    !
                IF(NENT(I).EQ.0)THEN
                 
                    MENT(I,I)=M(I)
                    QENT(I,I)=Q(NK)-EP(I)*CLW(I)
                    UENT(I,I)=U(NK)
                    VENT(I,I)=V(NK)
              !DO J=1,NTRA
                    TRAENT(I,I,1:NTRA)=TRA(NK,1:NTRA)
              !END DO
                    ELIJ(I,I)=CLW(I)
                    SIJ(I,I)=1.0
                    
                END IF 
                 
             end do do_a_hundred_seventy
             
             SIJ(INB,INB)=1.0
             
         end subroutine copes_flux
         
end module AMFLUX
