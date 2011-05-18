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

module TENDENCI2
    use VarDeclared
    implicit none
    
    contains
        !
        !   ***  CALCULATE TENDENCIES OF POTENTIAL TEMPERATURE AND MIXING RATIO  ***
        !   ***               AT LEVELS ABOVE THE LOWEST LEVEL                   ***
        !
        !   ***  FIRST FIND THE NET SATURATED UPDRAFT AND DOWNDRAFT MASS FLUXES  ***
        !   ***                      THROUGH EACH LEVEL                          ***
        
        subroutine copes_tendenci2()
            integer :: I,J,K
            real :: DPINV,CPINV,AMP1,AD,AWAT

            do_five_hundred: DO I=2,INB
            
                DPINV=0.01/(PH(I)-PH(I+1))
                CPINV=1.0/CPN(I)
                AMP1=0.0
                AD=0.0
                
                IF(I.GE.NK)THEN
                
                    DO  K=I+1,INB+1
                    AMP1=AMP1+M(K)
                    end do 
                    
                END IF
                
                do_four_hundred_fifty: DO  K=1,I
                
                    DO J=I+1,INB+1
                        AMP1=AMP1+MENT(K,J)
                    ENDDO
                    
                ENDDO do_four_hundred_fifty
                
                
                IF((2.*G*DPINV*AMP1).GE.DELTI)IFLAG=4
                
                DO  K=1,I-1
                
                    DO J=I,INB
                        AD=AD+MENT(J,K)
                    ENDDO 
                    
                ENDDO 
                    
                FT(I)=FT(I)+G*DPINV*(AMP1*(T(I+1)-T(I)+(GZ(I+1)-GZ(I))* &
                        CPINV)-AD*(T(I)-T(I-1)+(GZ(I)-GZ(I-1))*CPINV))  &
                        -SIGD*LVCP(I)*EVAP(I)
                            
                FT(I)=FT(I)+G*DPINV*MENT(I,I)*(HP(I)-H(I)+&
                        T(I)*(CPV-CPD)*(Q(I)-QENT(I,I)))*CPINV
                          
                FT(I)=FT(I)+SIGD*WT(I+1)*(CL-CPD)*WATER(I+1)* &
                        (T(I+1)-T(I))*DPINV*CPINV
                          
                FQ(I)=FQ(I)+G*DPINV*(AMP1*(Q(I+1)-Q(I))- &
                        AD*(Q(I)-Q(I-1)))
                          
                FU(I)=FU(I)+G*DPINV*(AMP1*(U(I+1)-U(I))- &
                        AD*(U(I)-U(I-1)))
                            
                FV(I)=FV(I)+G*DPINV*(AMP1*(V(I+1)-V(I))- &
                        AD*(V(I)-V(I-1)))
                          
                DO K=1,NTRA
                    FTRA(I,K)=FTRA(I,K)+G*DPINV*(AMP1*(TRA(I+1,K)- &
                                TRA(I,K))-AD*(TRA(I,K)-TRA(I-1,K)))
                END DO
                    
                do_four_hundred_eighty: DO  K=1,I-1
                    
                    AWAT=ELIJ(K,I)-(1.-EP(I))*CLW(I)
                    AWAT=MAX(AWAT,0.0)
                    FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-AWAT-Q(I))
                    FU(I)=FU(I)+G*DPINV*MENT(K,I)*(UENT(K,I)-U(I))
                    FV(I)=FV(I)+G*DPINV*MENT(K,I)*(VENT(K,I)-V(I))
                        
                    DO J=1,NTRA
                        FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)- &
                                    TRA(I,J))
                    END DO
                        
                END DO do_four_hundred_eighty
                    
                do_four_hundred_ninety: DO  K=I,INB
                    
                    FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-Q(I))
                    FU(I)=FU(I)+G*DPINV*MENT(K,I)*(UENT(K,I)-U(I))
                    FV(I)=FV(I)+G*DPINV*MENT(K,I)*(VENT(K,I)-V(I))
                        
                    DO J=1,NTRA
                        
                        FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)-&
                                    TRA(I,J))
                    END DO
                        
                END DO do_four_hundred_ninety
                    
                FQ(I)=FQ(I)+SIGD*EVAP(I)+G*(MP(I+1)* &
                          (QP(I+1)-Q(I))-MP(I)*(QP(I)-Q(I-1)))*DPINV
                          
                FU(I)=FU(I)+G*(MP(I+1)*(UP(I+1)-U(I))-MP(I)* &
                        (UP(I)-U(I-1)))*DPINV
                          
                FV(I)=FV(I)+G*(MP(I+1)*(VP(I+1)-V(I))-MP(I)* &
                        (VP(I)-V(I-1)))*DPINV
                          
                DO J=1,NTRA
                    FTRA(I,J)=FTRA(I,J)+G*DPINV*(MP(I+1)*(TRAP(I+1,J)-TRA(I,J))- &
                                MP(I)*(TRAP(I,J)-TRA(I-1,J)))
                END DO
            END DO do_five_hundred
        
        end subroutine copes_tendenci2
        
end module TENDENCI2

