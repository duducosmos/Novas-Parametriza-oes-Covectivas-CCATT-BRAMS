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

module NENTRE
    use VarDeclared
    implicit none
    contains
        !
        !   ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***
        !   ***              PROBABILITIES OF MIXING                     ***
        !
        subroutine copes_normentr_pmix()
        
            do_tow_hundred: DO I=ICB+1,INB
            
                if_tow_hundred: IF(NENT(I).NE.0)THEN
                
                    QP1=Q(NK)-EP(I)*CLW(I)
                    ANUM=H(I)-HP(I)-LV(I)*(QP1-QS(I))
                    DENOM=H(I)-HP(I)+LV(I)*(Q(I)-QP1)
                    IF(ABS(DENOM).LT.0.01)DENOM=0.01
                    SCRIT=ANUM/DENOM
                    ALT=QP1-QS(I)+SCRIT*(Q(I)-QP1)
                    IF(ALT.LT.0.0)SCRIT=1.0
                    SCRIT=MAX(SCRIT,0.0)
                    ASIJ=0.0
                    SMIN=1.0
                    
                    do_seventy_five: DO J=ICB,INB
                    
                        if_seventy_five: IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
                        
                            IF(J.GT.I)THEN
                            
                                SMID=MIN(SIJ(I,J),SCRIT)
                                SJMAX=SMID
                                SJMIN=SMID
                                
                                IF(SMID.LT.SMIN.AND.SIJ(I,J+1).LT.SMID)THEN
                                
                                    SMIN=SMID
                                    SJMAX=MIN(SIJ(I,J+1),SIJ(I,J),SCRIT)
                                    SJMIN=MAX(SIJ(I,J-1),SIJ(I,J))
                                    SJMIN=MIN(SJMIN,SCRIT)
                                    
                                END IF
                                
                            ELSE
                            
                                SJMAX=MAX(SIJ(I,J+1),SCRIT)
                                SMID=MAX(SIJ(I,J),SCRIT)
                                SJMIN=0.0
                                IF(J.GT.1)SJMIN=SIJ(I,J-1)
                                SJMIN=MAX(SJMIN,SCRIT)
                                
                            END IF
                            
                            DELP=ABS(SJMAX-SMID)
                            DELM=ABS(SJMIN-SMID)
                            ASIJ=ASIJ+(DELP+DELM)*(PH(J)-PH(J+1))
                            MENT(I,J)=MENT(I,J)*(DELP+DELM)*(PH(J)-PH(J+1))
                            
                        END IF if_seventy_five
                        
                    end do do_seventy_five
                    
                    ASIJ=MAX(1.0E-21,ASIJ)
                    ASIJ=1.0/ASIJ
                    
                    !DO 180 J=ICB,INB
                    MENT(I,ICB:INB)=MENT(I,ICB:INB)*ASIJ
                    !180    CONTINUE
                    
                    BSUM=0.0
                    
                    DO  J=ICB,INB !do 190
                        BSUM=BSUM+MENT(I,J)
                    end do
                    
                    IF(BSUM.LT.1.0E-18)THEN
                    
                        NENT(I)=0
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
                    
                END IF if_tow_hundred
                
            end do do_tow_hundred
        end subroutine copes_normentr_pmix
end module NENTRE
