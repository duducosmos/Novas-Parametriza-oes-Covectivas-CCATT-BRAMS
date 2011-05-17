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

module PRECIDOWN
    use VarDeclared
    implicit none
    
    contains
    
        subroutine copes_precip()
            !
            !   ***  CHECK WHETHER EP(INB)=0, IF SO, SKIP PRECIPITATING    ***
            !   ***             DOWNDRAFT CALCULATION                      ***
            !
            if_four_hundred_five: IF(EP(INB).LT.0.0001)then
            
                continue
                
            ELSE
                !
                !   ***  INTEGRATE LIQUID WATER EQUATION TO FIND CONDENSED WATER   ***
                !   ***                AND CONDENSED WATER FLUX                    ***
                !
                JTT=2
                !
                !    ***                    BEGIN DOWNDRAFT LOOP                    ***
                !
                do_four_hundred: DO I=INB,1,-1
            
                    !
                    !    ***              CALCULATE DETRAINED PRECIPITATION             ***
                
                    WDTRAIN=G*EP(I)*M(I)*CLW(I)
                    
                    IF(I.GT.1)THEN
                
                        do_three_hundred_twenty:DO 320 J=1,I-1
                            AWAT=ELIJ(J,I)-(1.-EP(I))*CLW(I)
                            AWAT=MAX(0.0,AWAT)
                            WDTRAIN=WDTRAIN+G*AWAT*MENT(J,I)
                        end do do_three_hundred_twenty
                        
                    END IF
                    !
                    !    ***    FIND RAIN WATER AND EVAPORATION USING PROVISIONAL   ***
                    !    ***              ESTIMATES OF QP(I)AND QP(I-1)             ***     
                    !
                    !  ***  Value of terminal velocity and coefficient of evaporation for snow   ***
                    ! 
                    COEFF=COEFFS
                    WT(I)=OMTSNOW
                    !      
                    !  ***  Value of terminal velocity and coefficient of evaporation for rain   ***
                    !
                    IF(T(I).GT.273.0)THEN
                        COEFF=COEFFR
                        WT(I)=OMTRAIN
                    END IF
                
                    QSM=0.5*(Q(I)+QP(I+1))
                    AFAC=COEFF*PH(I)*(QS(I)-QSM)/(1.0E4+2.0E3*PH(I)*QS(I))
                    AFAC=MAX(AFAC,0.0)
                    SIGT=SIGP(I)
                    SIGT=MAX(0.0,SIGT)
                    SIGT=MIN(1.0,SIGT)
                    B6=100.*(PH(I)-PH(I+1))*SIGT*AFAC/WT(I)
                    C6=(WATER(I+1)*WT(I+1)+WDTRAIN/SIGD)/WT(I)
                    REVAP=0.5*(-B6+SQRT(B6*B6+4.*C6))
                    EVAP(I)=SIGT*AFAC*REVAP
                    WATER(I)=REVAP*REVAP
                    !
                    !    ***  CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER     ***
                    !    ***              HYDROSTATIC APPROXIMATION                 ***
                    !   
                    IF(I.EQ.1)then
                    
                        continue 
                        
                    ELSE
                        DHDP=(H(I)-H(I-1))/(P(I-1)-P(I))
                        DHDP=MAX(DHDP,10.0)
                        MP(I)=100.*GINV*LV(I)*SIGD*EVAP(I)/DHDP
                        MP(I)=MAX(MP(I),0.0)
                        !
                        !   ***   ADD SMALL AMOUNT OF INERTIA TO DOWNDRAFT              ***
                        !
                        FAC=20.0/(PH(I-1)-PH(I))
                        MP(I)=(FAC*MP(I+1)+MP(I))/(1.+FAC)
                        !   
                        !    ***      FORCE MP TO DECREASE LINEARLY TO ZERO                 ***
                        !    ***      BETWEEN ABOUT 950 MB AND THE SURFACE                  ***
                        !
                        IF(P(I).GT.(0.949*P(1)))THEN
                            JTT=MAX(JTT,I)
                            MP(I)=MP(JTT)*(P(1)-P(I))/(P(1)-P(JTT))
                        END IF              
                    END IF
                    !
                    !    ***       FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT     ***
                    !
                    IF(I.EQ.INB) THEN
                
                        continue
                    
                    Else
                
                        IF(I.EQ.1)THEN
                            QSTM=QS(1)
                        ELSE
                            QSTM=QS(I-1)
                        END IF
                
                        IF(MP(I).GT.MP(I+1))THEN
                    
                            RAT=MP(I+1)/MP(I)
                            QP(I)=QP(I+1)*RAT+Q(I)*(1.0-RAT)+100.*GINV*SIGD*(PH(I)-PH(I+1))*(EVAP(I)/MP(I))
                            UP(I)=UP(I+1)*RAT+U(I)*(1.-RAT)
                            VP(I)=VP(I+1)*RAT+V(I)*(1.-RAT)
                            DO J=1,NTRA
                                TRAP(I,J)=TRAP(I+1,J)*RAT+TRAP(I,J)*(1.-RAT)
                            END DO
                        
                        ELSE
                            IF(MP(I+1).GT.0.0)THEN
                        
                                QP(I)=(GZ(I+1)-GZ(I)+QP(I+1)*(LV(I+1)+T(I+1)*(CL-CPD))&
                                      +CPD*(T(I+1)-T(I)))/(LV(I)+T(I)*(CL-CPD))
                                UP(I)=UP(I+1)
                                VP(I)=VP(I+1)
                                DO J=1,NTRA
                                    TRAP(I,J)=TRAP(I+1,J)
                                END DO
                            
                            END IF
                        END IF
                    
                        QP(I)=MIN(QP(I),QSTM)
                        QP(I)=MAX(QP(I),0.0)
                    
                    End IF
            End do do_four_hundred
      
      
            !  ***  CALCULATE SURFACE PRECIPITATION IN MM/DAY     ***
            !
            PRECIP=PRECIP+WT(1)*SIGD*WATER(1)*3600.*24000./(ROWL*G)
    
        End If if_four_hundred_five
    end subroutine copes_precip
    
end module PRECIDOWN
