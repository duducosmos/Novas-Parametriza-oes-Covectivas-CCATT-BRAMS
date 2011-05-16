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

!>
!>************************************************************************************************************************************************
!>
!> PERFORM DRY ADIABATIC ADJUSTMENT
!>
!>************************************************************************************************************************************************
!>
module ADI
    use VarDeclared
    implicit none
    contains 
    
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
end module ADI
