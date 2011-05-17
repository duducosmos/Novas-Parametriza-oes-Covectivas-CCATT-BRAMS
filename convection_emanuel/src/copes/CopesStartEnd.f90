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




!> This module contains the subroutines to initialize and finish 
!> the copes_convection subroutine at module copes
module CopesStartEnd
    use VarDeclared
    implicit none
    logical :: T1_EntryMatriz !> This Variable is used by copem_set_matriz_at1, if the input parameters are O.K. the code continue and this varyable is set as .TRUE.
     

    contains
    
        !>*****************************************************************************************************************************************
        !>The Subroutine copem_set_matriz_at1 make a test of the dimension of the arrays  T1,Q1,QS1,P1,PH1,FT1,FQ1
        !> And allocate the equivalente arrays T,Q,QS,P,PH,FT,FQ
        !> These array must have same dimension ND. If its OK, all matriz are dimensioned.
        !>*****************************************************************************************************************************************
        
        subroutine copes_set_matriz_at1(T1,Q1,QS1,U1,V1,TRA1,P1,PH1)
        integer :: N1,N2,N3,N4,N5,N6,N7,N8,N9
        real, dimension(:) :: T1,Q1,QS1,U1,V1,P1,PH1
        real, dimension(:,:) :: TRA1
        integer, dimension(2) :: NDNTRA
        integer :: ierr1,ierr2,ierr3,ierr4,ierr5,ierr6,ierr7,ierr8

            N1=size(T1); N2=size(Q1); N3=size(QS1); N4=size(U1); N5=size(V1); N6=size(P1); N7=size(PH1)
            NDNTRA = shape(TRA1)
            N8=NDNTRA(1)
            N9=NDNTRA(2)

            NTRA = N9

            
            if(N1 == N2 .and. N2 == N3 .and. N3 == N4 .and. N4 == N5 .and. N5 == N6 .and. N6 == N8) then
                ND= N6
                if(N7 == N6 +1 ) then
                    NA=ND+1
                    T1_EntryMatriz = .TRUE.
                    allocate(T(ND),stat=ierr1)
                    if(ierr1 /= 0) print*, 'Allocated error T'
                    allocate(Q(ND),stat=ierr2)
                    if(ierr2 /= 0) print*, 'Allocated error Q'
                    allocate(QS(ND),stat=ierr3)
                    if(ierr3 /= 0) print*, 'Allocated error QS'
                    allocate(U(ND),stat=ierr4)
                    if(ierr4 /= 0) print*, 'Allocated error U'
                    allocate(V(ND),stat=ierr5)
                    if(ierr5 /= 0) print*, 'Allocated error V'
                    allocate(P(ND),stat=ierr6)
                    if(ierr6 /= 0) print*, 'Allocated error P'
                    allocate(PH(NA),stat=ierr7)
                    if(ierr7 /= 0) print*, 'Allocated error U'
                    allocate(TRA(ND,NTRA),stat=ierr8)
                    if(ierr8 /= 0) print*, 'Allocated error TRA'

                    allocate(TH(NA))
                    allocate(FTRA(ND,NTRA))
                    allocate(FQ(ND)); allocate(FU(ND)); allocate(FV(ND))
                    allocate(FT(ND)); allocate(UENT(NA,NA)); allocate(VENT(NA,NA)); allocate(TRAENT(NA,NA,NTRA))
                    allocate(TRATM(NA)); allocate(UP(NA)); allocate(TRAP(NA,NTRA)); allocate(M(NA)); allocate(MP(NA))
                    allocate(MENT(NA,NA)); allocate(QENT(NA,NA)); allocate(ELIJ(NA,NA)); allocate(SIJ(NA,NA))
                    allocate(TVP(NA)); allocate(TV(NA)); allocate(WATER(NA)); allocate(QP(NA)); allocate(EP(NA))
                    allocate(WT(NA)); allocate(EVAP(NA)); allocate(CLW(NA)); allocate(SIGP(NA)); allocate(TP(NA))
                    allocate(TOLD(NA)); allocate(CPN(NA)); allocate(LV(NA)); allocate(LVCP(NA)); allocate(H(NA))
                    allocate(HP(NA)); allocate(GZ(NA)); allocate(HM(NA));allocate(NENT(NA))




                    if(ierr1 /= 0 .or. ierr2 /= 0 .or. ierr3 /= 0 .or. ierr4 /= 0 &
                       .or. ierr5 /= 0 .or. ierr6 /= 0 .or. ierr7 /= 0 .or. ierr8 /= 0 ) then

                        print*,"Allocated error, verify the code"
                        T1_EntryMatriz = .FALSE.
                        call copem_free_matriz

                    endif
                   
                else
                    print*,"The Dimension of PH array must be equal",ND,&
                           "+1, but it has dimension of ",N7
                    print*,"The  copem_convection will do nothing until this correction be made."
                    print*,"Verify your datas"
                    print*,"Good look for you"
                    T1_EntryMatriz = .FALSE.

                endif



            else
                print*,"The Dimension of arrays T,Q,QS,U,V and P, and the numbers of lines of TRA, must be equals,",&
                       "but they  have, respecively, the following dimensions: ",N1,N2,N3,N4,N5,N6,N8
                print*,"The  copem_convection will do nothing until this correction be made."
                print*,"Verify your datas."
                print*,"Good look for you."
                T1_EntryMatriz = .FALSE.

            endif

            

        end subroutine copes_set_matriz_at1
        
        !>
        !>************************************************************************************************************************************
        !>
        !> This subroutine is used to deallocate arrays T,Q,QS,P,PH,FT,FQ and matriz TRA
        !>
        !>************************************************************************************************************************************
        !>
        
        subroutine copes_free_matriz
            deallocate(T); deallocate(Q); deallocate(QS); deallocate(U)
            deallocate(V); deallocate(P); deallocate(PH); deallocate(TRA)
            deallocate(TH); deallocate(FV); deallocate(FU); deallocate(FQ)
            deallocate(FT); deallocate(FTRA)

            deallocate(TH); deallocate(FTRA); deallocate(FQ); deallocate(FU); deallocate(FV)
            deallocate(FT); deallocate(UENT); deallocate(VENT); deallocate(TRAENT)
            deallocate(TRATM); deallocate(UP); deallocate(TRAP); deallocate(M); deallocate(MP)
            deallocate(MENT); deallocate(QENT); deallocate(ELIJ); deallocate(SIJ)
            deallocate(TVP); deallocate(TV); deallocate(WATER); deallocate(QP); deallocate(EP)
            deallocate(WT); deallocate(EVAP); deallocate(CLW); deallocate(SIGP); deallocate(TP)
            deallocate(TOLD); deallocate(CPN); deallocate(LV); deallocate(LVCP); deallocate(H)
            deallocate(HP); deallocate(GZ); deallocate(HM); deallocate(NENT)
           

        end subroutine copes_free_matriz

end module CopesStartEnd
