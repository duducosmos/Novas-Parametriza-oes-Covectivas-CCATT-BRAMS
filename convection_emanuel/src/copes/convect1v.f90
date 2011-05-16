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

                  
module copes !>Convective parametrization based in K. A. Emanuel (1991,1999) scheme
    
     
    use VarDeclared
    use CopesStartEnd
    use GEO
    use TLIFT
    use INSTA
    use GEO
    use LEVEL
    use ADI
    implicit none
    private   
 
    
    public copes_convection,copes_TLIFT


    contains

!
!**************************************************************************************************
!**************************************************************************************************
!**************************************************************************************************
!

        subroutine copes_convection(T1,   Q1,    QS1,     U1,    V1,      TRA1,    P1,    PH1,&
                              NL, DELT, IFLAG,  FT,     FQ,   FU,&
                              FV,  FTRA, PRECIP, WD,   TPRIME, QPRIME, CBMF)    
            implicit none 
            real, dimension(:) :: T1,Q1,QS1,U1,V1,P1,PH1
            real, dimension(:,:) :: TRA1,FTRA
   
            real, dimension(:) :: FT,FQ,FU,FV
            
            integer :: NL,IFLAG
            real :: DELT,PRECIP,WD,TPRIME,QPRIME,CBMF


            integer :: IHMIN,NK
            real :: AHMIN, AHMAX
            integer :: i,j,k,l
           
            !> The first automatic test of this
            !> parametrization, difine if  
            !> T1_EntryMatriz will be .TRUE. 
            !> or .FALSE.
            call copes_set_matriz_at1(T1,Q1,QS1,U1,V1,TRA1,P1,PH1)  


            if( T1_EntryMatriz .eqv. .TRUE.) then
                DELTI=1.0/DELT

                !>           ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***

                FT=0.0;FQ=0.0;FU=0.0;FV=0.0;FTRA=0.0 !> Equivale a do 5 de convect43c.f

                do_seven: DO  I=1,NL+1
                              RDCP=(RD*(1.-Q(I))+Q(I)*RV)/(CPD*(1.-Q(I))+Q(I)*CPV)
                              TH(I)=T(I)*(1000.0/P(I))**RDCP !>Temperatura Potencial do ar umido
                          enddo do_seven

                PRECIP=0.0; WD=0.0; TPRIME=0.0; QPRIME=0.0; IFLAG=0

                if(IPBL /= 0 ) then
                   call copes_AdiabaticAdjustment()
                endif
                
                !CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY
                call copes_Geo_Heat_SEnergy() 


                !Find that model level below the level of minimum moist !
                call copes_Level() 

!>
!>
!>
!>
!>
                call copes_free_matriz!>
            else
                print*,'An erro in the input parameters'
                continue
            end if
            

        end subroutine  copes_convection

        
end module  copes


