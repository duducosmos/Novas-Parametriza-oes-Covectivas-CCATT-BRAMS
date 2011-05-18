program AskMe
    use copes_help
    integer :: choose
    
    print*,"Well come to the AskMe Program"
    print*,"I have some informations about the "
    print*,"copes package: Convective parametrization based in the Emanuel scheme."
    print*," "
    print*," "
    print*,"This package are composed by the following modules: "
    print*,"1) ADI, 2) ADJUSTTENDENCI,  3) copes, 4) CopesStartEnd, &
    5) INSTA,  6) NENTRE, 7) TENDENCI,8) VarDeclared, 9) ADJUST, 10) AMFLUX, &
    11) copes_help, 12) GEO, 13) LEVEL, 14) PRECIDOWN,15) TENDENCI2, 16) TLIFT,"

    print*,"type the equivalent number for information about the respective module"   
    print*," "
    print*," "            
    print*,"If you want to stop this help press 0"
    print*," "
    print*," " 
    read*, choose
    
    do_one: do while( choose /= 0)
        
        select case(choose)
            case (1)
                call help_AID
            case (2)
                call help_ADJUSTTENDENCI
            case (3)
                call help_copes
            case (4)
                call help_CopesStartEnd
            case (5)
                call help_INSTA
            case (6)
                call help_NENTRE
            case (7)
                call help_TENDENCI
            case (8 )
                call help_VarDeclared
            case (9) 
                call help_ADJUST
            case (10)
                call help_AMFLUX
            case (11)
                call help_copes_help
            case (12)
                call help_GEO
            case (13)
                call help_LEVEL
            case (14)
                call help_PRECIDOWN
            case (15) 
                call help_TENDENCI2
            case (16) 
                call help_TLIFT
            case DEFAULT
                print*,"You did not make a choose"
                print*," "
                print*," "
                print*,"1) ADI, 2) ADJUSTTENDENCI,  3) copes, 4) CopesStartEnd, &
                        5) INSTA,  6) NENTRE, 7) TENDENCI,8) VarDeclared, 9) ADJUST, 10) AMFLUX, &
                        11) copes_help, 12) GEO, 13) LEVEL, 14) PRECIDOWN,15) TENDENCI2, 16) TLIFT,"
    
                print*,"type the equivalent number for information about the respective module"
                print*," "
                print*," "            
                print*,"If you want to stop this help press 0"
                
            
        end select
        print*," "
        print*," " 
        print*,"1) ADI, 2) ADJUSTTENDENCI,  3) copes, 4) CopesStartEnd, &
                        5) INSTA,  6) NENTRE, 7) TENDENCI,8) VarDeclared, 9) ADJUST, 10) AMFLUX, &
                        11) copes_help, 12) GEO, 13) LEVEL, 14) PRECIDOWN,15) TENDENCI2, 16) TLIFT,"
    
        print*,"type the equivalent number for information about the respective module"
        print*," "
        print*," "    
        print*,"Tell me your choose."
        read*, choose
        
    end do do_one


end program AskMe
