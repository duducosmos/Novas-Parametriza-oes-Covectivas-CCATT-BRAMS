program teste2
    implicit none
    real a,b,c
    integer i
    b=0
    
    do i=1,10
        if(i == 1) then 
            continue
        else
            a=1
            b=b+a
            c=b*b
            print*,i
        end if
    end do
end program teste2
