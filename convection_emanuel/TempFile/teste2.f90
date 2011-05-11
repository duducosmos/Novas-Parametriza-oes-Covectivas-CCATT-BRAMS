module t1
    implicit none
    private
    real,parameter :: a1=45         
    real, parameter :: e1=a1,f1=3.0d+0
    public e1,f1
    
end module t1

module t2
    use t1
    implicit none
    contains
        subroutine Oi(x)
            real x
            print*,x
            print*,e1,f1
        end subroutine Oi
end module t2

module t5
    
    use t1
    use t2
    implicit none
    private
    real, parameter:: g1=e1,h1=1d-3
    public Oi2
    contains
        subroutine Oi2(x)
            real x
            print*,'Oi2'
            call Oi(x)
        end subroutine Oi2
    
end module t5

program teste2
    use t5
    implicit none
    real a,b,c
    integer i
    b=0
    call Oi2(b)
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
