       program teste
       real a,b,c
       integer i
        b=2
        c=3
        a=1
       do 10 i=1,5
           b=b*c
           write(*,*) a,b
           a=b
 10    continue
       end program     
