program main
    use systemparameter
    implicit none
    integer(4)::i
    !character(19) fmt
    !fmt = '(F7.2,"+",F7.2,"i")'

    open(14,file='data1.dat',status='replace')
    
    CO%A = 1
    CO%NA = 0.95
    CO%lambda = 532D-9
    CO%nr = 1
    CO%pi = dacos(-1.D0)
    CO%k = 2*CO%pi/CO%lambda
    CO%cali = 1D-6
    
    !CO%lb = -1D-6
    !CO%ub = 1D-6
    CO%No = 100
    CO%d = (ASIN(CO%NA/CO%nr) - 0)/CO%No
    
    CO%z = 0
    
    CO%Row = 101
    CO%Col = 101

    CO%m = -5
    CO%n = -1

    CO%w0 = CO%d*64/225
    CO%w1 = CO%d*(161+13*sqrt(17.5))/900
    CO%w2 = CO%d*(161-13*sqrt(17.5))/900

    CO%y1p = sqrt(245.0 - 14*sqrt(70.0))/21 
    CO%y1m = -sqrt(245.0 - 14*sqrt(70.0))/21  
    CO%y2p = sqrt(245.0 + 14*sqrt(70.0))/21  
    CO%y2m = -sqrt(245.0 + 14*sqrt(70.0))/21
    call plotting()
    !write(*,"(101(E15.7,E15.7),/)")(CO%matrix(:,i),i=1,CO%Col)
    write(14,"(101E15.7)")(CO%Intensity(:,i),i=1,CO%Col)
    close(14,status='keep')
end program
