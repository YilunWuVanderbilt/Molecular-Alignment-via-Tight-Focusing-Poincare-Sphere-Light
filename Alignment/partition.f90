subroutine distribution()
    use systemparameter
    implicit none
    integer(4)::J0,Tau0,M0,count,n,i,j,I1,J1!,spin
    integer(4)::nline,ioS,p
    real(8)::partition,DZC,spin!make change here
    
    real(8),allocatable,dimension(:)::WJ1(:),a(:)
    real(8),allocatable,dimension(:,:)::ZC1(:,:)
    open(10,file="out.txt")
    nline=0;
    open(101,file="wfunc.txt",status='old',action='read')
    do
        read(101,*,ioStat=ioS)
        if(ioS/=0)Exit
        nline=nline+1
    end do
    
    n = 2*(ME%NUM*(ME%NUM + 1)*(2*ME%NUM + 1))/3 + 2*ME%NUM**2 + 3*ME%NUM + 1
    count = nline/n
    !write ( *, * )n,count,nline
    allocate(ZC1(count,n))
    allocate(a(nline))
    close(101,status='keep')
    open(101,file="wfunc.txt",status='old',action='read')
    DO j = 1, n
        DO i = 1, count
            p = (j-1)*count + i
            !write ( *, * )i,j,p
            read(101,*,ioStat=ioS)a(p)
            ZC1(i,j) = a(p)
            !write ( *, '(2x,g14.6)' )a(p)
            if(ioS/=0)Exit
        END DO
    END DO
    close(101,status='keep')
    
    call partition_function(partition)
    DO 10 j = 1,count
    DZC = 0.0
    DO 20 J0 = 0,ME%NUM
    DO 30 Tau0 = -J0,J0
        call selfspin(J0,Tau0,spin)
        allocate(ME%UR((2*J0 + 1),(2*J0 + 1))) !2*J1 + 1
        allocate(WJ1(2*J0 + 1))
        call matrixelements(J0) !pass J1 to matrix***()
        call diagonalization(J0) !pass J1 to the matrix defined in diag***()
        deallocate(ME%matrix)
        DO J1 = 1,(2*J0 + 1)
            DO I1 = 1,(2*J0 + 1)
                ME%UR(I1,J1) = ME%VR(I1,J1)
                !write(*, "(2x,g14.6)")ME%VR(I1,ME%tau1+J1)
            END DO
            WJ1(J1) = ME%WJtau(J1)
        END DO
        deallocate(ME%VR)
        deallocate(ME%WJtau)
    DO 40 M0 = -J0,J0
        i = (2)*(J0-1)*J0*(2*J0-1)/3+2*(J0-1)**2+3*(J0-1)+1+(Tau0+J0)*(2*J0+1)+M0+J0+1
DZC = DZC + ZC1(j,i)*spin*exp(-WJ1(Tau0 + J0 + 1)*4.79924/(100*200))
        !write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' )ZC1(j,i)
    40 continue
        deallocate(ME%UR)
        deallocate(WJ1)
    30 continue
    20 continue
    DZC = DZC/partition
    write (10, '(2x,g14.6,2x,g14.6)' )DZC
    10 continue
    !deallocate(ZC1)
    close(10,status='keep')
end subroutine

subroutine partition_function(partition)
    use systemparameter
    implicit none
    integer(4)::a,b!,spin!make change here
    real(8)::partition,spin
    partition = 0.0D0
    
    Do 10 a = 0,ME%NUM
        ME%J1 = a !Do 10 ME%J1 = 0,ME%NUM
        call matrixelements(a) !pass J1 to matrix***()
        call diagonalization(a) !pass J1 to the matrix defined in diag***()
        deallocate(ME%matrix)
    Do 20 b = -ME%J1,ME%J1
        call selfspin(a,b,spin)
        ME%tau1 = b !Do 20 ME%tau1 = -ME%J1,ME%J1
partition = partition + (2*ME%J1 + 1)*spin*exp(-ME%WJtau(ME%tau1 + ME%J1 + 1)*4.79924/(100*200))
    20 continue
        write(*, "(2x,g14.6)")partition
        deallocate(ME%VR)
        deallocate(ME%WJtau)
    10 continue
    return
end subroutine

subroutine selfspin(J,Tau,g)
    use systemparameter
    implicit none
    integer(4)::J,Tau,Km1,Kp1!,g
    real(8)::g
    
    IF(Tau .LE. 0)THEN
        IF(mod(J + Tau,2) .EQ. 0)THEN
            Km1 = (J + Tau)/2
            Kp1 = (J + Tau)/2 - Tau
        ELSE
            Km1 = (J + Tau + 1)/2 + Tau
            Kp1 = (J + Tau + 1)/2
        END IF
    ELSE
        IF(mod(J - Tau,2) .EQ. 0)THEN
            Km1 = (J - Tau)/2
            Kp1 = (J - Tau)/2
        ELSE
            Km1 = (J - Tau + 1)/2 + Tau
            Kp1 = (J - Tau + 1)/2
        END IF
    END IF
    
    IF((mod(Km1,2) .EQ. 0).AND.(mod(Kp1,2) .EQ. 0))THEN
        g = 7
    ElSE
        g = 3
    END IF
end subroutine
