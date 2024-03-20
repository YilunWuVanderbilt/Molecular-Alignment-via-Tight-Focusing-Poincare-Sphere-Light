subroutine mean_cosin(t,m,Ut,CZC,CXC,CZA,CXA)
    use systemparameter
    implicit none

    integer(4)::m,a,b,c,d,e,f,i,j
    integer(4)::I1,J1
    real(8)::t
    complex(8)::Ut(m)
    complex(8)::CZC,CXC,CZA,CXA
    !real(8)::cos_zc,cos_xc,cos_za,cos_xa
    real(8),allocatable,dimension(:)::WJ1(:)
    real(8)::cos_zc1,cos_xc1,cos_za1,cos_xa1
    !real(8)::partition
    complex(8), parameter::ZI = (0.0, 1.0)
    CZC = (0.0D0,0.0D0)
    CXC = (0.0D0,0.0D0)
    CZA = (0.0D0,0.0D0)
    CXA = (0.0D0,0.0D0)
    
    !call partition_function(partition)
    
    Do 10 a = 0,ME%NUM
        ME%J1 = a
        allocate(ME%UR((2*a + 1),(2*a + 1))) !2*J1 + 1
        allocate(WJ1(2*a + 1))
        call matrixelements(a) !pass J1 to matrix***()
        call diagonalization(a) !pass J1 to the matrix defined in diag***()
        deallocate(ME%matrix)
        DO J1 = 1,(2*a + 1)
            DO I1 = 1,(2*a + 1)
                ME%UR(I1,J1) = ME%VR(I1,J1)
                !write(*, "(2x,g14.6)")ME%VR(I1,ME%tau1+J1)
            END DO
            WJ1(J1) = ME%WJtau(J1)
        END DO
        deallocate(ME%VR)
        deallocate(ME%WJtau)
    Do 20 b = -ME%J1,ME%J1
        ME%tau1 = b
    Do 30 c = -ME%J1,ME%J1
        ME%M1 = c
    DO 40 d = 0,ME%NUM
        ME%J2 = d
        !allocate(UR2((2*d + 1),(2*d + 1)))
        call matrixelements(d) !matrixelements(ME%J2)
        call diagonalization(ME%J2)
        deallocate(ME%matrix)
    DO 50 e = -ME%J2,ME%J2
        ME%tau2 = e
    Do 60 f = -ME%J2,ME%J2
        ME%M2 = f
!*****************************************************************************************
!   Assign values to 'J1, tau1, M1, J2, tau2, M2' that will be used in 'KSUM'.
!   And with 'J1, tau1, M1' unchanged, 'KSUM' will be called for 35 times.
        call KSUM(cos_zc1,cos_xc1,cos_za1,cos_xa1)
        j = (2)*(ME%J1-1)*ME%J1*(2*ME%J1-1)/3+2*(ME%J1-1)**2+ &
            3*(ME%J1-1)+1+(ME%tau1+ME%J1)*(2*ME%J1+1)+ME%M1+ME%J1+1
        i = (2)*(ME%J2-1)*ME%J2*(2*ME%J2-1)/3+2*(ME%J2-1)**2+ &
            3*(ME%J2-1)+1+(ME%tau2+ME%J2)*(2*ME%J2+1)+ME%M2+ME%J2+1
        
        CZC = CZC + cos_zc1*REAL(CONJG(Ut(j))*Ut(i))* &
            exp(ZI*(WJ1(ME%tau1 + ME%J1 + 1)-ME%WJtau(ME%tau2 + ME%J2 + 1))*t/1)
        CXC = CXC + cos_xc1*REAL(CONJG(Ut(j))*Ut(i))* &
            exp(ZI*(WJ1(ME%tau1 + ME%J1 + 1)-ME%WJtau(ME%tau2 + ME%J2 + 1))*t/1)
        CZA = CZA + cos_za1*REAL(CONJG(Ut(j))*Ut(i))* &
            exp(ZI*(WJ1(ME%tau1 + ME%J1 + 1)-ME%WJtau(ME%tau2 + ME%J2 + 1))*t/1)
        CXA = CXA + cos_xa1*REAL(CONJG(Ut(j))*Ut(i))* &
            exp(ZI*(WJ1(ME%tau1 + ME%J1 + 1)-ME%WJtau(ME%tau2 + ME%J2 + 1))*t/1)
        
    60 continue
    50 continue
        deallocate(ME%VR)
        deallocate(ME%WJtau)
    40 continue
    30 continue
    20 continue
        deallocate(ME%UR)
        deallocate(WJ1)
    10 continue
end subroutine

subroutine KSUM(cos_zc1,cos_xc1,cos_za1,cos_xa1)
    use systemparameter
    implicit none
    integer(4)::c,d!,N1,N2
    real(8)::cos_zc1,cos_xc1,cos_za1,cos_xa1
    real(8)::cos_zc,cos_xc,cos_za,cos_xa
    cos_zc1 = 0.0D+00
    cos_xc1 = 0.0D+00
    cos_za1 = 0.0D+00
    cos_xa1 = 0.0D+00

    DO 30 c = -ME%J1,ME%J1
        ME%K1 = c !DO ME%K1 = -J1,J1
    Do 40 d = -ME%J2,ME%J2
        ME%K2 = d !Do ME%K2 = -J2,J2
        call cosine(cos_zc,cos_xc,cos_za,cos_xa)
!   To sum 'K1,K2' respectively in this loop
!   And leave 'L1, tau1, M1, L2, tau2, M2' to get values from the 'motion' program
        cos_zc1 = cos_zc1 + ME%UR(ME%K1+ME%J1+1,ME%tau1+ME%J1+1)* &
            ME%VR(ME%K2+ME%J2+1,ME%tau2+ME%J2+1)*cos_zc
        cos_xc1 = cos_xc1 + ME%UR(ME%K1+ME%J1+1,ME%tau1+ME%J1+1)* &
            ME%VR(ME%K2+ME%J2+1,ME%tau2+ME%J2+1)*cos_xc
        cos_za1 = cos_za1 + ME%UR(ME%K1+ME%J1+1,ME%tau1+ME%J1+1)* &
            ME%VR(ME%K2+ME%J2+1,ME%tau2+ME%J2+1)*cos_za
        cos_xa1 = cos_xa1 + ME%UR(ME%K1+ME%J1+1,ME%tau1+ME%J1+1)* &
            ME%VR(ME%K2+ME%J2+1,ME%tau2+ME%J2+1)*cos_xa
        !cos_zc1 = cos_zc1 + ME%UR(ME%tau1+ME%J1+1,ME%K1+ME%J1+1)* &
        !    ME%VR(ME%tau2+ME%J2+1,ME%K2+ME%J2+1)*cos_zc
        !cos_xc1 = cos_xc1 + ME%UR(ME%tau1+ME%J1+1,ME%K1+ME%J1+1)* &
        !    ME%VR(ME%tau2+ME%J2+1,ME%K2+ME%J2+1)*cos_xc
        !cos_za1 = cos_za1 + ME%UR(ME%tau1+ME%J1+1,ME%K1+ME%J1+1)* &
        !    ME%VR(ME%tau2+ME%J2+1,ME%K2+ME%J2+1)*cos_za
        !cos_xa1 = cos_xa1 + ME%UR(ME%K1+ME%J1+1,ME%tau1+ME%J1+1)* &
        !    ME%VR(ME%tau2+ME%J2+1,ME%K2+ME%J2+1)*cos_xa
    40 continue
    30 continue
end subroutine

subroutine cosine(cos_zc,cos_xc,cos_za,cos_xa)
    use systemparameter
    implicit none
    real(8)::cos_zc,cos_xc,cos_za,cos_xa
    real(8)::M0,K0,M2,M_2,K2,K_2
    real(8)::win
    integer(4)::cond

    IF ((ME%J1.EQ.ME%J2).AND.(ME%K1.EQ.ME%K2).AND.(ME%M1.EQ.ME%M2))THEN
        cond = 1
    Else
        cond = 0
    End IF

    call Clebsch_Gordan(ME%J1,2,ME%J2,ME%M1,0,ME%M2,win)
    M0 = win*((-1)**(ME%J1-2+ME%M2))*sqrt(REAL(2*ME%J2+1))
    
    call Clebsch_Gordan(ME%J1,2,ME%J2,ME%K1,0,ME%K2,win)
    K0 = win*((-1)**(ME%J1-2+ME%K2))*sqrt(REAL(2*ME%J2+1))
        
    call Clebsch_Gordan(ME%J1,2,ME%J2,ME%M1,2,ME%M2,win)
    M2 = win*((-1)**(ME%J1-2+ME%M2))*sqrt(REAL(2*ME%J2+1))
    
    call Clebsch_Gordan(ME%J1,2,ME%J2,ME%M1,-2,ME%M2,win)
    M_2 = win*((-1)**(ME%J1-2+ME%M2))*sqrt(REAL(2*ME%J2+1))
    
    call Clebsch_Gordan(ME%J1,2,ME%J2,ME%K1,2,ME%K2,win)
    K2 = win*((-1)**(ME%J1-2+ME%K2))*sqrt(REAL(2*ME%J2+1))
    
    call Clebsch_Gordan(ME%J1,2,ME%J2,ME%K1,-2,ME%K2,win)
    K_2 = win*((-1)**(ME%J1-2+ME%K2))*sqrt(REAL(2*ME%J2+1))
        
    cos_zc = (1.0/3.0)*cond + (2.0/3.0)*sqrt(REAL(2*ME%J1 + 1)/REAL(2*ME%J2 + 1))*M0*K0
    cos_xc = (1.0/3.0)*cond - sqrt(1.0/6.0)*sqrt(REAL(2*ME%J1 + 1)/REAL(2*ME%J2 + 1))*&
        K0*(sqrt(2.0/3.0)*M0 - (M2 + M_2))
    cos_za = (1.0/3.0)*cond - sqrt(1.0/6.0)*sqrt(REAL(2*ME%J1 + 1)/REAL(2*ME%J2 + 1))*&
        M0*(sqrt(2.0/3.0)*K0 - (K2 + K_2))
    cos_xa = (1.0/3.0)*cond + (1.0/4.0)*sqrt(REAL(2*ME%J1 + 1)/REAL(2*ME%J2 + 1))*&
        (sqrt(2.0/3.0)*M0 - (M2 + M_2))*(sqrt(2.0/3.0)*K0 - (K2 + K_2))
end subroutine

!subroutine partition_function(partition)
!    use systemparameter
!    implicit none
!    integer(4)::a,b,spin
!    real(8)::partition
!    partition = 0.0D0
!
!    Do 10 a = 0,ME%NUM
!        ME%J1 = a !Do 10 ME%J1 = 0,ME%NUM
!        call matrixelements(a) !pass J1 to matrix***()
!        call diagonalization(a) !pass J1 to the matrix defined in diag***()
!        deallocate(ME%matrix)
!    Do 20 b = -ME%J1,ME%J1
!        call selfspin(a,b,spin)
!        ME%tau1 = b !Do 20 ME%tau1 = -ME%J1,ME%J1
!        partition = partition + (2*ME%J1 + 1)*spin*exp(-ME%WJtau(ME%tau1 + ME%J1 + 1)/(1.3794*4.0D6))
!    20 continue
!        write(*, "(2x,g14.6)")partition
!        deallocate(ME%VR)
!        deallocate(ME%WJtau)
!    10 continue
!    return
!end subroutine
!
!subroutine distribution(ZC1,count,n)
!    use systemparameter
!    implicit none
!    integer(4)::J0,Tau0,M0,count,n,i,j,I1,J1,spin
!    real(8)::partition,DZC
!    real(8),allocatable,dimension(:)::WJ1(:)
!    !real(8),allocatable,dimension(:,:)::ZC1(:,:)
!    real(8)::ZC1(count,n)
!    
!    !write(*,*)count
!    call partition_function(partition)
!    DO 10 j = 1,count
!    DZC = 0.0
!    DO 20 J0 = 0,ME%NUM
!    DO 30 Tau0 = -J0,J0
!        call selfspin(J0,Tau0,spin)
!        allocate(ME%UR((2*J0 + 1),(2*J0 + 1))) !2*J1 + 1
!        allocate(WJ1(2*J0 + 1))
!        call matrixelements(J0) !pass J1 to matrix***()
!        call diagonalization(J0) !pass J1 to the matrix defined in diag***()
!        deallocate(ME%matrix)
!        DO J1 = 1,(2*J0 + 1)
!            DO I1 = 1,(2*J0 + 1)
!                ME%UR(I1,J1) = ME%VR(I1,J1)
!                !write(*, "(2x,g14.6)")ME%VR(I1,ME%tau1+J1)
!            END DO
!            WJ1(J1) = ME%WJtau(J1)
!        END DO
!        deallocate(ME%VR)
!        deallocate(ME%WJtau)
!    DO 40 M0 = -J0,J0
!    i = (2)*(J0-1)*J0*(2*J0-1)/3+2*(J0-1)**2+3*(J0-1)+1+(Tau0+J0)*(2*J0+1)+M0+J0+1
!    DZC = DZC + ZC1(j,i)*spin*exp(-WJ1(Tau0 + J0 + 1)/(1.3794*4.0D6))
!    !write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' )i,count
!    40 continue
!        deallocate(ME%UR)
!        deallocate(WJ1)
!    30 continue
!    20 continue
!    DZC = DZC/partition
!    write (10, '(2x,g14.6,2x,g14.6)' )DZC
!    10 continue
!    !deallocate(ZC1)
!end subroutine
!
!subroutine selfspin(J,Tau,g)
!    use systemparameter
!    implicit none
!    integer(4)::J,Tau,Km1,Kp1,g
!
!    IF(Tau .LE. 0)THEN
!        IF(mod(J + Tau,2) .EQ. 0)THEN
!            Km1 = (J + Tau)/2
!            Kp1 = (J + Tau)/2 - Tau
!        ELSE
!            Km1 = (J + Tau + 1)/2 + Tau
!            Kp1 = (J + Tau + 1)/2
!        END IF
!    ELSE
!        IF(mod(J - Tau,2) .EQ. 0)THEN
!            Km1 = (J - Tau)/2
!            Kp1 = (J - Tau)/2
!        ELSE
!            Km1 = (J - Tau + 1)/2 + Tau
!            Kp1 = (J - Tau + 1)/2
!        END IF
!    END IF
!    
!    IF((mod(Km1,2) .EQ. 0).AND.(mod(Kp1,2) .EQ. 0))THEN
!        g = 7
!    ElSE
!        g = 3
!    END IF
!end subroutine
