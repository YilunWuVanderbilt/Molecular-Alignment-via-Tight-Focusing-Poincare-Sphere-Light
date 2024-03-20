subroutine matrixelements(J12)
    use systemparameter
    implicit none
    integer(4)::row,col
    integer(4)::J12
    integer(4)::a,b
    row = 2*J12+1
    col = 2*J12+1
    allocate(ME%matrix(row,col))
    !write (*,"(I5)")J12
    Do 10 a = -J12,J12 !a is col, representing the quantum no.K'
    Do 20 b = -J12,J12 !b is row, representing the quantum no.K
        If(a.EQ.b)then
            ME%matrix(b+J12+1,a+J12+1) = ME%F*(J12*(J12+1)-&
                a**2)+ME%G*a**2
        !write(*,"(F9.3,F9.3)")ME%F,ME%matrix(b+J12+1,a+J12+1)
        !write(*,"(I5,I5)")a,b
        Else If(b.EQ.(a + 2))then
            ME%matrix(b+J12+1,a+J12+1) = ME%H*(0.25*(J12*(J12+1)-&
                a*(a+1))*(J12*(J12+1)-(a+1)*(a+2)))**0.5
            !IF((J12.EQ.2).AND.(ME%J2.EQ.2))&
            !write(*,"(5F9.3)")(ME%matrix(b+J12+1,a+J12+1))
            !write(*,"(I5,I5)")a,b
        Else If(b.EQ.(a - 2))then
            ME%matrix(b+J12+1,a+J12+1) = ME%H*(0.25*(J12*(J12+1)-&
                a*(a-1))*(J12*(J12+1)-(a-1)*(a-2)))**0.5
        Else
            ME%matrix(b+J12+1,a+J12+1) = 0
        End if
    20 continue
    10 continue
    !If (J12.EQ.2) then
    !    open(10,file="out.txt")
    !    write(10,"(5F9.3)")(ME%matrix(i,:),i=1,row)
    !    write(10,*)
    !End if
end subroutine

subroutine diagonalization(J12)
    use systemparameter
    implicit none
 
    integer(4)::J12
    integer(4)::N,I,J
    integer(4)::LDA,LDVL,LDVR
    integer(4)::LWMAX
    integer(4)::INFO,LWORK   
    !real(8)::hami,inter_term
    real(8),allocatable,dimension(:,:)::A(:,:),VL(:,:),VR(:,:),&
        WR(:),WI(:),WORK(:)
    N = 2*J12+1 
    LDA = N
    LDVL = N
    LDVR = N
    LWMAX = 1000
    allocate(A(LDA, N),VL(LDVL,N),VR(LDVR,N),WR(N),WI(N),&
        WORK(LWMAX))
    LWORK = -1
    allocate(ME%VR(LDVR,N))
    allocate(ME%WJtau(N))
    !call matrixelements(J12)
    DO J = 1,N
        DO I = 1,LDA
            A(I,J) = ME%matrix(I,J)
        END DO
    END DO

    CALL DGEEV('N','Vectors',N,A,LDA,WR,WI,VL,LDVL,&
        VR,LDVR,WORK,LWORK,INFO)
    LWORK = MIN(LWMAX,INT(WORK(1)))
    CALL DGEEV('N','Vectors',N,A,LDA,WR,WI,VL,LDVL,&
        VR,LDVR,WORK,LWORK,INFO)
    IF( INFO.GT.0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
    END IF
    DO J = 1,N
    ME%WJtau(J) = 0.5*(ME%A + ME%C)*J12*(J12 + 1) + 0.5*(ME%A - ME%C)*WR( J )
    !WRITE(*,"(2x,g14.6)") ME%WJtau(J)
    END DO
    !DO J = 1, N
    !    IF( WI( J ).EQ.0 ) THEN
    !        WRITE(*,9998) WR( J )
    !        write(*,'(I5)') J
    !        WRITE(*,9998) 0.5*(ME%A + ME%C)*J12*(J12 + 1)
    !    ELSE
    !        WRITE(*,9999) WR( J ), WI( J )
    !    END IF
    !END DO
    !WRITE(*,*)
    !9998 FORMAT( 11(:,1X,F14.9) )
    !9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
    DO J = 1,N
        DO I = 1,LDA
            ME%VR(I,J) = VR(J,I) !row vector -> column vector
        END DO
    END DO
    !DO I = 1, N
    !    J = 1
    !    DO WHILE( J.LE.N )
    !        WRITE(*,9998) VR( I, J )
    !        J = J + 1
    !    END DO
    !    WRITE(*,*)
    !END DO
    !9998 FORMAT( 11(:,1X,F6.2) )

    deallocate(A,VL,VR,WR,WI,WORK)
end subroutine

subroutine component(inter_term)
    use systemparameter
    implicit none
    integer(4)::c,d!,N1,N2
    real(8)::hami,inter_term
!********************************************************************************
!   When any number from 'L1, tau1, M1, L2, tau2, M2' is changed, the interaction
!   term will be setted zero. But before the reinitialization, this term would be
!   conveyed back to the 'motion' program and be summed for 1/10/35 times there.
    inter_term = 0.0D+00
    !N1 = 2*ME%J1 + 1
    !N2 = 2*ME%J2 + 1

    DO 30 c = -ME%J1,ME%J1
        ME%K1 = c !DO ME%K1 = -J1,J1
    Do 40 d = -ME%J2,ME%J2
        ME%K2 = d !Do ME%K2 = -J2,J2
!********************************************************************************
!   With 'J1, K1, M1, J2, K2, M2' assigned with values
!   'Interaction(hami)' will be called for 5*5(*35) times
        !write(*,'(I5,I5)')c,d
        call interaction(hami)
        !write(*,*)
!   To sum 'K1,K2' respectively in this loop
        inter_term = inter_term + ME%UR(ME%K1+ME%J1+1,ME%tau1+ME%J1+1)* &
            ME%VR(ME%K2+ME%J2+1,ME%tau2+ME%J2+1)*((-1)**(ME%M2 - & !-
            ME%K2))*hami*SQRT(REAL(2*ME%J1+1)*REAL(2*ME%J2+1))
        
        !inter_term = inter_term + ME%UR(ME%tau1+ME%J1+1,ME%K1+ME%J1+1)* &
        !    ME%VR(ME%tau2+ME%J2+1,ME%K2+ME%J2+1)*((-1)**(ME%M2 - &
        !    ME%K2))*hami*SQRT(REAL(2*ME%J1+1)*REAL(2*ME%J2+1))
!   And leave 'L1, tau1, M1, L2, tau2, M2' to get values from the 'motion' program
        !write(*, "(2x,g14.6,2x,g14.6)") ME%UR(ME%K1+ME%J1+1,ME%tau1+ME%J1+1)* &
        !ME%VR(ME%K2+ME%J2+1,ME%tau2+ME%J2+1)*hami*((-1)**(ME%M2 - &
        !ME%K2))*((2*ME%J1+1)*(2*ME%J2+1))**(0.5),inter_term
    40 continue
    30 continue
    !write(*, "(2x,g14.6)")inter_term
end subroutine
