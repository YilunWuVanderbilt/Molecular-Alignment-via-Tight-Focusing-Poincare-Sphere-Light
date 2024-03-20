subroutine ifactr(Factorial,num)
    implicit none
    integer(4)::Factorial,num,k,J
    J = 1
    If(num.EQ.0) GO TO 10
    Do k = 1,num
        J = J*k
    end do
    10 Factorial = J
    !write(*,"(I5,/)")Factorial
end subroutine

!Calculate Sum over t of (-1)^t/x
subroutine sum_t(j1,j2,j3,m1,m2,st)
    implicit none
    integer(4)::j1,j2,j3,m1,m2,t1,t2,t3,t4,t5
    integer(4)::tmin,tmax,t,Factorial,sumt
    real(8)::st
    sumt = 1
    st = 0.D0
    t1 = j2 - j3 - m1
    t2 = j1 - j3 + m2
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2
    tmin = min( 0,  min( t1, min(t2, min(t3, min(t4, t5)))))
    tmax = max( t1, max(t2, max(t3, max(t4, t5))))
    Do t = tmin,tmax
    IF(((t.LT.0).OR.(t - t1).LT.0).OR.((t - t2).LT.0).OR.((t3 - t).LT.0)&
            .OR.((t4 - t).LT.0).OR.((t5 - t).LT.0)) GO TO 11
        call ifactr(Factorial,t)
        sumt = sumt*Factorial
        call ifactr(Factorial,t - t1)
        sumt = sumt*Factorial
        call ifactr(Factorial,t - t2)
        sumt = sumt*Factorial
        call ifactr(Factorial,t3 - t)
        sumt = sumt*Factorial
        call ifactr(Factorial,t4 - t)
        sumt = sumt*Factorial
        call ifactr(Factorial,t5 - t)
        sumt = sumt*Factorial
        st = st + float((-1)**t)/float(sumt)
    11 continue
    end do
    !write(*,"(f9.3,/)")st
end subroutine

!Calculate the triangle_coeff idelta
subroutine idelta(j1,j2,j3,delta)
    implicit none
    integer(4)::j1,j2,j3,Factorial
    integer(4)::F1,F2,F3,F4
    real(8)::delta
    delta = 0.d0
    call ifactr(Factorial,j1 + j2 - j3)
    F1 = Factorial
    !write(*,"(I5,/)")Factorial
    call ifactr(Factorial,j1 - j2 + j3)
    F2 = Factorial
    call ifactr(Factorial,-j1 + j2 + j3)
    F3 = Factorial
    call ifactr(Factorial,j1 + j2 + j3 + 1)
    F4 = Factorial
    delta = sqrt(float(F1*F2*F3)/float(F4))
    !write(*,"(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)")delta
end subroutine

!Calculate the sign and the factorials in the square root with delta
subroutine remainJM(j1,j2,j3,m1,m2,m3,rem)
    implicit none
    integer(4)::I1,IY,Factorial
    integer(4)::j1,j2,j3,m1,m2,m3
    real(8)::rem
    I1 = (-1)**(j1 - j2 - m3) !the sign
    call ifactr(Factorial,j1 + m1)
    IY = Factorial
    call ifactr(Factorial,j1 - m1)
    IY = IY*Factorial
    call ifactr(Factorial,j2 + m2)
    IY = IY*Factorial  
    call ifactr(Factorial,j2 - m2)
    IY = IY*Factorial
    call ifactr(Factorial,j3 + m3)
    IY = IY*Factorial  
    call ifactr(Factorial,j3 - m3)
    IY = IY*Factorial
    !write(*,"(I5,/)")Factorial
    rem = I1*sqrt(float(IY))
end subroutine

subroutine winger(j1,j2,j3,m1,m2,m3,win)
    implicit none
    integer(4)::j1,j2,j3,m1,m2,m3
    real(8)::st,delta,rem,win
    win = 0.D0
    IF( IABS(M1) .GT. J1 ) GOTO 973
    IF( IABS(M2) .GT. J2 ) GOTO 973
    IF( IABS(M3) .GT. J3 ) GOTO 973
    IF( IABS(J1-J2) .GT. J3 ) GOTO 973
    IF( J3 .GT. J1+J2 ) GOTO 973
    IF( M1+M2+M3 .NE. 0) GOTO 973
    call sum_t(j1,j2,j3,m1,m2,st)
    call idelta(j1,j2,j3,delta)
    call remainJM(j1,j2,j3,m1,m2,m3,rem)
    win = rem*delta*st
    973 CONTINUE
end subroutine
