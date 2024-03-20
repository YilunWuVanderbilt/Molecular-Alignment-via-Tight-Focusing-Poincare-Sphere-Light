subroutine rk4vec_test()
    use systemparameter
    implicit none

    integer(4)::n !Total C_JtM
    integer(4)::i,j,count,CON
    real(8),parameter::dt = 0.0001!0.00005
    real(8)::norm
    external rk4vec_test_f
    real(8)::t0,t1,ti,t_start
    real(8), parameter::tmax = 0.06!0.06!0.2!0.025! !evolving time in ps =60ps
    complex(8),allocatable,dimension(:)::u0,ut
    !real(8),allocatable,dimension(:)::u0,u1
    complex(8)::CZC,CXC,CZA,CXA,CXB
    real(8)::CYB,CYA,CYC,CZB
    integer(4)::J0,Tau0,M0
    real(8),allocatable,dimension(:,:)::ZC(:,:)
   
    n = 2*(ME%NUM*(ME%NUM + 1)*(2*ME%NUM + 1))/3 + 2*ME%NUM**2 + 3*ME%NUM + 1
    !write (*,"(I5)")n
    allocate(u0(n),ut(n))
    t_start = -0.005!-0.005!-0.01! !initial time in ps
!****************************************************************************
!   count the total steps
!
    ti = t_start !initial time in ps
    count = 0
    do
        if ( tmax <= ti ) then
            exit
        end if
        count = count + 1
        t1 = ti + dt
        ti = t1
    end do
    !write(*,*)count
    CON = count
    allocate(ZC(1:count,1:n))
    !count = 0
    t0 = t_start
!****************************************************************************
!   inital wave function
!
    DO 10 J0 = 0,ME%NUM
    DO 20 Tau0 = -J0,J0
    DO 30 M0 = -J0,J0
        !J0 = 1
        !Tau0 = -1
        !M0 = -1
        i = (2)*(J0-1)*J0*(2*J0-1)/3+2*(J0-1)**2+3*(J0-1)+&
            1+(Tau0+J0)*(2*J0+1)+M0+J0+1
        t0 = t_start
        count = 0
        u0(1:n) = (0.0D+00,0.0D+00)
        u0(i) = (1.0D+00,0.0D+00)
        write(*,*)i
    do
        !write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' )t0,u0(1)
        call mean_cosin(t0,n,u0,CZC,CXC,CZA,CXA)
        CYC = 1 - REAL(CXC) - REAL(CZC)
        CYA = 1 - REAL(CXA) - REAL(CZA)
        CYB = 1 - CYA - CYC
        CZB = 1 - REAL(CZC) - REAL(CZA)
        CXB = 1 - CYB - CZB
        write(*,'(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)')t0*1000,REAL(CZC)
        !t0*1000,(REAL(CYA)-REAL(CXA),(REAL(CXB)-REAL(CYB))
        !t0*1000,REAL(CXB),REAL(CYA),REAL(CZC)
        !t0*1000,REAL(CZC),REAL(CYB),REAL(CXA)
        !t0*1000,REAL(CZB),REAL(CYA),REAL(CXC)
        !REAL(CYA)-REAL(CXA),REAL(CXB)-REAL(CYB)

!   Stop if we've exceeded TMAX.
!
        if ( tmax <= t0 ) then
            exit
        end if
!
!   Otherwise, advance to time T1, and have RK4 estimate
!   the solution U1 there.
!
        count = count + 1
        !write(*,*)count
        ZC(count,i) = REAL(CZC)
        !write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' )ZC(count,i)
        t1 = t0 + dt
        call rk4vec(t0,n,u0,dt,rk4vec_test_f,ut)
        !write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' )t1,u0(1)
!   Shift the data to prepare for another step.
        norm = 0.0D+00
        DO j = 1,n
            norm = norm + REAL(CONJG(ut(j))*ut(j))
        END DO
        !write(*,*)norm
        t0 = t1
        u0(1:n) = ut(1:n)/norm
        !If(t0 >= tmax/100)stop
    end do
    30 continue
    20 continue
    10 continue
    write(*,*)count

    !DO i = 1,n
    !write(11,'(2x,g14.6)')ZC(:,i)
    !END DO
    !call distribution(ZC,count,n)
    deallocate(ZC)
    deallocate(u0,ut)
    return
end subroutine

subroutine rk4vec_test_f(t,n,u,uprime)
    use systemparameter
    implicit none

    integer(4)::n,a,b,c,d,e,f,i,j
    integer(4)::I1,J1
    real(8)::inter_term
    real(8)::t,T_tau
    complex(8), parameter::ZI = (0.0, 1.0)
    complex(8)::u(n)
    complex(8)::uprime(n)
    real(8),allocatable,dimension(:)::WJ1(:)
    
    uprime(1:n) = (0.0D+00,0.0D+00)
    If (t.LT.ME%t_zero) then
        T_tau = ME%tau_rise
    Else
        T_tau = ME%tau_fall
    End If
    Do 10 a = 0,ME%NUM
        ME%J1 = a !Do 10 ME%J1 = 0,ME%NUM
        allocate(ME%UR((2*a + 1),(2*a + 1))) !2*J1 + 1
        allocate(WJ1(2*a + 1))
        call matrixelements(a) !pass J1 to matrix***()
        call diagonalization(a) !pass J1 to the matrix defined in diag***()
        !write(*,*)"--------------------------------------------"
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
        ME%tau1 = b !Do 20 ME%tau1 = -ME%J1,ME%J1
    Do 30 c = -ME%J1,ME%J1
        ME%M1 = c !Do 30 ME%M1 = -ME%J1,ME%J1
    DO 40 d = 0,ME%NUM
        ME%J2 = d !DO 40 ME%J2 = 0,ME%NUM
        !allocate(UR2((2*d + 1),(2*d + 1)))
        call matrixelements(d) !matrixelements(ME%J2)
        call diagonalization(ME%J2)
        deallocate(ME%matrix)
        !DO J1 = 1,(2*d + 1)
        !   DO I1 = 1,(2*d + 1)
        !       UR2(I1,J1) = ME%VR(I1,J1)
        !   END DO
        !END DO
        !deallocate(ME%VR)
    DO 50 e = -ME%J2,ME%J2
        ME%tau2 = e !DO 50 ME%tau2 = -ME%J2,ME%J2
    Do 60 f = -ME%J2,ME%J2
        ME%M2 = f !Do 60 ME%M2 = -ME%J2,ME%J2
!*****************************************************************************************
!   Assign values to 'J1, tau1, M1, J2, tau2, M2' that will be used in 'component'.
!   And with 'J1, tau1, M1' unchanged, 'component' will be called for 35 times.
        !ME%epsx = 0 !新加的一行！！！！！！！!
        !call laser_tensor() !新加的一行！！！！！！！!
        call component(inter_term)
        j = (2)*(ME%J1-1)*ME%J1*(2*ME%J1-1)/3+2*(ME%J1-1)**2+ &
            3*(ME%J1-1)+1+(ME%tau1+ME%J1)*(2*ME%J1+1)+ME%M1+ME%J1+1
        i = (2)*(ME%J2-1)*ME%J2*(2*ME%J2-1)/3+2*(ME%J2-1)**2+ &
            3*(ME%J2-1)+1+(ME%tau2+ME%J2)*(2*ME%J2+1)+ME%M2+ME%J2+1
        !!uprime(j) = uprime(j) + inter_term*u(i)*(ME%eps0*exp(-((t-ME%t_zero)/(2*T_tau))**2)
        uprime(j) = uprime(j) + ((ZI)/1)*inter_term*u(i)* &
        exp(ZI*(WJ1(ME%tau1 + ME%J1 + 1) - ME%WJtau(ME%tau2 + ME%J2 + 1))*t/1)* &
        ((ME%eps0)*exp(-((t - ME%t_zero)/(T_tau))**2))
        !deallocate(ME%tensor)!新加的一行！！！！！！！!
        !ME%epsx = 1 !新加的一行！！！！！！！
        !call laser_tensor() !新加的一行！！！！！！！!
        !call component(inter_term)!新加的一行！！！！！！！!
        !uprime(j) = uprime(j) + ((ZI)/1)*inter_term*u(i)* &
        !exp(ZI*(WJ1(ME%tau1 + ME%J1 + 1) - ME%WJtau(ME%tau2 + ME%J2 + 1))*t/1)* &
        !    (ME%eps0*exp(-0.5*((t - 0.0085)/(0.00025))**2)) !新加的一行！！！！！！！!
        !deallocate(ME%tensor)!新加的一行！！！！！！！!
        !write ( *, '(I5,I5,I5,I5,I5,I5,I5,I5)' )j,a,b,c,i,d,e,f
        !write(*,*)"--------------------------------------------"
        !write ( *, '(2x,g9.3,2x,g9.3)' )((-ZI)/1)*inter_term*exp(ZI* &
        !    (WJ1(ME%tau1 + ME%J1 + 1)-ME%WJtau(ME%tau2 + ME%J2 + 1))* &
        !    t/1)*(ME%eps0*exp(-((t-ME%t_zero)/(2*T_tau))**2))**2
        !write ( *, *)
        !write ( *, '(g10.3,g10.3)' )exp(ZI*(WJ1(ME%tau1 + ME%J1 + &
        !    1)-ME%WJtau(ME%tau2 + ME%J2 + 1)))
    60 continue
    50 continue
        deallocate(ME%VR)
        deallocate(ME%WJtau)
    40 continue
        !write ( *,'(I5,I5,I5,I5,I5,I5,I5,I5)')j,a,b,c,i,d,e,f
        !write ( *, '(2x,g9.3,2x,g9.3)' )uprime(j)
        !stop
    30 continue
    20 continue
        deallocate(ME%UR)
        deallocate(WJ1)
    10 continue
    return
end subroutine

subroutine rk4vec(t0,m,u0,dt,f,u)
!************************************************************************
!  Parameters:
!
!    Input, real(8)::t0, the current time.
!
!    Input, integer(8)::m, the dimension of the system.
!
!    Input, real(8)::u0(m), the solution estimate at the current time.
!
!    Input, real(8)::dt, the time step.
!
!    Input, external f, a subroutine of the form
!      subroutine f(t,m,u,uprime)
!    which evaluates the derivative uprime(1:m) given the time T and
!    solution vector u(1:m).
!
!    Output, real(8)::u(m), the fourth-order Runge-Kutta solution
!    estimate at time t0+dt.
!
    use systemparameter
    implicit none

    integer(4)::m

    real(8)::dt
    external f
    complex(8)::f0(m),f1(m),f2(m),f3(m)
    !real(8)::f0(m),f1(m),f2(m),f3(m)
    real(8)::t0,t1,t2,t3
    complex(8)::u0(m),u1(m),u2(m),u3(m),u(m)
    !real(8)::u0(m),u1(m),u2(m),u3(m),u(m)
!
!   Get four sample values of the derivative.
!
    !write(*,"(I5)")m
    call f(t0,m,u0,f0)
    !write ( *, '(2x,g9.3)' )f0
    t1 = t0 + dt/2.0D+00
    u1(1:m) = u0(1:m) + dt*f0(1:m)/2.0D+00
    call f(t1,m,u1,f1)

    t2 = t0 + dt/2.0D+00
    u2(1:m) = u0(1:m) + dt*f1(1:m)/2.0D+00
    call f(t2,m,u2,f2)
!    write ( *, '(2x,g9.3)' )u2

    t3 = t0 + dt
    u3(1:m) = u0(1:m) + dt*f2(1:m)
    call f(t3,m,u3,f3)
!    write ( *, '(2x,g9.3)' )u3
!
!   Combine them to estimate the solution U at time T1.
!
    u(1:m) = u0(1:m) + (dt/6.0D+00)*( &
                   f0(1:m) &
       + 2.0D+00 * f1(1:m) &
       + 2.0D+00 * f2(1:m) &
       +           f3(1:m))

    return
end subroutine
