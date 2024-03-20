subroutine laser_tensor()
!************************************************************************
!  Discussion:
!   Input parameters: REAL(8)::ME%epsx
!
!   INTEGER(4)::ME%K and INTEGER(4)::ME%P are quantum number of the
!   laser field spherical tensor.
!   Only a few components of the tensor are none zero.
!
    use systemparameter
    implicit none
    !real(8),parameter::pi = 3.141592653589793D+00
    integer(4)::i,a,b,row
    row = (2 + 1)*(2 + 1)
    allocate(ME%tensor(row))
    Do 10,a = 0,2
        ME%k = a !Do 10,ME%k = 0,2
    Do 20,b = -ME%k,ME%k
        ME%p = b !Do 20 ME%p = -ME%k,ME%k
        i = ME%k**2 + ME%p + ME%k + 1
        If (ME%k == 2) then
            If (ME%p == -2) then
                ME%tensor(i) = 0.25*ME%epsx**2
                !write(*,'(2x,g14.6,2x,g14.6)')ME%k,ME%p
            Else If (ME%p == 0) then
                ME%tensor(i) = (0.5/(6**0.5))*(2-3*ME%epsx**2)
            Else If (ME%p == 2) then
                ME%tensor(i) = 0.25*ME%epsx**2
            !Else If (ME%p == 1) then
            !    ME%tensor(i) = -ME%epsx*SQRT(1-ME%epsx**2)
            !Else If (ME%p == -1) then
            !    ME%tensor(i) = ME%epsx*SQRT(1-ME%epsx**2)
            Else
                ME%tensor(i) = 0
            End If
        !Else if (ME%k == 0) then
        !    ME%tensor(i) = -0.5*3**(-0.5)
        Else
            ME%tensor(i) = 0
        End If
    20 continue
    10 continue
    !deallocate(ME%tensor)
    return
end subroutine

subroutine polar_tensor()
!************************************************************************
!  Discussion:
!   Input parameters: REAL(8)::ME%aa, ME%bb, ME%cc
!
!   INTEGER(4)::ME%K and INTEGER(4)::ME%Q are quantum number of the
!   polarization spherical tensor.
!   Only a few components of the tensor are none zero.
!
    use systemparameter
    implicit none
    !real(8),parameter::pi = 3.141592653589793D+00
    integer(4)::i,a,b,row
    row = 3*3
    allocate(ME%polar(row))
    !Do 10 ME%k = 0,2
    Do 10,a = 0,2
        ME%k = a
    !Do 20 ME%q = -ME%k,ME%k
    Do 20,b = -ME%k,ME%k
        ME%q = b
        i = ME%k**2 + ME%q + ME%k + 1
        If (ME%k == 2) then
            If (ME%q == -2)then
                ME%polar(i) = (ME%bb - ME%cc)/2!(ME%aa - ME%bb)/2 !
            Else If (ME%q == 0) then
                ME%polar(i) = (2*ME%aa - ME%cc - ME%bb)/(6**0.5) !(2*ME%cc - ME%bb - ME%aa)/(6**0.5)!
            Else If (ME%q == 2) then
                ME%polar(i) = (ME%bb - ME%cc)/2!(ME%aa - ME%bb)/2!
                !write ( *, '(2x,g14.6)')ME%polar(i)(ME%aa - ME%bb)/2 !
            Else
                ME%polar(i) = 0
            End If
        !Else if (ME%k == 0) then
        !    ME%polar(i) = - (ME%aa + ME%bb + ME%cc)/(3**0.5)
        Else
            ME%polar(i) = 0
        End If
    20 continue
    10 continue
    !deallocate(ME%polar)
    return
end subroutine

subroutine interaction1(i,passto)
!************************************************************************
!  Discussion:
!   Input parameter: INTEGER(4)::I
!   Output parameter: passto
!   Must call SUBROUTINE WINGER(...) and use the value of REAL(8)::WIN
!
!   To calculate the component of the SUBROUTINE INTERACTION(HAMI) that
!   gives the value of ME%K, ME%P, ...
!   It will repeated for 9 times when being called by subroutine interaction(hami)
!
    use systemparameter
    implicit none
    integer(4)::i
    real(8)::passto,win
    call Clebsch_Gordan(ME%J1,ME%k,ME%J2,ME%M1,ME%p,-ME%M2,win)
    !passto = ME%tensor(i)*win*((-1)**(ME%J1 - ME%k + ME%M2))* &
    !    ((2*ME%J2 + 1)**(0.5))
    passto = ME%tensor(i)*win
    !If (win.NE.0.0) then
    !write(*,'(I5,I5,I5,I5,I5,I5,F9.3,F9.3)')ME%J1,ME%k,ME%J2,ME%M1,ME%p,-ME%M2,&
        !ME%tensor(i),win
    !End If
end subroutine

subroutine interaction2(i,passto)
!************************************************************************
!  Discussion:
!   Input parameter: INTEGER(4)::I
!   Output parameter: passto
!
!   To calculate the component of the SUBROUTINE INTERACTION(HAMI)
!
    use systemparameter
    implicit none
    integer(4)::i
    real(8)::passto,win
    call Clebsch_Gordan(ME%J1,ME%k,ME%J2,ME%K1,ME%q,-ME%K2,win)
    !passto = ME%polar(i)*win*((-1)**(ME%J1 - ME%k + ME%K2))* &
    !    ((2*ME%J2 + 1)**(0.5))
    passto = ME%polar(i)*win
    !write(*,'(I5,I5,I5,I5,I5,I5,F9.3,F9.3)')ME%J1,ME%k,ME%J2,ME%K1,ME%q,-ME%K2,&
        !ME%polar(i),win
    !write(*,'(2x,g14.6)')ME%polar(i)
end subroutine

subroutine interaction(hami)
!************************************************************************
!  Discussion:
!   Input parameter: NONE
!   Output parameter: hami
!
!   To calculate the whole interaction term by giving the values of ME%K
!   and ME%P, ME%Q to SUBTOUTINE INTERACTION1, SUBTOUTINE INTERACTION1
!
    use systemparameter
    implicit none
    integer(4)::i,j
    integer(4)::a,b
    real(8)::hami1,hami2,hami,passto
!***************************************************************************
!   When 'K1' or 'K2' is changed, 'hami, hami1, hami2' will be setted zero.
!   But before the reinitialization, the hami term would be conveyed back to
!   the 'component' subroutine and be summed for 1/9/25 times there.
    hami1 = 0
    hami2 = 0
    hami = 0
    passto = 0
    !Do 10 ME%k = 0,2
    Do 10,a = 0,2
        ME%k = a
        !write(*,'(2x,g14.6)')ME%k
    Do 20,b = -ME%k,ME%k
        ME%p = b !Do 20 ME%p = -ME%k,ME%k
        i = ME%k**2 + ME%p + ME%k + 1 !+p
        call interaction1(i,passto)
        hami1 = hami1 + passto!*(-1)**(-ME%p)
        !write(*,'(2x,g14.6)')passto
    20 continue
    passto = 0
    Do 30,b = -ME%k,ME%k
        ME%q = b
        j = ME%k**2 - ME%q + ME%k + 1 !\alpha^k_{-q}
        call interaction2(j,passto)
        hami2 = hami2 + passto*(-1)**(-ME%q) !
        !write(*,'(2x,g14.3,2x,g9.3,2x,g9.3)')hami2,passto,ME%polar(j)
    30 continue
        hami = hami + hami1*hami2*(-1)**(ME%k)
        !write(*,'(2x,g14.6,g14.6)')hami1,hami2
    10 continue
    !write(*,'(2x,g14.6)')hami
end subroutine
