subroutine Clebsch_Gordan(L1,L2,L,m1,m2,m,CG_coeff)
!Program Clebsch_Gordan
  !  subroutine Clebsch_Gordan(L1,L2,L,m1,m2,m,CG_coeff)
  !
  !  Calculate the Wigner of Clebsch-Gordan 3J symbol values for all spherical harmonics up
  !  to L1=6, L2=6 in
  !
  !
  !
    integer :: i, L, k, L1, m1, L2, m2, r, m, L_upper, L_lower, phase(0:40)!, j,&
      !ii(30), kk(30), mm1(30), mm2(30), nterms, primes(77), iii, jjj, i1, power, itxt, jtxt
    double precision :: fact(0:40), norm, CG_coeff !CGCG_coeff(30), sum, sumc
    !character*1 :: neg
    !character :: line1*60, line2*60
    !data primes /2,3,5,7,11, 13,17,19,23,29, 31,37,41,43,47, 51, 53, 59, 61, 67,71, &
    !73, 79, 83, 87, 89, 97,101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, &
    !157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, &
    !241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, &
    !347, 349, 353, 359, 367, 373, 379/
      fact(0) = 1.d0
      phase(0) = 1
      do i = 1, 40
        fact(i) = fact(i - 1)*i
        phase(i) = -phase(i - 1)
      end do
      !open(unit=7, file="CG coefficients.txt")
      !open(unit=11, file="out.txt")
      !do L1 = 0, 6
      !  do L2 = 0, L1
!
!   Calculate all Clebsch-Gordan coefficients for spherical harmonics L1 and L2 up to 5 each
!
          !L1 = 6
          !L2 = 6
          !L = 11
          !m = -11
          !m1 = -5
          !m2 = -6

          L_upper = L1 + L2
          L_lower = Abs(L1 - L2)
          if((L.GE.L_lower).AND.(L.LE.L_upper).AND.((m1).EQ.(m - m2)))then
          !do L = L_lower, L_upper
          !  do m = -L, L
          !    nterms = 0
          !    do m2 = -L2, L2
          !      m1 = m - m2
          !      if (m1 > L1 .or. m1 < -L1) cycle
                norm = 0.d0
!
!  Now follows the Wigner formula for the 3J symbol
!
                do r = 0, 20
                  if (L2 + L - r - m1 < 0) cycle
                  if (L - m - r       < 0) cycle
                  if (L1 - m1 - r     < 0) cycle
                  if (L2 - L + m1 + r < 0) cycle
                  norm = norm + phase(L1 + r - m1)*fact(L1 + m1 + r)*fact(L2 + L - r - m1)/ &
                  (fact(r)*fact(L - m - r)*fact(L1 - m1 - r)*fact(L2 - L + m1 + r))
                end do
                CG_coeff = norm * &
                sqrt(fact(L + m)*fact(L - m)*fact(L1 - m1)*fact(L2 - m2)*fact(L1 + L2 - L)*(2*L + 1)/ &
                  (fact(L1 + m1)*fact(L2 + m2)*fact(L1 - L2 + L)*fact(L2 - L1 + L)*fact(L1 + L2 + L + 1)))
                !write(11,*)L1,L2,L,m,m1,m2,CG_coeff
                if (Abs(CG_coeff) > 1.d-10) then
                  norm = 1.d0/CG_coeff**2
                  do i = 1,245025
                    k = nint(norm*i)
                    if (Abs(k - norm*i) < 1.d-8) then
                      exit
                    end if
                  end do
                else
                  i = 0
                  k = 1
                end if
                !nterms = nterms + 1
                !mm1(nterms) = m1
                !mm2(nterms) = m2
                !ii(nterms) = i
                !kk(nterms) = k
                !CGCG_coeff(nterms) = CG_coeff
              !end do
          Else
              CG_coeff=0.0
          End If
          !write(11,*)L1,L2,L,m,m1,m2,CG_coeff
      !      end do
      !    end do
      !  end do
      !end do
  end subroutine
!end program Clebsch_Gordan
