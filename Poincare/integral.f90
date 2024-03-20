subroutine variation(X,Y,p_order)
    use systemparameter
    real(8)::X,Y
    integer(4)::p_order
    !write(*,"(I5,/)")p_order
    X = 0 + CO%d*(p_order + 0.5) + 0.5*CO%d*Y
end subroutine
      
subroutine functionF(Y,func,p_order)
    use systemparameter  
    real(8)::X,Y
    complex(8)::func
    integer(4)::p_order
    !write(*,"(I5,/)")p_order
    call variation(X,Y,p_order)
    func = sin(X)*sqrt(cos(X))*exp(-(sin(X)/CO%NA)**2)*(cos(X/2)**2)&
      *(cmplx(cos(CO%pi/2*(CO%m+1)),sin(CO%pi/2*(CO%m+1)),8)*cmplx(cos(CO%m*CO%phi),sin(CO%m*CO%phi),8)&
      *bessel_jn(CO%m,CO%k*CO%r*sin(X))+cmplx(cos(CO%pi/2*(CO%n+1)),sin(CO%pi/2*(CO%n+1)),8)&
      *cmplx(cos(CO%n*CO%phi),sin(CO%n*CO%phi),8)*bessel_jn(CO%n,CO%k*CO%r*sin(X)))&
      -(sin(X/2)**2)*(cmplx(cos(CO%pi/2*(CO%m+3)),sin(CO%pi/2*(CO%m+3)),8)&
      *cmplx(cos((CO%m+2)*CO%phi),sin((CO%m+2)*CO%phi),8)&
*bessel_jn(CO%m+2,CO%k*CO%r*sin(X))+cmplx(cos(CO%pi/2*(CO%n-1)),sin(CO%pi/2*(CO%n-1)),8)&
      *cmplx(cos((CO%n-2)*CO%phi),sin((CO%n-2)*CO%phi),8)*bessel_jn(CO%n-2,CO%k*CO%r*sin(X)))
end subroutine

subroutine functionF2(Y,func,p_order)
    use systemparameter
    real(8)::X,Y
    complex(8)::func
    integer(4)::p_order

    call variation(X,Y,p_order)
    func = sin(X)*sqrt(cos(X))*exp(-(sin(X)/0.95)**2)*(cos(X/2)**2)&
      *(cmplx(cos(CO%pi/2*(CO%m)),sin(CO%pi/2*(CO%m)),8)*cmplx(cos(CO%m*CO%phi),sin(CO%m*CO%phi),8)&
      *bessel_jn(CO%m,CO%k*CO%r*sin(X))-cmplx(cos(CO%pi/2*(CO%n)),sin(CO%pi/2*(CO%n)),8)&
      *cmplx(cos(CO%n*CO%phi),sin(CO%n*CO%phi),8)*bessel_jn(CO%n,CO%k*CO%r*sin(X)))&
      +(sin(X/2)**2)*(cmplx(cos(CO%pi/2*(CO%m+2)),sin(CO%pi/2*(CO%m+2)),8)&
      *cmplx(cos((CO%m+2)*CO%phi),sin((CO%m+2)*CO%phi),8)&
      *bessel_jn(CO%m+2,CO%k*CO%r*sin(X))-cmplx(cos(CO%pi/2*(CO%n-2)),sin(CO%pi/2*(CO%n-2)),8)&
      *cmplx(cos((CO%n-2)*CO%phi),sin((CO%n-2)*CO%phi),8)*bessel_jn(CO%n-2,CO%k*CO%r*sin(X)))
end subroutine

subroutine functionF3(Y,func,p_order)
    use systemparameter
    real(8)::X,Y
    complex(8)::func
    integer(4)::p_order

    call variation(X,Y,p_order)
    func = sin(X)*sqrt(cos(X))*exp(-(sin(X)/0.95)**2)&
      *(cmplx(cos(CO%pi/2*(CO%m+2)),sin(CO%pi/2*(CO%m+2)),8)&
      *cmplx(cos((CO%m+1)*CO%phi),sin((CO%m+1)*CO%phi),8)&
      *bessel_jn(CO%m+1,CO%k*CO%r*sin(X))+cmplx(cos(CO%pi/2*(CO%n)),sin(CO%pi/2*(CO%n)),8)&
      *cmplx(cos((CO%n-1)*CO%phi),sin((CO%n-1)*CO%phi),8)&
      *bessel_jn(CO%n-1,CO%k*CO%r*sin(X)))
end subroutine 

subroutine IntegralSum(I,I2,I3,func)
    use systemparameter
    complex(8)::I,I2,I3,func
    integer(4)::j
    I = cmplx(0,0)
    I2 = cmplx(0,0)
    I3 = cmplx(0,0)
    !write(*,"(f9.3,/)")p_r
    jloop: Do j = 1,CO%No
        call functionF(0.0d0,func,j)
        I = I + CO%w0*func
        call functionF2(0.0d0,func,j)
        I2 = I2 + CO%w0*func
        call functionF3(0.0d0,func,j)
        I3 = I3 + CO%w0*func

        call functionF(CO%y1p,func,j)
        I = I + CO%w1*func
        call functionF2(CO%y1p,func,j)
        I2 = I2 + CO%w1*func
        call functionF3(CO%y1p,func,j)
        I3 = I3 + CO%w1*func

        call functionF(CO%y1m,func,j)
        I = I + CO%w1*func
        call functionF2(CO%y1m,func,j)
        I2 = I2 + CO%w1*func
        call functionF3(CO%y1m,func,j)
        I3 = I3 + CO%w1*func

        call functionF(CO%y2p,func,j)
        I = I + CO%w2*func
        call functionF2(CO%y2p,func,j)
        I2 = I2 + CO%w2*func
        call functionF3(CO%y2p,func,j)
        I3 = I3 + CO%w2*func
        
        call functionF(CO%y2m,func,j)
        I = I + CO%w2*func
        call functionF2(CO%y2m,func,j)
        I2 = I2 + CO%w2*func
        call functionF3(CO%y2m,func,j)
        I3 = I3 + CO%w2*func
        
        !write(*,"(E15.7,E15.7)")I
        !if (j == 2) then
        !    exit jloop
        !end if
    end do jloop
end subroutine

subroutine plotting()
    use systemparameter
    integer(4)::j,l
    complex(8)::I,I2,I3,func
    real(8)::x,y
    lloop: Do l = 1,CO%Col
        jloop: Do j = 1,CO%Row
            x = (l - 51.d0)/(CO%Row - 1)*2.d0*CO%cali
            y = (j - 51.d0)/(CO%Col - 1)*2.d0*CO%cali
            CO%r = sqrt(x*x + y*y)
            if(x>0.and.y>=0)then
                CO%phi = ATAN(y/x)
            else if(x<0)then
                CO%phi = ATAN(y/x) + CO%pi
            else if(x>0.and.y<0)then
                CO%phi = ATAN(y/x) + 2*CO%pi
            else
                CO%phi = ATAN(y/x)
            end if
            !write(*,"(E15.7)")CO%r
            call IntegralSum(I,I2,I3,func)
            !write(*,"(E15.7,E15.7)")I
            CO%matrix(j,l) = sqrt(I*I+I2*I2+I3*I3)
            CO%Intensity(j,l) = (real(CO%matrix(j,l)))**2+(aimag(CO%matrix(j,l)))**2
            !if (j == 1) then
            !    exit lloop
            !end if
        end do jloop
    end do lloop
!write(*,"(E15.7,E15.7)")(CO%matrix(j,:),j=1,CO%Col)
end subroutine
