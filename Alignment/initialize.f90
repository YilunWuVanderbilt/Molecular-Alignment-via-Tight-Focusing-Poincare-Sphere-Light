subroutine initialize()
    use systemparameter
    implicit none
    real(8)::kappa
    ME%A = 144.7397987224*6.283!1.998*30*6.283!1.2291490778*6.283! ! GHz
    ME%B = 30.0092250458*6.283!1.998*30*6.283!0.8693981282*6.283!! GHz
    ME%C = 24.8228155224*6.283!0!3.1178415632*6.283!! GHz
    kappa = (2*ME%B - ME%A - ME%C)/(ME%A - ME%C)
    ME%F = 0.5*(kappa - 1)  !GHz
    ME%G = 1 !GHz
    ME%H = -0.5*(kappa + 1) !GHz
    ME%aa = 1.692!0!121.4*0.0148!*1.113!121.4*1.6487772754! !A**3
    ME%bb = 1.066!0!63.2*0.0148!*1.113! 63.2*1.6487772754! !A**3
    ME%cc = 1.491!1.710!163.9*0.0148!*1.113!*8.988 !A**3
    !ME%epsx = (0,-0.383)
    ME%epsx = 0.0!0.383!0.383!0.383! 0.66!
    ME%eps0 = 1483.0*12!*10!33.79*40!1214.749!2.749*200!
    ME%t_zero = 0.00 !
    ME%tau_rise = 0.0006!0.0006!0.06!600fs
    ME%tau_fall = 0.0006!0.06! !600fs
    ME%NUM = 2
    !write ( *, '(2x,g9.3,2x,g9.3,2x,g9.3)' )ME%tau_rise,ME%tau_fall,ME%t_zero
end subroutine
