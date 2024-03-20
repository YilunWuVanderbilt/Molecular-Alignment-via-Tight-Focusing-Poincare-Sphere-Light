module systemparameter
    implicit none
    type matrixegin
        integer(4)::J1,J2,K1,K2,M1,M2,tau1,tau2,k,p,q,NUM
        real(8)::eps0,aa,bb,cc,F,G,H,A,B,C
        real(8)::epsx
        real(8)::t_zero,tau_rise,tau_fall
        real(8),allocatable,dimension(:)::tensor,polar,WJtau
        real(8),allocatable,dimension(:,:)::matrix,VR,UR
    end type    
type(matrixegin)::ME
end module
