module systemparameter
implicit none
type com
    real(8)::lb,ub,lambda,cali,k,r,phi,z,y1p,y1m,y2p,y2m
    real(8)::pi,NA,nr
    integer(4)::order,No,m,n,Row,Col
    real(8)::d,d2
    real(8)::w0                                                       
    real(8)::w1
    real(8)::w2
    integer(4)::A
    complex(8)::matrix(101,101)
    real(8)::Intensity(101,101)
end type
type(com)::CO
end module
