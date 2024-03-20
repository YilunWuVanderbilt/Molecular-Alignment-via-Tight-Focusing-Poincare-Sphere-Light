program main
    use systemparameter
    implicit none
    open(10,file="out.txt")
    open(11,file="wfunc.txt")
    call initialize()
    call laser_tensor()
    call polar_tensor()
    !call diagonalization(inter_term)
    call rk4vec_test()
    !call distribution()
    deallocate(ME%tensor)
    deallocate(ME%polar)
    close(11,status='keep')
    close(10,status='keep')
end
