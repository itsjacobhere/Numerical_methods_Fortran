!______________________________________________________________________
!written by Jacob Turner
!Uses newtons method to approximate the root of a function within 
!a user specified precision.
!finished on 2/17/2018
!______________________________________________________________________
module dp
    integer, parameter :: wp = selected_real_kind(15)
end module dp
!______________________________________________________________________
module vars
    use dp
    !variable definitions:
    ! n = roughness 
    ! w = width (m)
    ! s = slope
    ! Q = flow rate (m^3/sec)
    real(wp),parameter :: n = 0.05, w = 1.0, s = 0.01, Q = 1.50
    real(wp),parameter :: c1 = ((w**5)*s**(3._wp/2._wp))/(n**3)
    real(wp),parameter :: c2 = -4*Q**3
    real(wp),parameter :: c3 = -(4*w)*(Q**3)
    real(wp),parameter :: c4 = -(w**2)*(Q**3)
end module vars
!______________________________________________________________________
    ! Estimate the zero of f(x) using Newton's method. 
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   x0: the initial guess
    !   debug: logical, prints iterations if debug=.true.
    ! Returns:
    !   the estimate x satisfying f(x)=0 (assumes Newton converged) 
    !   the number of iterations, iters
module newtonmod
    use dp
    use vars
    implicit none
    contains         
    subroutine Newton(f, fp, x0, x, maxiteration, iters, exitflag, epsi)
        implicit none
        !Variables passed between subroutine and main.
        real(wp), intent(in) :: x0
        real(wp), external :: f, fp
        integer, intent(in) :: exitflag
        real(wp), intent(out) :: x
        real(wp), intent(in) :: epsi
        integer, intent(out) :: iters
        integer, intent(inout) :: maxiteration
        ! Declare any local variables:
        real(wp) :: deltax, fx, fxprime
        integer :: k
        !-----------------------------------------
        ! initial guess
        x = x0
        if (exitflag == 0) then
            write(*,'(A,2x,f5.2)') "Initial guess:", x
        endif
        ! Newton iteration to find a zero of f(x) 
        do k=1,maxiteration
            ! evaluate function and its derivative:
            fx = f(x)
            fxprime = fp(x)
            if (abs(fx) < epsi) then
                exit  ! jump out of do loop
            endif
            ! compute Newton increment x:
            deltax = fx/fxprime
            ! update x:
            x = x - deltax
            if (exitflag == 0) then
                write(*,13) k,x
                13 format('After', i3, ' iterations, x = ', f15.10)
            endif
        enddo
        if (k > maxiteration) then
            ! might not have converged
            fx = f(x)
            if (abs(fx) > epsi) then
                write(*,*) "*** Warning: has not converged"
            endif
        endif 
        ! number of iterations taken:
        iters = k-1
    end subroutine Newton
end module newtonmod
!______________________________________________________________________ 
program rootfinder
    use dp
    use vars
    use newtonmod
    use functions
    implicit none
    integer :: maxiteration
    real(wp) :: x, x0, fx
    real(wp) :: epsi
    integer :: iters
    integer :: exitflag 
    !-------------------------------------
    11  format('x = ',f15.10,' after',i3,' iterations')
    12  format('f(x) = ', f15.10)
    
    write(*,*) "Computing root of function using Newton's method."
    exitflag = 0
    write(*,*) "Enter your initial guess."
    read(*,*) x0
    write(*,*) "Enter the maximum iterations allowed."
    read(*,*) maxiteration
    Write(*,*) "Enter the precision required to converge."
    read(*,*) epsi
    write(*,*)
    
    call Newton(f_sqrt, fprime_sqrt, x0, x, maxiteration,iters,exitflag,epsi)
    
    write(*,11) x, iters
    fx = f_sqrt(x) !store y value returned from function
    write(*,12) fx !write the y value at the root, aproxx zero
end program rootfinder
!______________________________________________________________________
module functions
    contains
    !begin the problem specific function 
    real(wp) function f_sqrt(x)
        use dp
        use vars
        implicit none
        real(wp), intent(in) :: x
        !function here
        f_sqrt = (c1*x**5) + (c2*x**2) + (c3*x) + (c4)
    end function f_sqrt
    !derivative of problem specific function
    real(wp) function fprime_sqrt(x)
        use dp
        use vars
        implicit none
        real(wp), intent(in) :: x
        !derivative here
        fprime_sqrt = (5*c1*x**4) + (2*c2*x) + (c3)
    end function fprime_sqrt
end module functions
