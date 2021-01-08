!______________________________________________________________________
!written by Jacob Turner
!This program wil perform secant root finding method
!finished on 3/20/2018
!______________________________________________________________________dp
module dp
    implicit none
    integer, parameter :: wp = selected_real_kind(15)
end module dp
!______________________________________________________________________vars
module vars
    use dp
    implicit none
    integer,parameter :: h = 10, L = 95,hp = 100
    real(wp),parameter :: v = 0.000001007,rho = 998.2, e = 0.6, g = 9.81, pi = 3.14,Q = 0.3, r = 0.0002591
    real(wp) :: A, Re, p, hf, Ep, ve
end module vars

!______________________________________________________________________newton
module roots
    implicit none
    contains         
    subroutine secantm(f,x1,x2,x,maxiteration,iters,exitflag,epsi)
    use dp
    use vars
		!use functions
    ! Estimate the zero of f(x) using Newton's method. 
    ! Input:
    !   f:  the function to find a root of
    !   x1 & x2: the initial guesses
    !   the estimate x satisfying f(x)=0 (assumes Newton converged) 
    !   the number of iterations iters
        implicit none
        real(wp), intent(inout) :: x1,x2
        !real(wp), external :: f
        integer, intent(in) :: exitflag
        real(wp), intent(out) :: x
        real(wp), intent(in) :: epsi
        integer, intent(out) :: iters
        integer, intent(inout) :: maxiteration
				
        interface
            function f(d)
            use dp
            use vars
                implicit none
                real(wp), intent(in) :: d
            real(wp) :: f
            end function f
        end interface
				
        ! Declare any local variables:
        real(wp) :: y1,y2
        integer :: k
        if (exitflag == 0) then
            write(*,'(A,2x,f5.2,A,f5.2)') "Initial guesses:", x1,", ", x2
        endif
        ! find a zero of f(x) 
        do k=1,maxiteration
            ! evaluate function and its derivative:
            y1 = f(x1)

            y2 = f(x2)

            if (abs(y2) < epsi) then !slow progress criteria?
                exit  ! jump out of do loop
            endif
            ! increment x:
            x = (x1*y2-x2*y1)/(y2-y1)
            if (exitflag == 0) then
                write(*,12) k,x
                12 format('After', i3, ' iterations, x = ', f15.10)
            endif
            x1 = x2
            x2 = x
        enddo
        if (k > maxiteration) then
            ! might not have converged
            y2 = f(x)
            if (abs(y2) > epsi) then
                write(*,*) "*** Diverging ***"
            endif
        endif 
        ! number of iterations taken:
        iters = k-1
    end subroutine secantm
end module roots
!______________________________________________________________________f
!module functions
!contains
    function f(d)
        use dp
        use vars
        implicit none
        real(wp), intent(in) :: d
        real(wp) :: f

                Ep = (76.04*e*hp)/(rho*Q)
                hf = Ep - h
                A = pi*((d**2)/4)
                ve = Q/A
                p = (2*d*g*hf)/(L*v**2)
                Re = (ve*d)/v
				!a = r/(3.7*d)
				!b = 2.51/(Re*sqrt(p))
				
                
                !f = 2*LOG10((r/(3.7*d))+(2.51/(Re*sqrt(p))))+(1/(sqrt(p)))
                !f = (d+2)*(d+2)
                !f = d**4 - 81
                
                f = -1/sqrt(p) - (2*log10(r/(3.7*(d)))) + (2.51/(Re*sqrt(p)))
                
                write(*,*) hp
                write(*,*) v
                write(*,*) g
                write(*,*) h
                write(*,*) r
                write(*,*) Q
            write(*,*) pi
                write(*,*) ve
                write(*,*) hf
                write(*,*) Re
                write(*,*) Ep
                write(*,*) rho
                write(*,*) e
                write(*,*) p
                !read(*,*) 
    end function f
!end module functions
!______________________________________________________________________ program:
program rootfinder
use dp
use vars
use roots
!use functions
implicit none
integer :: maxiteration, iters
real(wp) :: x, fx, x1, x2
real(wp) :: epsi
integer :: exitflag 

interface
    function f(d)
        use dp
        use vars
        implicit none
        real(wp), intent(in) :: d
        real(wp) :: f
    end function f
end interface

Re=0
hf=0
Ep=0
ve=0 
p=0
write(*,*) "Computing root of the function using Secant method."
exitflag = 0
write(*,*) "Enter your initial guesses."
read(*,*) x1,x2
!read(*,*) x2
write(*,*) "Enter maximum iterations."
read(*,*) maxiteration
Write(*,*) "Enter the precision of the root estimate."
read(*,*) epsi
write(*,*) ' '  ! blank line
call secantm(f,x1,x2,x,maxiteration,iters,exitflag,epsi)
write(*,11) x, iters
11  format('x = ', f15.10, ' after', i3, ' iterations')
fx = f(x)
write(*,12) fx
12  format('f(x) = ', f15.10)
end program rootfinder
!______________________________________________________________________