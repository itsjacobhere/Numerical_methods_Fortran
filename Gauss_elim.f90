!Written by Jacob Turner
!finished on 4-9-18
!this program will perform guassian elimination
!in a subroutine on a matrix of specified length
!n & ia are the dimensions of the matrix
!the dimensions and matrix will be read in 
!from an external file.
!_________________________________________________________________
program guass_elim
    !declarations
    !integer :: n, ia
    integer, parameter :: n=3, ia=3
    real, dimension (:,:), allocatable:: a 
    real, dimension (:), allocatable:: s,b,x
    integer, dimension(:), allocatable:: l
    integer :: ierror
    character(len=20) :: filename="matrix.dat"
    !---------------------------------------------------
    interface
        !reference to Subroutine that completes guass elimination
        subroutine gauss(n,a,ia,l,s)
            integer, intent(in) :: n,ia
            real, dimension (:,:), intent(inout)::a
            real, dimension (n)::s
            integer, dimension(n) , intent(out)::l
        end subroutine gauss
        subroutine solver(n,a,ia,l,b,x)
            integer, intent(in)::n,ia
            real, dimension(:,:), intent(in):: a
            real, dimension(:), intent(in) :: b
            integer, dimension(n), intent(in):: l
            real, dimension(n), intent(out) :: x
        end subroutine solver
    end interface
    !-----------------------------------------------------
    !uncomment for user input matrix dimensions
    !write(*,*) "What are the dimensions of the matrix?"
    !read(*,*) n, ia
    !initialize and allocate
    allocate(a(n,ia))
    allocate(s(n))
    allocate(b(n))
    allocate(x(n))
    allocate(l(n))
    a = 0
    b = 0

    file: do
        !    write(*,*) "What file would you like to open?"
        !    write(*,*)
        !    read(*,*) filename
        open(unit=10,file=filename,status="old",iostat=ierror)
        if (ierror/=0) then
            write(*,*) "File entered does not exist: ", filename
        else if (ierror==0) then
            write(*,*)
            write(*,*) "File, ",trim(filename)," was opened succesfully."
            write(*,*)
            exit
        end if 
    end do file
    !read a matrix from file
    read(10,*) ((a(i,j), j=1,n),i=1,ia)
    !read b vector
    read(10,*) (b(i), i = 1,n)
    !_________________________________________________________________
    write(*,*) "Matrix A:"
    b(1)= 1.0; b(2)= 4.0; b(3)= -1.0!; b(4)=26.0
    do i=1,n
        write(*,*) (a(i,j),j=1,n)
    end do
    write(*,*) "Vector b:"
    write(*,*) (b(i), i=1,n)
    !use guass subroutine
    call  gauss(n,a,n,l,s)
    
    write(*,*) "Matrix A after guass:"
    do i=1,n
    write(*,*) (a (i,j), j=1,n)
    end do
    !use solver subroutine
    call solver(n,a,ia,l,b,x)
    
    write(*,*) "x(i) vector is: " 
    write(*,*)
    write(*,*) (x(i), i=1,n)
    stop
end program guass_elim
!_________________________________________________________________
subroutine gauss(n,a,ia,l,s)    
    integer, intent(in) :: n, ia
    real, dimension (:,:), intent(inout):: a
    real, dimension (n):: s
    integer, dimension(n), intent(out)::l
    !-----------------------------------------
do  i = 1,n
    l(i) = i  
    smax = 0.0
    do j = 1,n
        smax = amax1(smax,abs(a(i,j)))
    end do
    s(i) = smax 
end do 

do  k = 1,n-1
    rmax = -1.0
    do i = k,n
        r = abs(a(l(i),k))/s(l(i))  
        if(r <= rmax)  exit    
        j = i   
        rmax = r
    end do 
    lk = l(j) 
    l(j) = l(k) 
    l(k) = lk  
    do i = k+1,n      
        xmult = a(l(i),k)/a(lk,k)   
        do j = k+1,n    
            a(l(i),j) = a(l(i),j) - xmult*a(lk,j) 
        end do 
        a(l(i),k) = xmult 
    end do 
end do   
end subroutine gauss 
!_________________________________________________________________
subroutine solver(n,a,ia,l,b,x)  
    integer, intent(in) :: n, ia !pass the dimensions
    real, dimension (:,:), intent(in):: a
    real, dimension (:), intent(inout) :: b
    real, dimension (n), intent(out):: x
    integer, dimension(n), intent(in)::l
    !local
    integer :: k, i, j
    !------------------------------------------------------
    do k = 1,n-1
        do i = k+1,n      
            b(l(i)) = b(l(i)) - a(l(i),k)*b(l(k)) 
        end do 
    end do 
    x(n) = b(l(n))/a(l(n),n)
    do i = n-1,1,-1     
        sum = b(l(i))       
        do j = i+1,n      
        sum = sum - a(l(i),j)*x(j)  
        end do 
        x(i) = sum/a(l(i),i)
    end do 
end subroutine solver
