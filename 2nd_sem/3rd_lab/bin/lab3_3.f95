module Ralston_method
contains


real*16 Function f(x, num) result(res)
real*16 :: x
integer :: num

if (num==1) then
    res=(x**5-5.2*x**3+2.5*x**2-7*x)*cos(2*x)
elseif (num==2) then
	res=abs((x**5-5.2*x**3+2.5*x**2-7*x)*cos(2*x))
end if

end function 


real*16 function ralston(a, b, n, num) result(result)
implicit none
real*16, intent(in) :: a, b
integer, intent(in) :: n, num
integer :: i
real*16, dimension(7) :: a_weights, t_knots
real*16 :: x_i, h, answ

if (n == 3) then
    a_weights = [0.555555555555556,0.888888888888889,0.555555555555556, 0.0,0.0,0.0,0.0]
    t_knots = [-0.774596669241483,0.0,0.774596669241483, 0.0,0.0,0.0,0.0]
else if (n == 4) then
    a_weights = [0.3478548451,0.6521451549,0.6521451549,0.3478548451, 0.0,0.0,0.0]
    t_knots = [-0.8611363116,-0.3399810436,0.3399810436,0.8611363116, 0.0,0.0,0.0]
else if (n == 5) then
    a_weights = [0.2369268851,0.4786286705,0.5688888889,0.4786286705,0.2369268851, 0.0,0.0]
    t_knots = [-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459, 0.0,0.0]
else if (n == 6) then
    a_weights=[0.171324492379170,0.360761573048139,0.467913934572691,0.467913934572691,0.360761573048139,0.171324492379170, 0.0]
    t_knots=[-0.932469514203152,-0.661209386466264,-0.238619186083197,0.238619186083197,0.661209386466264,0.932469514203152, 0.0]
else if (n == 7) then
    a_weights=[0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469,0.381830050505119,0.279705391489277,&
    0.129484966168870]
    t_knots=[-0.949107912342759,-0.741531185599394,-0.405845151377397,0.0,0.405845151377397,0.741531185599394,0.94910791234275]
else
    write(*,*) 'n error'
    return
end if
answ = 0
do i = 1, n
    x_i = (b+a)/2 + ((b-a)/2)*t_knots(i)
    answ = answ + a_weights(i) * f(x_i,num)
end do
result = answ * ((b-a)/2)

end function

 
recursive function calculate_quadrature(lower_limit, upper_limit, n, eps, q, i, num) result(result)
implicit none
real*16, intent(in) :: eps
integer :: n, q, i, num
real*16 :: lower_limit, upper_limit, result, result_left, result_right, integral, integral_double, error, middle

integral = ralston(lower_limit, upper_limit, n, num)
integral_double = ralston(lower_limit, upper_limit, n + 1, num)
error=abs(integral - integral_double) 
q=q+1
if (error <= 15*eps) then
    result = integral + error
    return
end if
middle = (upper_limit + lower_limit) / 2.0 
result_left = calculate_quadrature(lower_limit, middle, n, eps, q, i, num)
result_right = calculate_quadrature(middle, upper_limit, n, eps, q, i, num)
result = result_left + result_right

end function


character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)

end function

end module 


program lab2
use Ralston_method
implicit none
real*16 :: a,b,eps,integral,integral_check,integral_real(2), integral_richardson
integer :: num,i,q,n,p

integral_real(1)=(5.0*sin(2.0)+10.0*cos(2.0))/4.0 ! 0.096
integral_real(2)=3.17924893459238
a=-1.0
b=1.0
p=4

write(*,*)' 4th+6th+7th+9th researches'
do num=1,2
write(*,*); write(*,'(a5,$)') 'num: '; write(*,'(i1,$)') num; write(*,'(a2)') '/2' 
    
    do n=3,4
    write(*,*); write(*,'(a3,$)') 'n: '; write(*,'(i1,$)') n-2; write(*,'(a2)') '/2'
        
        open(19, file='bin/researches_res/3_'//trim(str(n))//'_err_on_eps_'//trim(str(num))//'.txt')
        open(20, file='bin/researches_res/3_'//trim(str(n))//'_iter_on_eps__'//trim(str(num))//'.txt')
        open(21, file='bin/researches_res/3_'//trim(str(n))//'_err_on_eps_richard.txt')
        eps=1
    
	    do i=1,12
            q=0
            eps=eps/10
            write(*,'(a5,$)') 'eps: '; write(*,'(i2,$)') i; write(*,'(a3)') '/12'
            integral = calculate_quadrature(a, b, n, eps, q, i, num)
            write(19, *) abs(integral-integral_real(num))/integral_real(num), eps, integral
            if (num==1) then 
                integral_richardson = integral + (integral - calculate_quadrature(a, b, n + 1, eps, q, i, num))/(2**p - 1)
                write(21,*) abs(integral_richardson-integral_real(num))/integral_real(num), eps, integral
            end if
           write(20,*) q, eps, n
        end do

        close(19)
        close(20)
        close(21)

	end do
    
end do

end program