module RK_Runge
contains


real*16 Function func(x,y) result(res)
real*16 :: x,y
integer :: num

res=2*x*(x**2+y)

end function 


integer Function Runge(a,b,eps) result(res)
implicit none	
real*16 :: a,b,eps

if (((1.0/15.0)*abs(b-a))<eps) then
	res=1
else
	res=0
end if 

end function 


subroutine RK(x_max,y_0,x_0,h,eps,check)
implicit none
integer :: check, i
real*16 ::x_0,y_0,eps,h,h_max,x_max,k1,k2,k3,delt_y,x_old,x,y

x=x_0
y=y_0
h_max=eps
if (check==0) write(21,*) x, h
if (check==1) write(22,*) x, y
if (check==2) write(23,*) x, y 
i=1
do while (x<x_max)
	k1=func(x,y)*h
	k2=func(x+h/2.0,y+k1/2.0)*h
	k3=func(x+h,y-k1+2*k2)*h
	delt_y=(k1+4*k2+k3)/6.0
	x=x+h
	y=y+delt_y
    if ((h<=h_max)) then 
		x=x+h
		y=y+delt_y
		if (check==0) write(21,*) x, h
		if (check==1) write(22,*) x, y
		if (check==2) write(23,*) x, y 
	elseif ((delt_y<=eps)) then 
		x=x+h
		y=y+delt_y
		if (check==0) write(21,*) x, h
		if (check==1) write(22,*) x, y
		if (check==2) write(23,*) x, y 
	else 
		h=h*(eps/delt_y)**0.5
	end if 
end do
if ((check==10) .or. (check==0)) write(20,*) eps, y

end subroutine 


character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)

end function

end module 


program lab4
use RK_Runge
implicit none
real*16 :: x,x0,y,eps,h,h_start,y1,y2
integer :: check=10,v

open(20, file = 'bin/researches_res/3_eps_error.txt')
open(21, file = 'bin/researches_res/2_tcoord_step.txt')
open(22, file = 'bin/researches_res/1_t_y_0_1.txt')
open(23, file = 'bin/researches_res/1_t_y_0_05.txt')

write(*,*)' 2nd + 3rd researches'
check=10
eps=1.0
do v=1,6
	eps=eps/10
	print*, v, '/12'
	if (v==6) check=0
	x0=1.0
	h=0.1
	y=exp(1.0)
	x=2.0
	call RK(x, y, x0, h, eps, check)
	check=10
end do
close(20)
close(21)

write(*,*)' 1st researches'

eps=0.00001
h = 0.01
do v=1,2
	print*, v,'/2'
	x0=1.0
	y=exp(1.0)
	x=2.0
	call RK(x, y, x0, h, eps, v)
	h=0.000001
end do
close(22)
close(23)

write(*,*)' end program'

end program