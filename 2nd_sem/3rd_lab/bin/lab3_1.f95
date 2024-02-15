module Simpson_Runge
contains


real*16 Function f(x, num) result(res)
real*16 :: x
integer :: num

if (num==1) then
    res=x**5-5.2*x**3+2.5*x**2-7*x
elseif (num==2) then
	res=abs(x**5-5.2*x**3+2.5*x**2-7*x)
end if

end function 


real*16 Function Simpson(a,b,n,num) result(res)
implicit none
integer :: num,nn,i
real*16 ::a,b,n,h,a_i,f_a,f_b,f_i,sum_2,sum_4

f_a=f(a, num)
f_b=f(b, num)

nn=int(n)
h=(b-a)/n
sum_2=0
sum_4=0
do i=0,nn-1,1
	a_i=a+h*i+h/2.0
	f_i=f(a_i, num)
	sum_4=sum_4+f_i
end do
do i=1,nn-1,1
	a_i=a+h*i
	f_i=f(a_i, num)
	sum_2=sum_2+f_i
end do

res=(h/6.0)*(f_a+f_b+2.0*sum_2+4.0*sum_4)

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


character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)

end function

end module 


program lab2
use Simpson_Runge
implicit none
real*16 :: a,b,eps,integral,integral_1,n,integral_real_1,integral_real_2, ab(10), integral_real(10)
integer :: num,i,rung,v,q

integral_real(1)=27270901.0/600000.0
integral_real(2)=15813251.0/600000.0
integral_real(3)=1030301.0/75000.0
integral_real(4)=85184.0/9375.0
integral_real(5)=3442951.0/600000.0
integral_real(6)=83349.0/25000.0
integral_real(7)=1030301.0/600000.0
integral_real(8)=6859.0/9375.0
integral_real(9)=44217.0/200000.0
integral_real(10)=2197.0/75000.0
integral_real_1=52331868094707.0/2000000000000.0
integral_real_2=132578467905293.0/2000000000000.0

ab(1)=3.01
ab(2)=2.51
ab(3)=2.02
ab(4)=1.76
ab(5)=1.51
ab(6)=1.26
ab(7)=1.01
ab(8)=0.76
ab(9)=0.51
ab(10)=0.26
a=-2.6
b=2.3

open(21, file='bin/researches_res/1_err_on_eps.txt')
open(22, file='bin/researches_res/1_iter_on_eps.txt')
open(23, file='bin/researches_res/1_err_on_eps_abs.txt')
open(24, file='bin/researches_res/1_iter_on_eps_abs.txt')
open(25, file='bin/researches_res/1_err_on_len.txt')

write(*,*)' 1st+2nd+4th researches'
eps=1.0
do v=1,12

	write(*,*)v,'/12'
	eps=eps/10
	q=0 ! 1st iteration
	n=1.0
	integral_1=Simpson(a,b,n,num)
	
	do num=1,2

		do

			q=q+1 ! next iteration
			n=n*2.0 
			integral=Simpson(a,b,n,num)
			rung=Runge(integral_1,integral,eps)
			if (rung==1) then
				if (num==1) then
					write(21,*) abs(integral-integral_real_1),eps
					write(22,*) q,eps
				elseif (num==2) then 
					write(23,*) abs(integral-integral_real_2),eps
					write(24,*) q,eps
				end if
				exit
			else 
				integral_1=integral
			end if

		end do

	end do

end do

write(*,*)' 3rd research'
num=1
eps=0.000000000001
do v=1,10

	write(*,*)v,'/10'
	a=-ab(v)
	b=ab(v)
	q=0 ! 1st iteration
	n=1.0
	integral_1=Simpson(a,b,n,num)
	do

		q=q+1 ! next iteration
		n=n*2.0 
		integral=Simpson(a,b,n,num)
		rung=Runge(integral_1,integral,eps)
		if (rung==1) then
			write(25,*) abs(integral-integral_real(v)),b-a
			exit
		else 
			integral_1=integral
		end if

	end do

end do

write(*,*)' end program'

end program