module functions_and_rand
contains

real*16 Function ranN(percent) result(res)
integer :: percent, j
real :: ranperc,r

call random_number(ranperc)
call random_number(r)
if (r>=0.5) then
	j=1
else
	j=0
end if
res=1+(percent/100.0)*((ranperc+(2*percent-1)/2.0)*0.01)*(-1)**j

end function ranN


real*16 Function f(x, num, ranfluct) result(res)
real*16 :: x,randomchik
integer :: num, ranfluct

if (num==1) then
    res=(x**5-5.2*x**3+2.5*x**2-7*x)*sin(x**4)
elseif (num==2) then
	if (ranfluct==0) then
		res=1.0/(x**2+4.0*x+3.0)
	else 
		randomchik=ranN(ranfluct)
		res=randomchik/(x**2+4.0*x+3.0)
	end if
elseif (num==3) then
	res=log(1.0-x)/x
end if

end function 

end module


module Simpson_adapt
contains

real*16 Function Simpson(a,b,n,num) result(res)
use functions_and_rand
implicit none
integer :: num,nn,i,ranfluct
real*16 ::a,b,n,h,a_i,f_a,f_b,f_i,sum_2,sum_4

f_a=f(a, num,ranfluct)
f_b=f(b, num,ranfluct)

nn=int(n)
h=(b-a)/n
sum_2=0
sum_4=0
do i=0,nn-1,1
	a_i=a+h*i+h/2.0
	f_i=f(a_i, num,ranfluct)
	sum_4=sum_4+f_i
end do
do i=1,nn-1,1
	a_i=a+h*i
	f_i=f(a_i, num,ranfluct)
	sum_2=sum_2+f_i
end do

res=(h/6.0)*(f_a+f_b+2.0*sum_2+4.0*sum_4)

end function 


recursive subroutine a_quad_asr(num,a,b,eps,n,res,culc,ranfluct,iterations)
use functions_and_rand
implicit none
real*16 :: a,b,eps,n,res,integral_x,integral_y,resl,resr
integer :: num,culc,ranfluct,iterations

integral_x=Simpson(a,b,n,num)
integral_y=Simpson(a,b,2*n,num)

culc=culc+1
if (culc>1000000000) stop

if (abs(integral_x-integral_y)<eps) then
	res=integral_y
	if (iterations==8) write(23,*)a,(b-a)
else
	call a_quad_asr(num,a,(a+b)/2.0,eps,n,resl,culc,ranfluct,iterations)
	call a_quad_asr(num,(a+b)/2.0,b,eps,n,resr,culc,ranfluct,iterations)
	res=resl+resr
end if

end subroutine


subroutine quad_asr(res,a,b,eps,num,n,culc,ranfluct,iterations) 
use functions_and_rand
implicit none
real*16 :: a,b,m,int,fm,fa,fb,eps,res,n
integer :: num,culc,ranfluct,iterations

culc=0

call a_quad_asr(num,a,b,eps,n,res,culc,ranfluct,iterations)

end subroutine


character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)

end function

end module 


program lab2
use Simpson_adapt
use functions_and_rand
implicit none
real*16 :: a,b,eps,integral,integral_realll(3),integral_real, aa(3),bb(3), n, pi
integer :: num,i,rung,v,q,culc,ranfluct

pi=4*atan(1.0)
aa(1)=-1.51
bb(1)=1.51
aa(2)=1.0
bb(2)=100000.0
aa(3)=0.0000000000001
bb(3)=0.9999999
integral_realll(1)=1.064406881186503
integral_realll(2)=log(2.0)/2.0
integral_realll(3)=-(pi**2.0/6.0)

n=4.0

write(*,*)' 1st+2nd+3rd researches'

ranfluct=0

do num=1,3

	write(*,*)num,'/ 3'
	a=aa(num)
	b=bb(num)
	integral_real=integral_realll(num)
	open(21, file='bin/researches_res/2_'//trim(str(num))//'_err_on_eps.txt')
	open(22, file='bin/researches_res/2_'//trim(str(num))//'_iter_on_eps.txt')

	eps=1
	do i=1,16

		eps=eps/10.0
		if (i==8) open(23, file='bin/researches_res/2_'//trim(str(num))//'_len_on_x.txt')
		call quad_asr(integral,a,b,eps,num,n,culc,ranfluct,i)
		write(21,*) abs(integral-integral_real),eps
		write(22,*) culc,eps
		if (i==8) close(23)

	end do


	close(21)
	close(22)

end do


num=2
a=aa(num)
b=bb(num)
integral_real=integral_realll(num)


write(*,*)' 4th research'

do ranfluct=1,3

	write(*,*)ranfluct,'/ 3'
	open(21, file='bin/researches_res/2_'//trim(str(num))//'_err_fluct_'//trim(str(ranfluct))//'.txt')
	eps=1

	do i=1,16

		eps=eps/10.0
		if (i==8) open(23, file='bin/researches_res/2_'//trim(str(num))//'_len_fluct_'//trim(str(ranfluct))//'.txt')
		call quad_asr(integral,a,b,eps,num,n,culc,ranfluct,i)
		write(21,*) abs(integral-integral_real),eps
		if (i==8) close(23)

	end do

	close(21)

end do

write(*,*) 'end program'

end program