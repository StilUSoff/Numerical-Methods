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
res=1+((ranperc+(2*percent-1)/2.0)*0.01)*(-1)**j
end function ranN

real*16 Function f(x, num, ranfluct) result(res)
real*16 :: x
integer :: num, ranfluct

if (num==1) then
        res=exp(-x**2)+1-x
elseif (num==2) then
		if (ranfluct==0) then
		    res=x**4-6*x**2+12*x-8
		else
			res=ranN(ranfluct)*x**4-6*x**2+12*x-8
        end if
elseif (num==3) then
		res=sin(x)-exp(-1/x)
else
		print*, "nope"
end if
end function 

end module 


module bisection_and_secant_methods
use functions_and_rand
integer, public :: iter
real*16, public :: root, E, delt
contains

subroutine bis(a1,b1,E,n, num, ranfluct, iterfluctcheck)
implicit none
real*16 :: E, a1, a, b1, b, c, f1, f2, check
integer :: i, n, num, ranfluct, iterfluctcheck

a=a1
b=b1
if (iterfluctcheck==1) open(10, file='graph/iterflcut_bis.txt')
i=0
f1=f(b, num, ranfluct)
do while(i<n)
	c=(a+b)/2.0
	f2=f(c, num, ranfluct)
	check=f1*f2
	if (check>0) then 
		b=c
		f1=f2
	elseif (check<0) then
		a=c
	elseif (check==0) then
		a=c
		b=c		
	end if
	i=i+1
	if (iterfluctcheck==1) write(10, *) c, i
	delt=abs(a-b)
	if (delt<E) exit 
end do
if (iterfluctcheck==1) close(10)
if (iterfluctcheck==2) root=c
iter=i
write(*,*) "bisection: ", c, delt, i
end subroutine bis

	
subroutine sec(a1,b1,E,n,num, ranfluct, iterfluctcheck)
implicit none
real*16 :: E, a1, a, b1, b, c, f1, f2
integer :: i, n, num, ranfluct, iterfluctcheck
delt=0	
a=a1
b=b1
if (iterfluctcheck==1) open(11, file='graph/iterflcut_sec.txt')
i=0
f1=f(a, num, ranfluct)
f2=f(b, num, ranfluct)
do while(i<n)
    c=b-f2*((b-a)/(f2-f1))
    f1=f2
	f2=f(c, num, ranfluct)
    a=b
    b=c
    i=i+1
	if (iterfluctcheck==1) write(11, *) c, i
    delt=abs(a-b)
    if (delt<E) exit
end do
if (iterfluctcheck==1) close(11)
if (iterfluctcheck==2) root=c
iter=i
write(*,*) "secant: ", c, delt, i
end subroutine sec

end module


program lab1
use bisection_and_secant_methods
implicit none
real*16 :: a1, b1, xsec, xbis
integer :: n=1000, num, k, iterfluctcheck, ranfluct

write(*,"(A15,$)") 'what function? '
read(*,*) num
if (num<1 .or. num>3) stop
if (num==1) then
    a1=1
    b1=2
else if (num==2) then
    a1=1
	b1=2
else if (num==3) then
	a1=2
	b1=3
end if
ranfluct=0 !number 2
iterfluctcheck=1
E=1.0e-16
call bis(a1, b1, E, n, num, ranfluct, iterfluctcheck)
call sec(a1, b1, E, n, num, ranfluct, iterfluctcheck)
iterfluctcheck=2 !number 3+4
E=1.0e-1
open(12, file='graph/roots_bis.txt')
open(13, file='graph/roots_sec.txt')
do k=1, 15
	call bis(a1, b1, E, n, num, ranfluct, iterfluctcheck)
	write(12, *) root, iter, E
	if (k==15) xbis=root
	call sec(a1, b1, E, n, num, ranfluct, iterfluctcheck)
	write(13, *) root, iter, E
	if (k==15) xsec=root
	E=E/10.0
end do
close(12)
close(13)
if (num==2) then !number 5
iterfluctcheck=2
E=1.0e-16
open(14, file='graph/randfluct_bis.txt')
open(15, file='graph/randfluct_sec.txt')
do ranfluct=1, 5
	do k=1,10
		call bis(a1, b1, E, n, num, ranfluct, iterfluctcheck)
		write(14, *) ranfluct, abs(xbis-root)/abs(root)
		call sec(a1, b1, E, n, num, ranfluct, iterfluctcheck)
		write(15, *) ranfluct, abs(xsec-root)/abs(root)
	end do
end do
close(14)
close(15)
end if
end program