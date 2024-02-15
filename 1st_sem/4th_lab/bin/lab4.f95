module matrix_fluct
contains 

function matrix_fluctA(X, L, iterfluctcheck) result(res)
implicit none 
real*16, allocatable, dimension (:,:) :: res
real*16, dimension (:,:) :: X
real*16 ran, plusminus, max
integer :: i,j,L, iterfluctcheck, iflag1, jflag1

allocate(res(L,L))
max=X(1,1)
do i = 1, L-1 ! length of sequence
    do j = 1, L-1
		if (max<X(i,j)) then
			max=X(i,j)
			iflag1=i
			jflag1=j
		end if
    end do
end do
X(iflag1,jflag1)=X(iflag1,jflag1)*ranN(iterfluctcheck)
do i=1,L
	do j=1,L
		res(i,j) = X(i,j)
	end do
end do

end function

function matrix_fluctB(X, L, iterfluctcheck) result(res)
implicit none 
real*16, allocatable, dimension (:) :: res
real*16, dimension (:) :: X
real*16 ran, plusminus
integer :: i,j, L, iterfluctcheck

allocate(res(L))
do i = 1, L ! length of sequence
    res(i) = X(i)*ranN(iterfluctcheck)
end do

end function

real*16 Function ranN(percent) result(res)
integer :: perc, ransign, percent
real*16 :: ranperc,ran

call random_number(ranperc)
call random_number(ran)
if (ran>=0.5) then
	ransign=1
else
	ransign=0
end if
res=1+((ranperc+(2*percent-1)/2.0)*0.01)*(-1)**ransign

end function

end module


module yakobi_iter_method
contains

subroutine yakobi_iter(A,B,X,eps,q,ERRORi,Xreal)
implicit none
integer :: k,q,i,j,n,f,g,h,p=100
real*16, allocatable, dimension(:) :: E,R,MAXi,Rsum,erro
real*16, dimension(:,:) :: A,ERRORi !на вход. А
real*16, dimension(:) :: B,X,Xreal !на вход. В
real*16 :: eps,max,error

n=size(B)
allocate(R(n))
allocate(Rsum(n))
allocate(E(n))
allocate(MAXi(n))
allocate(erro(n))
do i=1,n
	X(i)=(B(i)/A(i,i)) !матрица X(0)
end do
do i=1,n
	E(i)=X(i) !матрица Х за этот ход
end do
do i=1,n
    R(i)=X(i) !матрица Х за прошлый ход
end do

q=0
do k=1,p
	q=q+1
	do i=1,n
		Rsum(i)=0
	end do
	error=0
	do i=1,n
		do j=1,n
			Rsum(i)=Rsum(i)+A(i,j)*R(j)
		end do
	end do
	do i=1,n
		Rsum(i)=Rsum(i)-A(i,i)*R(i)
		E(i)=(B(i)-Rsum(i))/A(i,i)
		MAXi(i)=abs(E(i)-R(i))
	end do
	max=MAXi(1)
	do i = 1, n
		if (max<MAXi(i)) then
			max=MAXi(i)
		end if
    end do
	do i=1,n
    	R(i)=E(i)
		erro(i)=abs(Xreal(i)-R(i))
		error=error+erro(i)
	end do
	ERRORi(k,2)=error
	ERRORi(k,1)=q
	if (max<eps) exit
end do
write(*,'(a3,$)')'i= '
write(*,*)q
do i=1,n
    X(i)=E(i)
end do

end subroutine


character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)

end function

end module 


program lab4
use matrix_fluct
use yakobi_iter_method
implicit none
real*16, allocatable, dimension (:,:) :: A, A_fluct,Z,ERRORi
real*16, allocatable, dimension (:) :: B, X, B_fluct,T,erro,fluct,Xreal,cond
integer ::  iterfluctcheck,i,j,k,n,v,opsi,q=0
real*16 ran,error,time,start,finish,fluctuat,eps
character :: form1, form2


open(31, file='bin/matrix_res/size.txt')
read(31, *) n
close(31)
write(*,'(a16,$)') 'Size of matrix: ' 
write(*,'(i3)') n
allocate(cond(6))
allocate(A(n,n))
allocate(B(n))
allocate(fluct(n))
allocate(Z(n,n))
allocate(T(n))
allocate(X(n))
open(30, file='bin/matrix_res/cond.txt')
do i=1,6
	read(30, *) cond(i)
end do
close(30)
allocate(Xreal(n))
allocate(A_fluct(n,n))
allocate(B_fluct(n))
allocate(erro(n))
allocate(ERRORi(100,2))
do i=1,n !real X
	Xreal(i)=1
end do

eps=0.00000000000001
! 1st+2nd researches
open(17, file='bin/researches_res/error.txt')
open(18, file='bin/researches_res/time.txt')
write(*,*)'start error+time...'
do v=1,6
	do i=n,1,-1
		X(i)=0
	end do
	call cpu_time(start)
	write(*,'(a18,$)') 'Condition number: '
	write(*,'(i1)') v
    open(10+v, file='bin/matrix_res/'//trim(str(v))//'matrix.txt')
	do i=1,n
		do j=1,n
			read(10+v, *) A(i,j)
		end do
	end do
	do i=1,n
		B(i)=0
	end do
	do i=1,n
		do j=1,n
    			B(i)=B(i)+A(i,j)*Xreal(j)
    		end do
	end do
	do i=1,n
		do j=1,n
		    Z(i,j)=A(i,j)
		end do
	end do
	do i=1,n
		T(i)=B(i)
	end do
	call yakobi_iter(Z,T,X,eps,q,ERRORi,Xreal)
	error=0
	do i=1,n
		erro(i)=0
	end do
	do i=1,n
		erro(i)=abs(Xreal(i)-X(i))
		error=error+erro(i)
	end do
	opsi=cond(v)
	write(17,*) opsi, error
	close(10+v)
	call cpu_time(finish)
	time=finish-start
	write(*,'(a6,$)') 'Time: ' 
	write(*,'(f9.7)') time
	write(*,*)
	write(18,*) opsi, time
end do
close(17)
close(18)
write(*,*)'end error+time...'

! 3rd research
open(19, file='bin/researches_res/Bfluct.txt')
write(*,*)'start fluct for B...'
do iterfluctcheck=1,3
	do i=n,1,-1
		X(i)=0
	end do
	write(*,'(a13,$)') 'Fluctuation: '
	write(*,'(i1)') iterfluctcheck
    open(11, file='bin/matrix_res/6matrix.txt')
	do i=1,n
		do j=1,n
			read(11, *) A(i,j)
		end do
	end do
	do i=1,n
		B(i)=0
	end do
	do i=1,n
		do j=1,n
    		B(i)=B(i)+A(i,j)*Xreal(j)
    	end do
	end do
	B_fluct=matrix_fluctB(B,n,iterfluctcheck)
	do i=1,n
		do j=1,n
		    Z(i,j)=A(i,j)
		end do
	end do
	do i=1,n
		T(i)=B_fluct(i)
	end do
	call yakobi_iter(Z,T,X,eps,q,ERRORi,Xreal)
	fluctuat=0
	do i=1,n
		fluct(i)=0
	end do
	do i=1,n
		fluct(i)=100*(abs(Xreal(i)-X(i)))/abs(Xreal(i))
		fluctuat=fluctuat+fluct(i)
	end do
	write(19,*) iterfluctcheck, fluctuat
	close(11)
	write(*,*)
end do
close(19)
write(*,*)'end fluct for B...'

write(*,*)'start fluct for A...'
open(20, file='bin/researches_res/Afluct.txt')
do iterfluctcheck=1,3
	do i=n,1,-1
		X(i)=0
	end do
	write(*,'(a13,$)') 'Fluctuation: '
	write(*,'(i1)') iterfluctcheck
    open(21, file='bin/matrix_res/6matrix.txt')
	do i=1,n
		do j=1,n
			read(21, *) A(i,j)
		end do
	end do
	do i=1,n
		B(i)=0
	end do
	do i=1,n
		do j=1,n
    		B(i)=B(i)+A(i,j)*Xreal(j)
    	enddo
	end do
	A_fluct=matrix_fluctA(A, n, iterfluctcheck)
	do i=1,n
		do j=1,n
		    Z(i,j)=A_fluct(i,j)
		end do
	end do
	do i=1,n
		T(i)=B(i)
	end do
	call yakobi_iter(Z,T,X,eps,q,ERRORi,Xreal)
	fluctuat=0
	do i=1,n
		fluct(i)=0
	end do
	do i=1,n
		fluct(i)=100*(abs(Xreal(i)-X(i)))/abs(Xreal(i))
		fluctuat=fluctuat+fluct(i)
	end do
	write(20,*) iterfluctcheck, fluctuat
	close(21)
	write(*,*)
end do
close(20)
write(*,*)'end fluct for A...'


! 4th + 5th researches
open(22, file='bin/researches_res/error_eps.txt')
open(23, file='bin/researches_res/iter_eps.txt')
write(*,*)'start error+iter...'
eps=1.0
do v=1,16
	eps=eps/10
	do i=n,1,-1
		X(i)=0
	end do
	write(*,'(a9,$)') 'Epsilon: '
	write(*,'(F18.16)') eps
    open(21, file='bin/matrix_res/6matrix.txt')
	do i=1,n
		do j=1,n
			read(21, *) A(i,j)
		end do
	end do
	do i=1,n
		B(i)=0
	end do
	do i=1,n
		do j=1,n
    		B(i)=B(i)+A(i,j)*Xreal(j)
    	enddo
	end do
	do i=1,n
		do j=1,n
		    Z(i,j)=A(i,j)
		end do
	end do
	do i=1,n
		T(i)=B(i)
	end do
	call yakobi_iter(Z,T,X,eps,q,ERRORi,Xreal)
	error=0
	do i=1,n
		erro(i)=0
	end do
	do i=1,n
		erro(i)=abs(Xreal(i)-X(i))
		error=error+erro(i)
	end do
	write(22,*) eps, error
	write(23,*) eps, q
	close(21)
	write(*,*)
end do
close(22)
close(23)
write(*,*)'end error+iter...'


! 6th researches
eps=9.99999999999999999999999999999999858E-0017
open(24, file='bin/researches_res/error_iter.txt')
write(*,*)'start error on iter...'

do i=n,1,-1
	X(i)=0
end do
open(21, file='bin/matrix_res/1matrix.txt')
do i=1,n
	do j=1,n
		read(21, *) A(i,j)
	end do
end do
do i=1,n
	B(i)=0
end do
do i=1,n
	do j=1,n
		B(i)=B(i)+A(i,j)*Xreal(j)
   	enddo
end do
do i=1,n
	do j=1,n
	    Z(i,j)=A(i,j)
	end do
end do
do i=1,n
	T(i)=B(i)
end do
call yakobi_iter(Z,T,X,eps,q,ERRORi,Xreal)
do i=1,3
	write(24,*) ERRORi(i,1), ERRORi(i,2)
end do
close(21)
write(*,*)
close(24)
write(*,*)'start error on iter...'



end program