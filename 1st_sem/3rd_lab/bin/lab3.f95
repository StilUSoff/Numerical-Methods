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


module ortogonal_method
contains

subroutine ort_meth(A,B)
implicit none
integer :: k,z,i,j,n,f,g,h
real*16, allocatable, dimension(:,:) :: M,D,M1
real*16, allocatable, dimension(:) :: C,S,E,R,E1
real*16, dimension(:,:) :: A !на вход. А
real*16, dimension(:) :: B !на вход. В

n=size(B)
allocate(M(n,n))
allocate(M1(n,n))
allocate(D(n,n))
allocate(R(n))
allocate(E(n))
allocate(E1(n))
do i=1,n
	do j=1,n
	    D(i,j)=A(i,j)
	end do
end do
do i=1,n
	R(i)=B(i)
end do
do i=1,n
	do j=1,n
	    M(i,j)=1
	end do
end do
do i=1,n
    E(i)=1
end do
allocate(S(n))
allocate(C(n))
do k=n,2,-1
	z=n+1-k
	f=n-k
	do i=n,1,-1
		if (n-k==i) exit
		do j=n,1,-1
			if (n-k==j) exit
			M(i-f,j-f)=D(i,j)
		end do
	end do
	do i=n,1,-1
		if (n-k==i) exit
		E(i-f)=R(i)
	end do
	do i=1,k-1
		do g=1,n
			E1(g)=E(g)
			do h=1,n
				M1(g,h)=M(g,h)
			end do
		end do
		C(i)=(M1(1,1))/(sqrt((M1(1,1))**2+(M(i+1,1))**2))
		S(i)=(M(i+1,1))/(sqrt((M1(1,1))**2+(M(i+1,1))**2))
		do j=1,k
			M(1,j)=C(i)*M1(1,j)+S(i)*M(i+1,j)
			M(i+1,j)=-S(i)*M1(1,j)+C(i)*M(i+1,j)
		end do
		E(1)=C(i)*E1(1)+S(i)*E(i+1)
		E(i+1)=-S(i)*E1(1)+C(i)*E(i+1)
	end do
	do i=n,1,-1
		if (n-k==i) exit
		do j=n,1,-1
			if (n-k==j) exit
			D(i,j)=M(i-f,j-f)
		end do
	end do
	do i=n,1,-1
		if (n-k==i) exit
	    R(i)=E(i-f)
	end do
end do
do i=1,n
	do j=1,n
		if (abs(D(i,j))<0.0000000000000001) then
			D(i,j)=0.0000000000000000
		end if
		A(i,j)=D(i,j)
	end do
end do
do i=1,n
    B(i)=R(i)
end do

end subroutine


function find_X(Z,T) result(X)
real*16, dimension (:,:) :: Z
real*16, dimension (:) :: T
real*16, allocatable, dimension (:) :: X,Xmeth
integer :: n,i,j

n=size(T)
allocate(Xmeth(n))
allocate(X(n))

do i=n,1,-1
	Xmeth(i)=0
	X(i)=0
end do
Xmeth(n)=T(n)
X(n)=T(n)/Z(n,n)
do i=n-1,1,-1
	do j=n,1,-1
		if (i==j)  exit
		Xmeth(i)=Xmeth(i)-Z(i,j)*X(j)
	end do
	Xmeth(i)=Xmeth(i)+T(i)
	X(i)=Xmeth(i)/Z(i,i)
end do

end function

character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)

end function

end module 


program lab3
use matrix_fluct
use ortogonal_method
implicit none
real*16, allocatable, dimension (:,:) :: A, A_fluct,Z
real*16, allocatable, dimension (:) :: B, X, B_fluct,Xmeth,X0,T,erro,fluct
integer :: L, iterfluctcheck,i,j,k, ob,n,v
real*16 ran, error, time, start, finish,fluctuat
character :: form1, form2

ob=1
L=10
allocate(A(L,L))
allocate(B(L))
allocate(fluct(L))
allocate(Z(L,L))
allocate(T(L))
allocate(X(L))
allocate(Xmeth(L))
allocate(X0(L))
allocate(A_fluct(L,L))
allocate(B_fluct(L))
allocate(erro(L))
n=size(B)

do i=1,n
	X0(i)=i
end do

! 1st+2nd researches
open(17, file='bin/researches_res/error.txt')
open(18, file='bin/researches_res/time.txt')
do v=1,6
	call cpu_time(start)
	write(*,'(a18,$)') 'Condition number: '
	write(*,'(i1)') v
    open(10+v, file='bin/matrix_res/'//trim(str(v))//'matrix.txt')
	do i=1,n
		read(10+v, *) A(i,1), A(i,2), A(i,3), A(i,4), A(i,5), A(i,6), A(i,7), A(i,8), A(i,9), A(i,10)
	end do
	do i=1,n
		B(i)=0
	end do
	do i=1,n
		do j=1,n
    		B(i)=B(i)+A(i,j)*X0(j)
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
	call ort_meth(Z,T)
	X=find_X(Z,T)
	error=0
	do i=1,n
		erro(i)=0
	end do
	do i=1,n
		erro(i)=abs(X0(i)-X(i))
		error=error+erro(i)
	end do
	write(17,*) 10**(1-v), error
	close(10+v)
	call cpu_time(finish)
	time=finish-start
	write(*,'(a6,$)') 'Time: ' 
	write(*,'(f9.7)') time
	write(*,*)
	write(18,*) 10**(1-v), time
end do
close(17)
close(18)

! 3rd research
open(19, file='bin/researches_res/Bfluct.txt')
do iterfluctcheck=1,3
	write(*,'(a13,$)') 'Fluctuation: '
	write(*,'(i1)') iterfluctcheck
    open(11, file='bin/matrix_res/1matrix.txt')
	do i=1,n
		read(11, *) A(i,1), A(i,2), A(i,3), A(i,4), A(i,5), A(i,6), A(i,7), A(i,8), A(i,9), A(i,10)
	end do
	do i=1,n
		B(i)=0
	end do
	do i=1,n
		do j=1,n
    		B(i)=B(i)+A(i,j)*X0(j)
    	end do
	end do
	B_fluct=matrix_fluctB(B,L,iterfluctcheck)
	do i=1,n
		do j=1,n
		    Z(i,j)=A(i,j)
		end do
	end do
	do i=1,n
		T(i)=B_fluct(i)
	end do
	call ort_meth(Z,T)
	X=find_X(Z,T)
	fluctuat=0
	do i=1,n
		fluct(i)=0
	end do
	do i=1,n
		fluct(i)=100*(abs(X0(i)-X(i)))/abs(X0(i))
		fluctuat=fluctuat+fluct(i)
	end do
	write(19,*) iterfluctcheck, fluctuat
	close(11)
	write(*,*)
end do
close(19)


open(20, file='bin/researches_res/Afluct.txt')
do iterfluctcheck=1,3
	call cpu_time(start)
	write(*,'(a13,$)') 'Fluctuation: '
	write(*,'(i1)') iterfluctcheck
    open(21, file='bin/matrix_res/1matrix.txt')
	do i=1,n
		read(21, *) A(i,1), A(i,2), A(i,3), A(i,4), A(i,5), A(i,6), A(i,7), A(i,8), A(i,9), A(i,10)
	end do
	do i=1,n
		B(i)=0
	end do
	do i=1,n
		do j=1,n
    		B(i)=B(i)+A(i,j)*X0(j)
    	enddo
	end do
	A_fluct=matrix_fluctA(A, L, iterfluctcheck)
	do i=1,n
		do j=1,n
		    Z(i,j)=A_fluct(i,j)
		end do
	end do
	do i=1,n
		T(i)=B(i)
	end do
	call ort_meth(Z,T)
	X=find_X(Z,T)
	fluctuat=0
	do i=1,n
		fluct(i)=0
	end do
	do i=1,n
		fluct(i)=100*(abs(X0(i)-X(i)))/abs(X0(i))
		fluctuat=fluctuat+fluct(i)
	end do
	write(20,*) iterfluctcheck, fluctuat
	close(21)
	write(*,*)
end do
close(20)

end program