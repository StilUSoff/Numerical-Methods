module matrix_fluct
contains 


function matrix_fluct_B(X, L, iterfluctcheck) result(res)
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
res=1+((ranperc*percent)*0.01)*(-1)**ransign

end function

end module


module gradient_method_for_preconditioned_system
contains


subroutine precond(P,A,n,pup, L, U)
implicit none 
integer :: n, i, j, k, pup
real*16, dimension(n,n) :: P, A, L, U,  P_1, P_final, lev
real*16, dimension(n) :: y, itr
real*16  sum

do i=1,n
	do j=1,n
		U(i,j)=0
		L(i,j)=0
		P(i,j)=0
	end do
	L(i,i)=1
end do

do i=1,n
	do j=1,n 
		if (A(i,j)/=0) then
			lev(i,j)=0
		else
			lev(i,j)=100000
		end if
	end do
end do
do i=2,n
	do k=1,i-1
		if (lev(i,k) <= pup) then
			a(i,k)=a(i,k)/a(k,k)
			do j=k+1,n
				a(i,j)=a(i,j)-a(i,k)*a(k,j)
				lev(i,j)=min(lev(i,j),lev(i,k)+lev(k,j)+1)
			end do
		end if
	end do
	do j=2,n 
		if (lev(i,j) > pup) then
			a(i,j)=0
		end if
	end do
end do

do i=1,n
	do j=1,n
		if (j>=i) then
			U(i,j)=a(i,j)
		else
			L(i,j)=a(i,j)
		end if
	end do
end do
do j=1,n
    do i=1,n
    	sum = 0.0  ! enables ACC parallelism for k-loop
        do k=1,n
        	sum = sum + L(i,k) * U(k,j)
        end do
        P(i,j) = sum
    end do
end do

end subroutine


function z_find(b,n,L,U) result(x)
implicit none	
real*16 :: b(n),L(n,n),U(n,n),x(n), y(n)
real*16 :: sum
integer :: n, i, p
		
y=0
x=0
do i=1,n
	sum=0
	do p=1,i-1 
		sum=sum+L(i,p)*y(p)
	end do
	y(i)=b(i)-sum 
end do 

x(n)=y(n)/U(n,n)
do i=1,n-1
	sum=0
	do p=0,i-1 
		sum=sum+U(n-i,n-p)*x(n-p)
	end do
	x(n-i)=(y(n-i)-sum)/U(n-i,n-i)
end do 

end function


subroutine gradient_method(A,B,X,eps,q,ERRORi,Xreal,n, kol, check_error,pup)
implicit none
integer :: k,i,j,n,q, kol, check_error, pup
real*16, allocatable, dimension(:) :: erro, zk, zk_1, rk, rk_1, pk, qk, A_X, A_pk, rk_1_P, rk_copy, rk_1_copy
real*16, allocatable, dimension(:,:) :: P, A_copy, P_copy, L, U
real*16, dimension(:,:) :: A, ERRORi
real*16, dimension(:) :: B,X,Xreal
real*16 :: eps, error, norm_r=0, norm_r0=0, ak=0, bk=0, rkT_zk, pkT_qk, zk_1T_rk_1, zkT_rk

allocate(zk(n))
allocate(zk_1(n))
allocate(rk(n))
allocate(rk_1(n))
allocate(pk(n))
allocate(qk(n))
allocate(A_X(n))
allocate(A_pk(n))
allocate(rk_1_P(n))
allocate(erro(n))
allocate(P(n,n))
allocate(A_copy(n,n))
allocate(L(n,n))
allocate(U(n,n))

A_copy=A
call precond(P,A_copy,n,pup,L,U)

do i=1,n 		!initial conditions
	A_X(i)=0
	do j=1,n
		A_X(i)=A_X(i)+A(i,j)*X(j)
	end do
	rk(i)=B(i)-A_X(i)
	norm_r0=norm_r0+rk(i)**2
end do
zk=z_find(rk,n,L,U)
pk=zk

q=0
do k=1,kol
	q=q+1
	norm_r=0
	error=0

	do i=1,n
		qk(i)=0
		do j=1,n
			qk(i)=qk(i)+A(i,j)*pk(j)
		end do
	end do
	rkT_zk=0
	pkT_qk=0
	do i=1,n
		rkT_zk=rkT_zk+rk(i)*zk(i)
		pkT_qk=pkT_qk+pk(i)*qk(i)
	end do
	ak=rkT_zk/pkT_qk
	do i=1,n
		X(i)=X(i)+ak*pk(i)
		rk_1(i)=rk(i)-ak*qk(i)
	end do
	zk_1=z_find(rk_1,n,L,U)
	zk_1T_rk_1=0
	zkT_rk=0
	do i=1,n
		zk_1T_rk_1=zk_1T_rk_1+zk_1(i)*rk_1(i)
		zkT_rk=zkT_rk+zk(i)*rk(i)
	end do
	bk=zk_1T_rk_1/zkT_rk
	do i=1,n
		pk(i)=zk_1(i)+bk*pk(i)
	end do

	do i=1,n	!saving for next iter
		rk(i)=rk_1(i)
		zk(i)=zk_1(i)
	end do

	if (check_error==1) then
		do i=1,n
			erro(i)=abs(Xreal(i)-X(i))
			error=error+erro(i)
		end do
		ERRORi(q,2)=error
		ERRORi(q,1)=q
	end if
	do i=1,n
		norm_r=norm_r+rk(i)**2
	end do
	if (((sqrt(norm_r))/(sqrt(norm_r0)))<eps) exit
end do
write(*,'(a3,$)')'i= '
write(*,*)q

deallocate(zk)
deallocate(zk_1)
deallocate(rk)
deallocate(rk_1)
deallocate(pk)
deallocate(qk)
deallocate(A_X)
deallocate(A_pk)
deallocate(rk_1_P)
deallocate(erro)

end subroutine


character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)

end function

end module 


program lab4
use matrix_fluct
use gradient_method_for_preconditioned_system
implicit none
real*16, allocatable, dimension (:,:) :: A, M, A_fluct,ERRORi, size_m
real*16, allocatable, dimension (:) :: B, X, B_fluct,erro,fluct,Xreal
integer ::  iterfluctcheck,i,j,n,k,v,q=0, kol=1000, number, pup, check_error
real*16 ran,error,time,start,finish,fluctuat,eps

eps=0.000000000001
allocate(size_m(49,2))
open(60, file='bin/matrix_res/matrix_size.txt')
do i=1,49
	read(60, *) size_m(i,1), size_m(i,2)
end do
close(60)

! 1st+2nd researches
open(61, file='bin/researches_res/error.txt')
open(62, file='bin/researches_res/time.txt')
write(*,*)'start error+time...'
pup=0
check_error=0
do v=5,53
	n=size_m(v-4,2)
	allocate(B(n))
	allocate(X(n))
	allocate(Xreal(n))
	allocate(erro(n))
	allocate(ERRORi(kol,2))

	do i=1,n
		X(i)=0 !initial approximation
		Xreal(i)=1 !real X
	end do
	number=size_m(v-4,1)
	write(*,'(a6,$)') 'Size: '
	write(*,'(i4)') number
    open(6+v, file='bin/matrix_res/A_matrix_'//trim(str(v))//'.txt')
    allocate(A(n,n))
	do i=1,n
        read(6+v,*)(A(i,j),j=1,n)
    end do
	do i=1,n
		B(i)=0
		do j=1,n
    		B(i)=B(i)+A(i,j)*Xreal(j)
    	end do
	end do
	call cpu_time(start)
	call gradient_method(A,B,X,eps,q,ERRORi,Xreal,n, kol,check_error,pup)
	call cpu_time(finish)
	time=finish-start
	error=0
	do i=1,n
		erro(i)=0
	end do
	do i=1,n
		erro(i)=abs(Xreal(i)-X(i))
		error=error+erro(i)
	end do
	write(61,*) v, error
	write(*,'(a6,$)') 'Time: ' 
	write(*,'(f14.7)') time
	write(*,*)
	write(62,*) v, time

	deallocate(A)
	deallocate(B)
	deallocate(X)
	deallocate(Xreal)
	deallocate(erro)
	deallocate(ERRORi)
end do
close(61)
close(62)
write(*,*)'end error+time...'

!3rd research
v=35
n=size_m(v-4,2)

allocate(A(n,n))
allocate(B(n))
allocate(fluct(n))
allocate(X(n))
allocate(Xreal(n))
allocate(B_fluct(n))
allocate(erro(n))
allocate(ERRORi(kol,2))

do i=1,n
	Xreal(i)=1 !real X
end do
write(*,*)'start fluct...'
open(63, file='bin/matrix_res/A_matrix_'//trim(str(v))//'.txt')
    do i=1,n
        read(63,*)(A(i,j),j=1,n)
    end do
close(63)
open(65, file='bin/researches_res/fluct.txt')
do iterfluctcheck=1,3
    do i=1,n
        X(i)=0
    end do
    do i=1,n
		B(i)=0
        do j=1,n
            B(i)=B(i)+A(i,j)*Xreal(j)
        end do
    end do
    write(*,'(a13,$)') 'Fluctuation: '
    write(*,'(i1)') iterfluctcheck
    B_fluct=matrix_fluct_B(B, n, iterfluctcheck)
    call gradient_method(A,B_fluct,X,eps,q,ERRORi,Xreal,n,kol,check_error,pup)
    fluctuat=0
    do i=1,n
		fluct(i)=0
        fluct(i)=100*(abs(Xreal(i)-X(i)))/abs(Xreal(i))
        fluctuat=fluctuat+fluct(i)
    end do
    write(65,*) iterfluctcheck, fluctuat
	write(*,*)fluctuat
    write(*,*)

end do
close(65)
write(*,*)'end fluct...'
deallocate(A)
deallocate(B)
deallocate(fluct)
deallocate(X)
deallocate(Xreal)
deallocate(B_fluct)
deallocate(erro)
deallocate(ERRORi)

! 4th + 5th researches

allocate(A(n,n))
allocate(B(n))
allocate(X(n))
allocate(Xreal(n))
allocate(erro(n))
allocate(ERRORi(kol,2))

do i=1,n
	Xreal(i)=1 !real X
end do
open(63, file='bin/matrix_res/A_matrix_'//trim(str(v))//'.txt')
    do i=1,n
        read(63,*)(A(i,j),j=1,n)
    end do
close(63)
open(66, file='bin/researches_res/error_eps.txt')
open(67, file='bin/researches_res/iter_eps.txt')
write(*,*)'start error+iter on eps...'
eps=1.0
do v=1,16
    eps=eps/10
    do i=1,n
        X(i)=0
    end do
    do i=1,n
		B(i)=0
        do j=1,n
            B(i)=B(i)+A(i,j)*Xreal(j)
        end do
    end do
    write(*,'(a9,$)') 'Epsilon: '
    write(*,'(F18.16)') eps
    call gradient_method(A,B,X,eps,q,ERRORi,Xreal,n,kol,check_error,pup)
    error=0
    do i=1,n
        erro(i)=0
    end do
    do i=1,n
        erro(i)=abs(Xreal(i)-X(i))
        error=error+erro(i)
    end do
    write(66,*) eps, error
    write(67,*) eps, q
    write(*,*)
end do
close(66)
close(67)
write(*,*)'end error+iter on eps...'
deallocate(A)
deallocate(B)
deallocate(X)
deallocate(Xreal)
deallocate(erro)
deallocate(ERRORi)

! 6th researches
check_error=1
v=35
n=size_m(v-4,2)

allocate(A(n,n))
allocate(B(n))
allocate(X(n))
allocate(Xreal(n))
allocate(erro(n))
allocate(ERRORi(kol,2))

do i=1,n
	Xreal(i)=1 !real X
end do
open(63, file='bin/matrix_res/A_matrix_'//trim(str(v))//'.txt')
    do i=1,n
        read(63,*)(A(i,j),j=1,n)
    end do
close(63)
open(68, file='bin/researches_res/error_iter.txt')
write(*,*)'start error on iter...'
eps=0.000000000001
do i=1,n
    X(i)=0
end do
do i=1,n
	B(i)=0
    do j=1,n
        B(i)=B(i)+A(i,j)*Xreal(j)
    end do
end do
call gradient_method(A,B,X,eps,q,ERRORi,Xreal,n,kol,check_error,pup)
do i=1,q
    write(68,*) ERRORi(i,1), ERRORi(i,2)
end do
close(68)
write(*,*)'end error on iter...'
deallocate(A)
deallocate(B)
deallocate(X)
deallocate(Xreal)
deallocate(erro)
deallocate(ERRORi)

! 7th researches
v=35
n=size_m(v-4,2)
check_error=0
allocate(A(n,n))
allocate(B(n))
allocate(X(n))
allocate(Xreal(n))
allocate(erro(n))
allocate(ERRORi(kol,2))

eps=0.00000000000001
do i=1,n
	Xreal(i)=1 !real X
end do
open(63, file='bin/matrix_res/A_matrix_'//trim(str(v))//'.txt')
    do i=1,n
        read(63,*)(A(i,j),j=1,n)
    end do
close(63)
open(68, file='bin/researches_res/error_p.txt')
write(*,*)'start error on p...'
do pup=0,15,3
    do i=1,n
        X(i)=0
    end do
    do i=1,n
		B(i)=0
        do j=1,n
            B(i)=B(i)+A(i,j)*Xreal(j)
        end do
    end do
    write(*,'(a3,$)') 'p= '
    write(*,'(i2)') pup
    call gradient_method(A,B,X,eps,q,ERRORi,Xreal,n,kol,check_error,pup)
    error=0
    do i=1,n
        erro(i)=0
    end do
    do i=1,n
        erro(i)=abs(Xreal(i)-X(i))
        error=error+erro(i)
    end do
    write(68,*) pup, error
    write(*,*) error
	write(*,*)
end do
close(68)
write(*,*)'end error on p...'
deallocate(A)
deallocate(B)
deallocate(X)
deallocate(Xreal)
deallocate(erro)
deallocate(ERRORi)

end program