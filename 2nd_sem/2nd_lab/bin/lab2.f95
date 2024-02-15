module matrix_fluct
contains 


function matrix_fluct_A(X, L, iterfluctcheck) result(res)
implicit none 
real*16, allocatable, dimension (:,:) :: res
real*16, dimension (:,:) :: X
real*16 rand
integer :: i,j,L, iterfluctcheck

allocate(res(L,L))
do j=1,L
	X(j,6)=X(j,6)*ranN(iterfluctcheck)
end do
res = X

end function


real*16 Function ranN(percent) result(res)
integer :: perc, ransign, percent
real*16 :: ranperc,rand

call random_number(ranperc)
call random_number(ran)
if (rand>=0.5) then
	ransign=1
else
	ransign=0
end if
res=1+((ranperc*percent)*0.01)*(-1)**ransign

end function

end module


module ort_Givens__reverse_iter__LR
contains


subroutine LU(A,n,L,U)
implicit none 
integer :: n, i, j, k, pup
real*16, dimension(n,n) :: P, A, L, U,  P_1, P_final, lev
real*16, dimension(n) :: y, itr
real*16  sum

do i=1,n
	do j=1,n
		U(i,j)=0
		L(i,j)=0
	end do
end do
do i = 1,n
	do j = 1, n
		U(1,i) = A(1,i)
        L(i,1) = A(i,1)/U(1,1)
		sum = 0
        do k = 1,i
            sum = sum + L(i, k) * U(k, j)
		end do
        U(i, j) = A(i, j) - sum
        if (i > j) then
			L(j, i) = 0
		else
            sum = 0
            do k = 1,i
				sum = sum + L(j, k) * U(k, i)
			end do 
			L(j, i) = (A(j, i) - sum) / U(i, i)
        end if
	end do
	L(i,i)=1
end do

end subroutine


function y_find(b,n,L,U) result(x)
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


function norm(y,n) result(x)
implicit none	
real*16 :: y(n)
real*16 :: sum, x
integer :: n, i
		
sum=0
do i=1,n
	sum = sum+(y(i)**2)
end do 
x=sqrt(sum)

end function


subroutine reverse_iter(A_input,x,eps,lamb,n, q, kol) 
implicit none
integer :: k,i,j,n,q, kol, kk, p
real*16, allocatable, dimension(:) :: erro
real*16, dimension(:,:) :: A_input
real*16, dimension(:) :: x
real*16 :: L(n,n),U(n,n), A(n,n), x_1(n)
real*16 :: eps, lamb, lamb_1, x_1x_1, xx_1, start, x_norm

q=0
A=A_input
start=lamb
do i=1,n 
	A(i,i)=A(i,i)-lamb
end do
call LU(A,n,L,U)
do kk=1,kol
	lamb_1=lamb
	q=q+1
	x_1=x
	x=y_find(x,n,L,U)
	P=0
	xx_1=0
	x_1x_1=0
	do i=1,n 
		xx_1=xx_1+x(i)*x_1(i)
		x_1x_1=x_1x_1+x_1(i)*x_1(i)
	end do
	lamb=(xx_1)/(x_1x_1)
	x_norm=norm(x,n)
	x=x/x_norm
	if (abs(lamb-lamb_1)<eps) exit
end do
lamb=1/lamb+start

end subroutine


subroutine reverse_iter_shift(A_input,x,eps,lamb,n, q, kol) 
implicit none
integer :: k,i,j,n,q, kol, kk, p
real*16, allocatable, dimension(:) :: erro
real*16, dimension(:,:) :: A_input
real*16, dimension(:) :: x
real*16 :: L(n,n),U(n,n), A(n,n), y(n)
real*16 :: eps, lamb, lamb_1, y_norm, X_Y

q=0
do kk=1,kol
	q=q+1
	lamb_1=lamb
	A=A_input
	do i=1,n 
		A(i,i)=A(i,i)-lamb_1
	end do
	call LU(A,n,L,U)
	y=y_find(x,n,L,U)
	P=0
	X_Y=0
	do i=1,n
		if (abs(y(i))>eps) then
			P=P+1
			X_Y=X_Y+x(i)/y(i)
		end if
	end do
	lamb=lamb_1+X_Y/P
	y_norm=norm(y,n)
	x=y/y_norm
	if (abs(lamb-lamb_1)<eps) exit
end do

end subroutine


subroutine iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q,obs, iterfluctcheck)
use matrix_fluct
real*16, allocatable, dimension (:,:) :: A, A_fluct, L, U, VV
real*16, allocatable, dimension (:) :: X,erro,erro2, D, VV_1
real*16, dimension(:) :: obs
integer ::  research_chek, end_iter, jok, iterfluctcheck,i,j,n,k,v,lamb_iter,ii,shift_check, q, iter
real*16 ran,eps,lamb,VVx,x_len,VV_len, keker

allocate(A(n,n), L(n,n), U(n,n), VV(n,n),X(n), D(n), VV_1(n))

if (research_chek==52) then
	open(11, file='bin/matrix_res/A_matrix_sim_'//trim(str(jok))//'.txt')
	do i=1,n
		read(11,*)(A(i,j),j=1,n)
	end do
		open(12, file='bin/matrix_res/vect_sim_'//trim(str(jok))//'.txt')
	do i=1,n
		read(12,*)(VV(i,j),j=1,n)
	end do
	open(13, file='bin/matrix_res/numbers_sim_'//trim(str(jok))//'.txt')
	do i=1,1
		read(13,*)(D(j),j=1,n)
	end do
elseif (research_chek==51) then
	open(11, file='bin/matrix_res/A_matrix_unsim_'//trim(str(jok))//'.txt')
	do i=1,n
		read(11,*)(A(i,j),j=1,n)
	end do
		open(12, file='bin/matrix_res/vect_unsim_'//trim(str(jok))//'.txt')
	do i=1,n
		read(12,*)(VV(i,j),j=1,n)
	end do
	open(13, file='bin/matrix_res/numbers_unsim_'//trim(str(jok))//'.txt')
	do i=1,1
		read(13,*)(D(j),j=1,n)
	end do
else
	open(11, file='bin/matrix_res/A_matrix_'//trim(str(jok))//'.txt')
	do i=1,n
		read(11,*)(A(i,j),j=1,n)
	end do
	open(12, file='bin/matrix_res/vect_'//trim(str(jok))//'.txt')
	do i=1,n
		read(12,*)(VV(i,j),j=1,n)
	end do
	open(13, file='bin/matrix_res/numbers_'//trim(str(jok))//'.txt')
	do i=1,1
		read(13,*)(D(j),j=1,n)
	end do
end if
close(11)
close(12)
close(13)


kol=30000
erro=10
erro2=10

if (research_chek==2) then
 	A=matrix_fluct_A(A,n,iterfluctcheck)
end if

end_iter=int((maxval(D))*10000)
iter=int(((end_iter-10000)/140))

do lamb_iter=10000,end_iter, iter
	lamb=real(lamb_iter)/10000
	do i=1,n
		x(i)=sqrt(1.0/n) !initial approximation
	end do
	if (shift_check==0) then
		call reverse_iter(A,x,eps,lamb,n, q, kol)
	else
		call reverse_iter_shift(A,x,eps,lamb,n, q, kol)
	end if
	ii=1
	do i=1,n 
		if (abs(lamb-D(i))<(obs(jok)/20.0)) exit 
		ii= ii+1
	end do
	if (ii<=10) then
		if ((abs(D(ii)-lamb))<(erro(ii))) erro(ii)=abs(D(ii)-lamb)
		do j=1,n 
			VV_1(j)=VV(j,ii)
		end do
		x_len=norm(x,n)
		VV_len=norm(VV_1,n)
		VVx=0
		do j=1,n
			VVx=VVx+VV_1(j)*x(j)
		end do
		if ((abs(sin(acos((VVx)/(VV_len*x_len)))))<(erro2(ii))) erro2(ii)=abs(sin(acos((VVx)/(VV_len*x_len))))
	end if
end do
deallocate(A, L, U, VV,X, D, VV_1)

end subroutine


subroutine LR(A,eps,n, kol)
implicit none
integer :: k,i,j,n, kol, kk
real*16, allocatable, dimension(:) :: erro
real*16, dimension(:,:) :: A
real*16 :: L(n,n),U(n,n), P_1(n,n)
real*16 :: eps, check

do kk=1,kol
	check=A(1,1)
	call LU(A,n,L,U)
	do i=1,n
		do j=1,n
 			A(i,j)=0
 			do k=1,n
 				A(i,j) = A(i,j) + U(i,k)*L(k,j)
			end do
		end do
	end do 
	if ((abs(A(1,1)-check))<eps) exit
end do

end subroutine


subroutine LR_researches(research_chek, jok, erro, erro2, eps, n, q,obs, iterfluctcheck)
use matrix_fluct
real*16, allocatable, dimension (:,:) :: A, A_fluct, A_copy,A_copy_1, L, U, VV, VV_LR,P
real*16, allocatable, dimension (:) :: X,erro,erro2, D, VV_1, VV_LR_1, D_LR, b
real*16, dimension(:) :: obs
integer ::  research_chek, end_iter, jok, iterfluctcheck,i,j,n,k,v,lamb_iter,ii,shift_check, q, iter
real*16 ran,eps,lamb,VVx,x_len,VV_len,keker

allocate(A(n,n), L(n,n), U(n,n), VV(n,n), D(n), b(n), VV_1(n),X(n),A_copy_1(n,n),A_copy(n,n),D_LR(n),VV_LR_1(n),VV_LR(n,n), P(n,n))

if (research_chek==52) then
	open(11, file='bin/matrix_res/A_matrix_sim_'//trim(str(jok))//'.txt')
	do i=1,n
		read(11,*)(A(i,j),j=1,n)
	end do
		open(12, file='bin/matrix_res/vect_sim_'//trim(str(jok))//'.txt')
	do i=1,n
		read(12,*)(VV(i,j),j=1,n)
	end do
	open(13, file='bin/matrix_res/numbers_sim_'//trim(str(jok))//'.txt')
	do i=1,1
		read(13,*)(D(j),j=1,n)
	end do
elseif (research_chek==51) then
	open(11, file='bin/matrix_res/A_matrix_unsim_'//trim(str(jok))//'.txt')
	do i=1,n
		read(11,*)(A(i,j),j=1,n)
	end do
		open(12, file='bin/matrix_res/vect_unsim_'//trim(str(jok))//'.txt')
	do i=1,n
		read(12,*)(VV(i,j),j=1,n)
	end do
	open(13, file='bin/matrix_res/numbers_unsim_'//trim(str(jok))//'.txt')
	do i=1,1
		read(13,*)(D(j),j=1,n)
	end do
else
	open(11, file='bin/matrix_res/A_matrix_'//trim(str(jok))//'.txt')
	do i=1,n
		read(11,*)(A(i,j),j=1,n)
	end do
	open(12, file='bin/matrix_res/vect_'//trim(str(jok))//'.txt')
	do i=1,n
		read(12,*)(VV(i,j),j=1,n)
	end do
	open(13, file='bin/matrix_res/numbers_'//trim(str(jok))//'.txt')
	do i=1,1
		read(13,*)(D(j),j=1,n)
	end do
end if
close(11)
close(12)
close(13)

kol=10000

if (research_chek==2) then
 	A=matrix_fluct_A(A,n,iterfluctcheck)
end if

A_copy=A

call LR(A_copy,eps,n, kol)

do i=1,n

	ii=1
	D_LR(i)=A_copy(i,i)
	do j=1,n 
		if (abs(D(j)-D_LR(i))<(0.005)) exit
		ii=ii+1
	end do

	if (ii<=10) then

		erro(i)=abs(D(ii)-D_LR(i))
		A_copy_1=A
		do j=1,n
			A_copy_1(j,j)=A_copy_1(j,j)-D_LR(i)
		end do
		call LU(A_copy_1,n,L,U)
		b=eps
		VV_LR_1=y_find(b,n,L,U)
		do j=1,n 
			VV_1(j)=VV(j,ii)
		end do
		x_len=norm(VV_LR_1,n)
		VV_len=norm(VV_1,n)
		VVx=0
		do j=1,n
			VVx=VVx+VV_1(j)*VV_LR_1(j)
		end do
		erro2(i)=abs(sin(acos((VVx)/(VV_len*x_len))))

	end if

end do

deallocate(A, L, U, VV,X, D, VV_1)


end subroutine


character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)

end function

end module 


program lab2
use matrix_fluct
use ort_Givens__reverse_iter__LR
implicit none
real*16, allocatable, dimension (:,:) :: A, A_fluct, L, U, VV
real*16, allocatable, dimension (:) :: X,erro,erro2, D, VV_1, obs
integer ::  iterfluctcheck,i,j,n,k,v,q, kol=1000, number, pup, check_error, lamb_iter,ii, jok, end_iter,research_chek,shift_check
real*16 ran,eps,lamb,VVx,x_len,VV_len

n=10

write(*,*)'start reverse iter'
open(100, file='bin/matrix_res/obs.txt')
allocate(obs(10))
do i=1,1
	read(100,*)(obs(j),j=1,10)
end do
close(100)


write(*,*)'error on lamb...' ! 1st
research_chek=1

open(21, file='bin/researches_res/error_no_shift.txt')
open(22, file='bin/researches_res/error_shift.txt')
open(41, file='bin/researches_res/error_LR.txt')
eps=0.000000000001

do jok=1,8 ! 1st

	write(*,*)jok,'/ 8'

	shift_check=0
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(21,*) maxval(erro), maxval(erro2), obs(jok)
	deallocate(erro,erro2)

	shift_check=1
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(22,*) maxval(erro), maxval(erro2), obs(jok)
	deallocate(erro,erro2)

	allocate(erro(n),erro2(n))
	call LR_researches(research_chek, jok, erro, erro2, eps, n, q,obs, iterfluctcheck)
	write(41,*) maxval(erro), maxval(erro2), obs(jok)
	deallocate(erro,erro2)

end do
close(22)
close(21)
close(41)
write(*,*)'end error on lamb...'


write(*,*)'error on fluct...' ! 2nd
research_chek=2

open(23, file='bin/researches_res/fluct_no_shift.txt')
open(24, file='bin/researches_res/fluct_shift.txt')
open(42, file='bin/researches_res/fluct_LR.txt')
jok=5

do iterfluctcheck=1,5

	write(*,*)iterfluctcheck,'/ 5'

	shift_check=0
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(23,*) maxval(erro), maxval(erro2)
	deallocate(erro,erro2)
	
	shift_check=1
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(24,*) maxval(erro), maxval(erro2)
	deallocate(erro,erro2)

	allocate(erro(n),erro2(n))
	call LR_researches(research_chek, jok, erro, erro2, eps, n, q,obs, iterfluctcheck)
	write(42,*) maxval(erro), maxval(erro2)
	deallocate(erro,erro2)

end do
close(23)
close(24)
close(42)
write(*,*)'end error on fluct...' 


write(*,*)'work on simmetric...' ! 3th

open(31, file='bin/researches_res/error_sim_and_unsim_no_shift.txt')
open(32, file='bin/researches_res/error_sim_and_unsim_shift.txt')
open(43, file='bin/researches_res/error_sim_and_unsim_LR.txt')
eps=0.000000000001

do jok=1,2

	write(*,*)jok, ' / 2'

	shift_check=0
	research_chek=52
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(31,*) maxval(erro), maxval(erro2), obs(jok+8), research_chek-51
	deallocate(erro,erro2)
	research_chek=51
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(31,*) maxval(erro), maxval(erro2), obs(jok+8), research_chek-51
	deallocate(erro,erro2)

	shift_check=1
	research_chek=52
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(32,*) maxval(erro), maxval(erro2), obs(jok+8), research_chek-51
	deallocate(erro,erro2)
	research_chek=51
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(32,*) maxval(erro), maxval(erro2), obs(jok+8), research_chek-51
	deallocate(erro,erro2)

	research_chek=52
	allocate(erro(n),erro2(n))
	call LR_researches(research_chek, jok, erro, erro2, eps, n, q,obs, iterfluctcheck)
	write(43,*) maxval(erro), maxval(erro2), obs(jok+8), research_chek-51
	deallocate(erro,erro2)
	research_chek=51
	allocate(erro(n),erro2(n))
	call LR_researches(research_chek, jok, erro, erro2, eps, n, q,obs, iterfluctcheck)
	write(43,*) maxval(erro), maxval(erro2), obs(jok+8), research_chek-51
	deallocate(erro,erro2)

end do

close(31)
close(32)
close(43)

write(*,*)'work on simmetric...'


write(*,*)'error on eps...' ! 4rd
research_chek=1

open(25, file='bin/researches_res/eps_no_shift.txt')
open(26, file='bin/researches_res/eps_shift.txt')
jok=6
eps=0.1

do v=1,8

	write(*,*)v,'/ 8'
    eps=eps/10

	shift_check=0
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(25,*) maxval(erro), maxval(erro2), eps
	deallocate(erro,erro2)

	shift_check=1
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(26,*) maxval(erro), maxval(erro2), eps
	deallocate(erro,erro2)

end do
close(25)
close(26)
write(*,*)'end error on eps...' 


write(*,*)'iter on eps...' ! 5th

open(27, file='bin/researches_res/bad_iter_no_shift.txt')
open(28, file='bin/researches_res/bad_iter_shift.txt')
open(29, file='bin/researches_res/good_iter_no_shift.txt')
open(30, file='bin/researches_res/good_iter_shift.txt')
eps=0.1

do v=1,8

	write(*,*)v,'/ 8'
	eps=eps/10

	jok=3

	shift_check=0
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(27,*) q, eps
	deallocate(erro,erro2)

	shift_check=1
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(28,*) q, eps
	deallocate(erro,erro2)

	jok=6
	
	shift_check=0
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(29,*) q, eps
	deallocate(erro,erro2)

	shift_check=1
	allocate(erro(n),erro2(n))
	call iter_researches(research_chek, jok, shift_check, erro, erro2, eps, n, q, obs, iterfluctcheck)
	write(30,*) q, eps
	deallocate(erro,erro2)

end do
close(27)
close(28)
close(29)
close(30)
write(*,*)'end iter on eps...'


end program