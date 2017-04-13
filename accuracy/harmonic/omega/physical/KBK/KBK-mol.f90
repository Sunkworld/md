module init
  implicit none
  real(8), parameter :: beta = 1d0                                           !needs to be modified
  real(8), parameter :: h = 0.1d0      
  real(8), parameter :: k = 1d0
  real(8), parameter :: m = 1d0
  integer, parameter :: eqstep=2d7/h
  integer, parameter :: tsstep=1d7/h
  integer, parameter :: sample=20
    real*8, parameter :: bound=3.0d0 ! plotting distribution density function from -bound to bound
    real*8, parameter :: wihh=0.1d0 ! wihh of each bin
    integer, parameter :: nbin=2*bound/wihh ! number of bins
    real*8 :: sigma(2,2),CC(2,2)
end module init

program main
  call molphys
end program main
subroutine calForce(fn, x)
  use init, only: k
  implicit none
  real(8) :: fn, x
  fn = -k*x                                                                              !needs to be modified
end subroutine calForce

subroutine molphys
  use init
  use random
  implicit none
  real(8) :: rand1,rand2, qn, pn, qnp1, pnp1, fn, fnp1, a, b
  integer :: i, j
  real*8 :: eptmp, ektmp, ep(sample), ek(sample),sigmax(sample),sigmap(sample),sigmaxp(sample),x2ave,xave,p2ave,pave,xpave,ep_ave,ek_ave,ep_std,ek_std,wxx,wyy,wxy
  real*8 :: cor_ep(tsstep/2),cor_ek(tsstep/2)
  real*8 :: t, cortimep, cortimek,cortimep_std,cortimek_std
    real*8 :: eigval(2),eigvec(2,2)
  integer :: n!, nh   
  character(30) :: c
  real*8 :: ettmp,et(sample),et_ave,et_std
  real*8 :: ekk, epp,ett

  real*8 :: gamma = 1d0
  real*8 :: c1
  c1=exp(-gamma*h*0.5d0)

    sigma(1,1) = 1d0/(m*beta*gamma**2)*(gamma*h-3d0+4d0*c1-exp(-gamma*h))
    sigma(1,2) = 1d0/(beta*gamma)*(1d0-c1)**2
    sigma(2,1) = 1d0/(beta*gamma)*(1d0-c1)**2
    sigma(2,2) = m/beta*(1d0-exp(-gamma*h))
    call mat_diag(eigval, eigvec, sigma, 2)

    CC(1,1) = eigvec(1,1)*sqrt(eigval(1))
    CC(1,2) = eigvec(1,2)*sqrt(eigval(2))
    CC(2,1) = eigvec(2,1)*sqrt(eigval(1))
    CC(2,2) = eigvec(2,2)*sqrt(eigval(2))

  open(22,file='result.maindat')
  ep(:)=0
  ek(:)=0
    write(*,*) 'gamma=',gamma, 'dt=', h
    do j=1, sample
       write(c,'(I2)') j
       write(*,*) 'Sample=', j
       if ((j==1) .and. (h==0.1d0 .or. h==0.2d0)) then
       open(33,file=trim('miao_'//adjustl(c)))
       endif
       call random_normal(rand1)
       pn = rand1
       call random_number(rand1)
       qn = 4d0*(rand1-0.5d0)
    !   write(*,*) 'sample=', j
      do i=1, eqstep
          call random_normal(rand1)
	  call random_normal(rand2)

	  qn = qn + pn/(m*gamma)*(1d0-c1) + CC(1,1)/sqrt(m)*Rand1 + CC(1,2)/sqrt(m)*Rand2
          pn = c1*pn + CC(2,1)*sqrt(m)*Rand1 + CC(2,2)*sqrt(m)*Rand2    !needs to be modified
	
          call calForce(fn, qn)
          pn = pn + h*fn
          call random_normal(rand1)
	  call random_normal(rand2)

	  qn = qn + pn/(m*gamma)*(1d0-c1) + CC(1,1)/sqrt(m)*Rand1 + CC(1,2)/sqrt(m)*Rand2
          pn = c1*pn + CC(2,1)*sqrt(m)*Rand1 + CC(2,2)*sqrt(m)*Rand2    !needs to be modified
	!         if (mod(i, tsstep/10+1) .eq. 0) then
!             write(*,*) real(i)/real(tsstep)*100, '%'
!             write(*,*) qn, pn
!         end if
         eptmp = 0.5*k*qn**2                                !needs to be modified
         ektmp = 0.5*pn**2/m
         ettmp = eptmp + ektmp
!         ep(j) = ep(j)+eptmp/tsstep
!         ek(j) = ek(j) + ektmp/tsstep
if ((j==1) .and. (h==0.1d0 .or. h==0.2d0)) then
	           epp = eptmp/i+epp*(i-1)/i
           ekk = ektmp/i+ekk*(i-1)/i
           ett = ettmp/i+ett*(i-1)/i
           if (mod(i,200) .eq. 1) then
                     write(33,'(F24.8,F16.8,F16.8,F16.8)') i*h,epp,ekk,ett
           endif
           endif
enddo

      do i=eqstep+1, eqstep+tsstep
          call random_normal(rand1)
	  call random_normal(rand2)

	  qn = qn + pn/(m*gamma)*(1d0-c1) + CC(1,1)/sqrt(m)*Rand1 + CC(1,2)/sqrt(m)*Rand2
          pn = c1*pn + CC(2,1)*sqrt(m)*Rand1 + CC(2,2)*sqrt(m)*Rand2    !needs to be modified
	
          call calForce(fn, qn)
          pn = pn + h*fn
          call random_normal(rand1)
	  call random_normal(rand2)

	  qn = qn + pn/(m*gamma)*(1d0-c1) + CC(1,1)/sqrt(m)*Rand1 + CC(1,2)/sqrt(m)*Rand2
          pn = c1*pn + CC(2,1)*sqrt(m)*Rand1 + CC(2,2)*sqrt(m)*Rand2    !needs to be modified
	        if (mod((i-eqstep), tsstep/10+1) .eq. 0) then
             write(*,*) real((i-eqstep))/real(tsstep)*100, '%'
             write(*,*) qn, pn
         end if
         eptmp = 0.5*k*qn**2                                !needs to be modified
         ektmp = 0.5*pn**2/m
         ettmp = eptmp + ektmp
         ep(j) = ep(j)+eptmp/tsstep
         ek(j) = ek(j) + ektmp/tsstep
         et(j) = et(j) + ettmp/tsstep
	           epp = eptmp/i+epp*(i-1)/i
           ekk = ektmp/i+ekk*(i-1)/i
           ett = ettmp/i+ett*(i-1)/i
            if ((j==1) .and. (h==0.1d0 .or. h==0.2d0)) then
           if (mod(i,200) .eq. 1) then
                     write(33,'(F24.8,F16.8,F16.8,F16.8)') i*h,epp,ekk,ett
           endif
endif
        enddo
            if ((j==1) .and. (h==0.1d0 .or. h==0.2d0)) then
      close(33)
    endif

   enddo

   ep_ave = sum(ep)/sample
   ep_std = sqrt(sum((ep-ep_ave)**2)/(sample-1)/sample)
   ek_ave = sum(ek)/sample
   ek_std = sqrt(sum((ek-ek_ave)**2)/(sample-1)/sample)
    et_ave = sum(et)/sample
   et_std = sqrt(sum((et-et_ave)**2)/(sample-1)/sample)
   write(22,'(F8.2,F8.2,6F16.8)') h,gamma,ep_ave,ep_std,ek_ave,ek_std,et_ave,et_std
end subroutine molphys

subroutine mat_diag(eigval, eigvec, A, n)
    implicit none

    integer, intent(in) :: n
    real(8), intent(in) :: A(n,n)
    real(8), intent(out) :: eigval(n)
    real(8), intent(out) :: eigvec(n,n)

    integer :: lwork, info
    real(8) :: work(n*(3+n/2))
    lwork=n*(3+n/2)

    eigvec = A
    ! use lapack routine DSYEV to get eigenvalues and eigenvectors
    call dsyev('V','U',n,eigvec,n,eigval,work,lwork,info)

    if (info .gt. 0) then
        write(*,*)'Error: matrix diagonalization failed'
        stop
    end if

    return
end subroutine mat_diag

