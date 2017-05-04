module init
  implicit none
  real(8), parameter :: kT = 1d0
  real(8), parameter :: h = 0.1d0                                                       !needs to be modified
  real(8), parameter :: k = 1d0
  real(8), parameter :: m = 1d0
  integer, parameter :: eqstep=3d7/h
  integer, parameter :: tsstep=1d7/h
  integer, parameter :: sample=20
end module init

program main
  call molphys
end program main
subroutine calForce(fn, x)
  use init, only: k
  implicit none
  real(8) :: fn, x
  fn = -x
end subroutine calForce

subroutine molphys
  use init
  use random
  implicit none
  real(8) :: rand, qn, pn, qnp1, pnp1, fn, fnp1, a, b
  integer :: i, j
  real*8 :: eptmp, ektmp, ep(sample), ek(sample),ep_ave,ek_ave,ep_std,ek_std
  real*8 :: ettmp,et(sample),et_ave,et_std
  integer :: n!, ndt
  character(30) :: c
  real*8 :: epp, ekk,ett
  real(8) :: gamma
  if (h .lt. 0.85) then
  gamma =1d0/h*log((-128d0+32d0*h**4-12d0*h**6+h**8+h*(h*h-8d0)*sqrt(1024d0-1280d0*h*h+512d0*h**4-48d0*h**6-8d0*h**8+h**10))/(-128d0+256d0*h*h-96d0*h**4+8d0*h**6))
else
    gamma = 6.42174990024850d0
endif
  open(22,file='result.maindat')
  ep(:)=0
  ek(:)=0
  ekk=0
  epp=0
    a = -exp(-gamma*h)
    write(*,*) 'gamma=',gamma, 'dt=', h
    do j=1, sample
       write(c,'(I2)') j
       write(*,*) 'Sample=', j
!       open(33,file=trim('traj_'//adjustl(c)))
       if ((j==1) .and. (h==0.1d0)) then
           open(33,file='miao')
       endif
       call random_normal(rand)
       pn = rand
       call random_number(rand)
       qn = 4d0*(rand-0.5d0)
       call calForce(fn, qn)
    !   write(*,*) 'sample=', j
        do i=1, eqstep



                    call random_normal(rand)
          pn = a*pn + sqrt((1-a*a)*kT)*sqrt(m)*rand 

          pn = pn + 0.5d0*h*fn

          qn = qn + h*pn/m
                                              
          call calForce(fn, qn)

          pn = pn + 0.5d0*h*fn

          call random_normal(rand)
          pn = a*pn + sqrt((1-a*a)*kT)*sqrt(m)*rand 


         eptmp = 0.5d0*qn**2
         ektmp = 0.5d0*pn**2/m
         ettmp = eptmp + ektmp
!         ep(j) = ep(j)+eptmp/eqstep
!         ek(j) = ek(j) + ektmp/eqstep
         if ((j==1) .and. (h==0.1d0)) then
           epp = eptmp/i+epp*(i-1)/i
           ekk = ektmp/i+ekk*(i-1)/i
           ett = ettmp/i+ett*(i-1)/i
           if (mod(i,50) .eq. 1) then
                     write(33,'(F24.8,F16.8,F16.8,F16.8)') i*h,epp,ekk,ett
           endif
         endif
        enddo
       do i=eqstep+1, eqstep+tsstep




                    call random_normal(rand)
          pn = a*pn + sqrt((1-a*a)*kT)*sqrt(m)*rand 

          pn = pn + 0.5d0*h*fn

          qn = qn + h*pn/m
                                              
          call calForce(fn, qn)

          pn = pn + 0.5d0*h*fn

          call random_normal(rand)
          pn = a*pn + sqrt((1-a*a)*kT)*sqrt(m)*rand 




         if (mod((i-eqstep), tsstep/10+1) .eq. 0) then
             write(*,*) real((i-eqstep))/real(tsstep)*100, '%'
             write(*,*) qn, pn
         end if
         eptmp = 0.5d0*qn**2
         ektmp = 0.5d0*pn**2/m
         ettmp = eptmp + ektmp
         ep(j) = ep(j)+eptmp/tsstep
         ek(j) = ek(j) + ektmp/tsstep
         et(j) = et(j) + ettmp/tsstep
           epp = eptmp/i+epp*(i-1)/i
           ekk = ektmp/i+ekk*(i-1)/i
           ett = ettmp/i+ett*(i-1)/i
       if ((j==1) .and. (h==0.1d0)) then
           if (mod(i,50) .eq. 1) then
                     write(33,'(F24.8,F16.8,F16.8,F16.8)') i*h,epp,ekk,ett
           endif
       endif
        enddo
       if ((j==1) .and. (h==0.1d0)) then
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

