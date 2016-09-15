module init
  implicit none
  real(8), parameter :: kT = 1d0
  real(8), parameter :: gamma = 10d0
  real(8), parameter :: h=0.75d0
  real(8), parameter :: k = 1d0
  integer, parameter :: Nstep = 1d7
!  real(8), parameter :: q0 = -1d0
!  real(8), parameter :: p0 = -0.1d0
  real(8), parameter :: m = 1d0
  integer, parameter :: eqstep=1d5/h
  integer, parameter :: tsstep=1d7/h
  integer, parameter :: sample=20
  real(8), parameter :: bound=3d0
  real(8), parameter :: width=0.1d0
  integer, parameter :: nbin=2*bound/width
  real(8) :: x, dx, dx0
end module init

program main



  call molphys
!  call BAOAB
!  call vv
end program main
subroutine calForce(fn, x)
  use init, only: k
  implicit none
  real(8) :: fn, x
  fn = -k*x**3
end subroutine calForce
subroutine BAOAB
  use init
  use random
  implicit none
  real(8) :: rand, qn, pn, qnp1, pnp1, fn, ph, qh, c
  integer :: i, j, tmp
  integer :: bin(nbin)=0
  !  real(8) :: eps(sample)=0d0, eks(sample)=0d0
  real(8) :: ep(sample)=0d0, ek(sample)=0d0
  real(8) :: ep_std=0d0, ek_std=0d0, ep_ave=0d0, ek_ave=0d0
  character(len=30) :: aa
  
  write(aa,'(F8.2)') h
  
  open(12, file='dt='//trim(adjustl(aa))//'-baoab-x.txt')
  open(11, file='dt='//trim(adjustl(aa))//'-baoab-e.txt')

  
  call random_normal(rand)
  pn = rand
  call random_number(rand)
  qn = 4d0*(rand-0.5d0)
  call calForce(fn, qn)
!  open(2333, file='q_baoab.txt')
  do j=1, sample
     write(*,*) 'sample=', j
     do i=1, eqstep
        ph = pn + 0.5*h*fn
        qh = qn + 0.5*h*ph/m
        c = exp(-gamma*h)
        call random_normal(rand)
        ph = c*ph + sqrt((1-c*c)*kT)*sqrt(m)*rand
        qnp1 = qh + 0.5*h*ph/m
        call calForce(fn, qnp1)
        pnp1 = ph + 0.5*h*fn
        pn = pnp1
        qn = qnp1
        !     write(2333,*) qn

     end do
     do i=1, tsstep
        ph = pn + 0.5*h*fn
        qh = qn + 0.5*h*ph/m
        c = exp(-gamma*h)
        call random_normal(rand)
        ph = c*ph + sqrt((1-c*c)*kT)*sqrt(m)*rand
        qnp1 = qh + 0.5*h*ph/m
        call calForce(fn, qnp1)
        pnp1 = ph + 0.5*h*fn
        pn = pnp1
        qn = qnp1
        if (mod(i, tsstep/10+1) .eq. 0) then
           write(*,*) real(i)/real(tsstep)*100, '%'
           write(*,*) qn, pn
        end if
        ep(j)=ep(j)+0.25*k*qn**4
        ek(j)=ek(j)+0.5*pn**2/m
        if (qn>=-bound+width .and. qn<bound+width) then
           tmp=(qn+bound)/width
           bin(tmp)=bin(tmp)+1
        end if


     end do
     ep(j) = ep(j)/tsstep
     ek(j) = ek(j)/tsstep

  end do
  ep_ave = sum(ep)/sample
  ep_std = sqrt(sum((ep-ep_ave)**2)/(sample-1)/sample)
  ek_ave = sum(ek)/sample
  ek_std = sqrt(sum((ek-ek_ave)**2)/(sample-1)/sample)
  write(11,111) h,'ep=',ep_ave,ep_std,'ek=',ek_ave,ek_std
111 format(F7.3,2x,A5,2x,F16.8,2x,F16.8,2x,A5,2x,F16.8,2x,F16.8)
  do i=1, nbin
     x=-bound+width*i+width/2
     dx=dble(bin(i))/dble(tsstep)/dble(sample)/width ! numerical distribution of x
     dx0=1d0/2.56369d0*exp(-k*x**4/4/kT) ! exact distribution of x
     write(12,*) x, dx, dx0
  end do

end subroutine BAOAB
  
subroutine molphys
  use init
  use random
  implicit none
  real(8) :: rand, qn, pn, qnp1, pnp1, fn, fnp1, a, b
  integer :: i, j, tmp
  integer :: bin(nbin)=0
  !  real(8) :: eps(sample)=0d0, eks(sample)=0d0
  real(8) :: ep(sample)=0d0, ek(sample)=0d0
  real(8) :: ep_std=0d0, ek_std=0d0, ep_ave=0d0, ek_ave=0d0
  character(len=30) :: aa

  write(aa,'(F8.2)') h
  
  open(12, file='dt='//trim(adjustl(aa))//'-mp-x.txt')
  open(11, file='dt='//trim(adjustl(aa))//'-mp-e.txt')

!  qn = q0
!  pn = p0

!  open(233, file='q_mp.txt')
!  open(666, file='p_mp.txt')

  a = (1-gamma*h/2/m)/(1+gamma*h/2/m)
  b = 1/(1+gamma*h/2/m)
  
  call random_normal(rand)
  pn = rand
  call random_number(rand)
  qn = 4d0*(rand-0.5d0)
  call calForce(fn, qn)
  do j=1, sample
     write(*,*) 'sample=', j
     do i=1, eqstep
        call random_normal(rand)
        qnp1 = qn + b*h*pn/m + b*h**2/2/m*fn + b*h/2/m*rand*sqrt(2*gamma*kT*h)
        call calForce(fnp1, qnp1)
        pnp1 = a*pn/m + h/2/m*(a*fn+fnp1) + b/m*rand*sqrt(2*gamma*kT*h)
        qn = qnp1
        pn = pnp1
        fn = fnp1
        !     write(233, *) qn
        !     write(666, *) pn

     end do
     
     do i=1, tsstep
        call random_normal(rand)
        qnp1 = qn + b*h*pn/m + b*h**2/2/m*fn + b*h/2/m*rand*sqrt(2*gamma*kT*h)
        call calForce(fnp1, qnp1)
        pnp1 = a*pn/m + h/2/m*(a*fn+fnp1) + b/m*rand*sqrt(2*gamma*kT*h)
        qn = qnp1
        pn = pnp1
        fn = fnp1
        !     write(233, *) qn
        !     write(666, *) pn
        if (mod(i, tsstep/10+1) .eq. 0) then
           write(*,*) real(i)/real(tsstep)*100, '%'
           write(*,*) qn, pn
        end if
        ep(j)=ep(j)+0.25*k*qn**4
        ek(j)=ek(j)+0.5*pn**2/m
        if (qn>=-bound+width .and. qn<bound+width) then
           tmp=(qn+bound)/width
           bin(tmp)=bin(tmp)+1
        end if

     end do
     ep(j) = ep(j)/tsstep
     ek(j) = ek(j)/tsstep
     
  end do
  ep_ave = sum(ep)/sample
  ep_std = sqrt(sum((ep-ep_ave)**2)/(sample-1)/sample)
  ek_ave = sum(ek)/sample
  ek_std = sqrt(sum((ek-ek_ave)**2)/(sample-1)/sample)
  write(11,111) h,'ep=',ep_ave,ep_std,'ek=',ek_ave,ek_std
  111 format(F7.3,2x,A5,2x,F16.8,2x,F16.8,2x,A5,2x,F16.8,2x,F16.8)
  do i=1, nbin
     x=-bound+width*i+width/2
     dx=dble(bin(i))/dble(tsstep)/dble(sample)/width ! numerical distribution of x
     dx0=1d0/2.56369d0*exp(-k*x**4/4/kT) ! exact distribution of x
     write(12,*) x, dx, dx0
  end do


end subroutine molphys







