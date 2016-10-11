module init
  implicit none
  real(8), parameter :: kT = 1d0
!  real(8), parameter :: gamma=0.5d0
!  real(8), parameter :: h = 1.1d0
  integer, parameter :: num = 3
  real(8), parameter :: k = 1d0
!  integer, parameter :: Nstep = 1d7 !Nonsense
!  real(8), parameter :: q0 = -1d0
!  real(8), parameter :: p0 = -0.1d0
  real(8), parameter :: m = 1d0
!  integer, parameter :: eqstep=1d5/h
!  integer, parameter :: tsstep=1d7/h
  integer, parameter :: sample=20
  real(8), parameter :: bound=3d0
  real(8), parameter :: width = 0.1d0
  integer, parameter :: nbin=2*bound/width
!  integer,parameter :: ndt=1.2d2/dt
  integer,parameter :: tt0=5d4
  integer, parameter :: tottime = 6d6
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
  fn = -x+0.3*x**2-0.4*x**3
end subroutine calForce

subroutine molphys
  use init
  use random
  implicit none
  real(8) :: rand, qn, pn, qnp1, pnp1, fn, fnp1, a, b
  integer :: i, j, tmp, dt, eqstep, tsstep, l
  integer :: bin(nbin)=0
  !  real(8) :: eps(sample)=0d0, eks(sample)=0d0
  real(8) :: ep(sample)=0d0, ek(sample)=0d0, ep2(sample)=0d0, ek2(sample)=0d0
  real(8) :: ep_std=0d0, ek_std=0d0, ep_ave=0d0, ek_ave=0d0
  character(len=30) :: aa, bb, cc, dd
  real(8) :: h, gamma,lambda
  integer :: ii,jj,kk

  integer :: nt0=1
  real*8 :: eptmp, ektmp, aveep, aveek, aveek2,aveep2
  real*8 :: ep0(tt0),ek0(tt0)
  real*8,allocatable :: cor_ep(:,:), cor_ek(:,:)
  real*8,allocatable :: corep(:), corek(:)
  real*8 :: t, cortimep, cortimek
  integer :: n, ndt
  do jj=-6,6
  a = 0.1d0*jj
  write(bb,'(F8.2)') a
  !open(22,file='a='//trim(adjustl(bb))//'-result.maindat')

  lambda = (2-2*a)/(a+1)
  b = (a+1)/2
  do ii=3,10
    ep(:)=0
    ek(:)=0
    ep2(:)=0
    ek2(:)=0
    open(22,file='result.maindat',position='append')
    h = 0.1d0*ii
    gamma = lambda/h
    ndt = ceiling(120d0/h)
    write(*,*) h,ndt
    allocate(cor_ep(0:ndt-1,sample),cor_ek(0:ndt-1,sample),corep(0:ndt-1),corek(0:ndt-1))
    cor_ep(:,:)=0d0
    cor_ek(:,:)=0d0
!    h = 0.1d0 * dt
    write(*,*) 'a=',a, 'dt=', h
    eqstep = 3d5/h
    tsstep = tottime/h
  !  tt0 = tsstep/ndt
  !  open(100,file="ei.txt")
  !  qn = q0
  !  pn = p0

  !  open(233, file='q_mp.txt')
  !  open(666, file='p_mp.txt')

  !  a = (1-gamma*h/2/m)/(1+gamma*h/2/m)
  !  b = 1/(1+gamma*h/2/m)

    do j=1, sample
       write(aa,'(I2)') j
       nt0=0
!       open(12, file=trim(adjustl(aa))//'-x.txt')
    !   open(11, file=trim(adjustl(aa))//'-e.txt')

       call random_normal(rand)
       pn = rand
       call random_number(rand)
       qn = 4d0*(rand-0.5d0)
       call calForce(fn, qn)
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
          eptmp = 0.5*qn**2-0.1*qn**3+0.1*qn**4
          ektmp = 0.5*pn**2/m
          if(mod(i-1,ndt) .eq. 0) then
            nt0=nt0+1
      !      write(*,*) i-1
            ep0(nt0)=eptmp
            ek0(nt0)=ektmp
            cor_ep(0,j) = cor_ep(0,j)+eptmp**2
            cor_ek(0,j) = cor_ek(0,j)+ektmp**2
    !        write(*,*) cor_ep(0,j)
  !          write(*,*) nt0,eptmp
          else
            n=mod(i-1,ndt)
            cor_ep(n,j)=cor_ep(n,j)+eptmp*ep0(nt0)
            cor_ek(n,j)=cor_ek(n,j)+ektmp*ek0(nt0)

          endif

          ep(j)=ep(j)+eptmp
          ek(j)=ek(j)+ektmp
    !      if (qn>=-bound+width .and. qn<bound+width) then
  !           tmp=(qn+bound)/width
    !         bin(tmp)=bin(tmp)+1
    !      end if
          ep2(j)=ep2(j)+(eptmp)**2
          ek2(j)=ek2(j)+(ektmp)**2
          !write(11,'(I8,F16.8,F16.8)') i, 0.5*qn**2-0.1*qn**3+0.1*qn**4, 0.5*pn**2/m
       end do

    !   close(11)
       ep(j) = ep(j)/tsstep
       ek(j) = ek(j)/tsstep
       ep2(j) = ep2(j)/tsstep
       ek2(j) = ek2(j)/tsstep

       cor_ep(:,j)=(cor_ep(:,j)/tt0-ep(j)**2)/(ep2(j)-ep(j)**2)
       cor_ek(:,j)=(cor_ek(:,j)/tt0-ep(j)**2)/(ep2(j)-ep(j)**2)
  !     write(100,'(I8,F16.8,F16.8,F16.8,F16.8)') j, ep(j), ek(j), ep2(j), ek2(j)
    write(*,*) cor_ep(0,j)

    end do
    !close(100)
    write(cc,'(F8.2)') h
!    open(11,file="a="//trim(adjustl(bb))//"_dt="//trim(adjustl(cc))//"_cor_energy.dat")
    do i=0,ndt-1
        corep(i)=sum(cor_ep(i,:))/sample
        corek(i)=sum(cor_ek(i,:))/sample

!  write(11,"(I5,3(2x,F20.12))") i, i*h, corep(i), corek(i)
    end do
  !  close(11)
    cortimep = sum(corep)
    cortimek = sum(corek)


    ep_ave = sum(ep)/sample
    ep_std = sqrt(sum((ep-ep_ave)**2)/(sample-1)/sample)
    ek_ave = sum(ek)/sample
    ek_std = sqrt(sum((ek-ek_ave)**2)/(sample-1)/sample)
    write(22,'(F8.2,F8.2,F8.2,F16.8,F16.8)') a,h,gamma, ep_ave, cortimep

!    write(22,111) h,'ep=',ep_ave,ep_std,'ek=',ek_ave,ek_std
!    111 format(F7.3,2x,A5,2x,F16.8,2x,F16.8,2x,A5,2x,F16.8,2x,F16.8)
    deallocate(cor_ep,cor_ek,corep,corek)

	close(22)
  enddo

enddo
end subroutine molphys
