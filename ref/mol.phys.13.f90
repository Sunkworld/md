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
  integer,parameter :: ndt=1.2d2/dt
  integer,parameter :: tt0=5d5
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

  real*8 :: eptmp, ektmp, aveep, aveek, aveek2,aveep2
  real*8 :: ep0(tt0),ek0(tt0)
  real*8 :: cor_ep(0:ndt-1,sample)=0.0d0, cor_ek(0:ndt-1,sample)=0.0d0
  real*8 :: corep(0:ndt-1), corek(0:ndt-1)
  real*8 :: t, cortimep, cortimek
  integer :: n
  do jj=-5,5
  a = -0.1d0*jj
  write(bb,'(F8.2)') a
  open(22,file='a='//trim(adjustl(bb))//'-result.txt')

  lambda = (2-2*a)/(a+1)
  b = (a+1)/2
  do ii=2,7
    h = 0.1d0*ii
    gamma = lambda/h



!    h = 0.1d0 * dt

    eqstep = 3d5/h
    tsstep = 6d6/h

    open(100,file="ei.txt")
  !  qn = q0
  !  pn = p0

  !  open(233, file='q_mp.txt')
  !  open(666, file='p_mp.txt')

  !  a = (1-gamma*h/2/m)/(1+gamma*h/2/m)
  !  b = 1/(1+gamma*h/2/m)

    do j=1, sample
       write(aa,'(I2)') j

!       open(12, file=trim(adjustl(aa))//'-x.txt')
       open(11, file=trim(adjustl(aa))//'-e.txt')

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
          ep(j)=ep(j)+0.5*qn**2-0.1*qn**3+0.1*qn**4
          ek(j)=ek(j)+0.5*pn**2/m
          if (qn>=-bound+width .and. qn<bound+width) then
             tmp=(qn+bound)/width
             bin(tmp)=bin(tmp)+1
          end if
          ep2(j)=ep2(j)+(0.5*qn**2-0.1*qn**3+0.1*qn**4)**2
          ek2(j)=ek2(j)+(0.5*pn**2/m)**2
          write(11,'(I8,F16.8,F16.8)') i, 0.5*qn**2-0.1*qn**3+0.1*qn**4, 0.5*pn**2/m
       end do
       close(11)
       ep(j) = ep(j)/tsstep
       ek(j) = ek(j)/tsstep
       ep2(j) = ep2(j)/tsstep
       ek2(j) = ek2(j)/tsstep
       write(100,'(I8,F16.8,F16.8,F16.8,F16.8)') j, ep(j), ek(j), ep2(j), ek2(j)


    end do
    close(100)


    ep_ave = sum(ep)/sample
    ep_std = sqrt(sum((ep-ep_ave)**2)/(sample-1)/sample)
    ek_ave = sum(ek)/sample
    ek_std = sqrt(sum((ek-ek_ave)**2)/(sample-1)/sample)
!    write(22,111) h,'ep=',ep_ave,ep_std,'ek=',ek_ave,ek_std
!    111 format(F7.3,2x,A5,2x,F16.8,2x,F16.8,2x,A5,2x,F16.8,2x,F16.8)

    open(100,file="ei.txt")
    do i=1,sample
       write(aa,*) i

       nt0=1

       read(100,*) l, aveep, aveek, aveep2, aveek2

       open(99,file=trim(adjustl(aa))//'-e.txt')
       read(99,*) t, ep0(1), ek0(1)
       write(*,*)trim(adjustl(aa))//'-e.txt'

       cor_ep(0,i)=cor_ep(0,i)+(ep0(1)-aveep)**2
       cor_ek(0,i)=cor_ek(0,i)+(ek0(1)-aveek)**2

       do j=1,tsstep-1
          read(99,*) t, eptmp, ektmp

          if(mod(j,ndt) .eq. 0) then
             nt0=nt0+1
             ep0(nt0)=eptmp
             ek0(nt0)=ektmp
       cor_ep(0,i)=cor_ep(0,i)+(eptmp-aveep)**2
             cor_ek(0,i)=cor_ek(0,i)+(ektmp-aveek)**2
          else

             n=mod(j,ndt)
       cor_ep(n,i)=cor_ep(n,i)+(eptmp-aveep)*(ep0(nt0)-aveep)
             cor_ek(n,i)=cor_ek(n,i)+(ektmp-aveek)*(ek0(nt0)-aveek)
    end if
       end do

       close(99)
  write(*,*)'file closed'
       cor_ep(:,i)=cor_ep(:,i)/tt0/(aveep2-aveep**2)
       cor_ek(:,i)=cor_ek(:,i)/tt0/(aveek2-aveek**2)
    enddo
    close(100)
    write(cc,'(F8.2)') h
    open(11,file="a="//trim(adjustl(bb))//"_dt="//trim(adjustl(cc))//"cor_energy.dat")

    do j=0,ndt-1
        corep(j)=sum(cor_ep(j,:))/sample
        corek(j)=sum(cor_ek(j,:))/sample
  write(11,"(I5,3(2x,F20.12))") j, j*dt, corep(j), corek(j)
    end do

    close(11)
    cortimep = sum(corep)
    cortimek = sum(corek)
    write(22,'(F8.2,F16.8,F16.8)') h, ep_ave, cortimep
!    do i=1, nbin
!       x=-bound+width*i+width/2
!       dx=dble(bin(i))/dble(tsstep)/dble(sample)/width ! numerical distribution of x
!       dx0=1d0/2.17433d0*exp(-(0.5*x**2-0.1*x**3+0.1*x**4)/kT) ! exact distribution of x
!       write(12,*) x, dx, dx0
!    end do
  enddo
  close(22)
enddo

end subroutine molphys
