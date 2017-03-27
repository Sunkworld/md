module init
  implicit none
  real(8), parameter :: kT = 1d0
!  real(8), parameter :: gamma=0.5d0
  real(8), parameter :: h = 0.05d0                                                       !needs to be modified
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
  integer,parameter :: ndt=200000!needs to be modified
  integer, parameter :: tottime = 1d7
  integer,parameter :: tt0=floor(tottime/ndt/h)
  real(8) :: x, dx
end module init

program main

! test3

  call molphys
!  call BAOAB
!  call vv
end program main
subroutine calForce(fn, x)
  use init, only: k
  implicit none
  real(8) :: fn, x
  fn =  -x
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
  real(8) :: et(sample)=0d0, et2(sample)=0d0
  real(8) :: ep_std=0d0, ek_std=0d0, et_std=0d0
  real(8) :: ep_ave=0d0, ek_ave=0d0, et_ave=0d0
  character(len=30) :: aa, bb, cc, dd
  real(8) :: gamma,lambda,gamma2!,h
  integer :: ii,jj,kk

  integer :: nt0=1
  real*8 :: eptmp, ektmp, ettmp, aveep, aveek, aveet, aveek2,aveep2, aveet2
  real*8 :: ep0(tt0),ek0(tt0), et0(tt0)
  real*8,allocatable :: cor_ep(:,:), cor_ek(:,:), cor_et(:,:)
  real*8,allocatable :: corep(:), corek(:), coret(:)
  real*8 :: t, cortimep, cortimek,cortimet,cortimep_std,cortimek_std,totalp,totalk,totalt
  integer :: n!, ndt
  do jj=6,6
!  a = 0.1d0*jj
  write(bb,'(F8.2)') a
  !open(22,file='a='//trim(adjustl(bb))//'-result.maindat')

!  lambda = (2-2*a)/(a+1)
!  b = (a+1)/2
  do ii=1,1
    ep(:)=0
    ek(:)=0
    ep2(:)=0
    ek2(:)=0
    open(22,file='result.maindat',position='append')
!    h = 0.02d0
    gamma2 = 0.4d0                                                                      !needs to be modified
    a = -exp(-gamma2*h)
    lambda = (2-2*a)/(a+1)
    b = (a+1)/2
    gamma = lambda/h
!    ndt = ceiling(tottime/tt0/h)
    write(*,*) 'Assumed correlation time=',ndt
    allocate(cor_ep(0:ndt-1,sample),cor_ek(0:ndt-1,sample),corep(sample),corek(sample))
    allocate(cor_et(0:ndt-1,sample),coret(sample))
    cor_ep(:,:)=0d0
    cor_ek(:,:)=0d0
!    h = 0.05d0 * dt
    write(*,*) 'gamma=',gamma2, 'dt=', h
    eqstep = 1d5/h
    tsstep = tottime/h
  !  tt0 = tsstep/ndt
  !  open(100,file="ei.txt")
  !  qn = q0
  !  pn = p0

    open(233, file='ct.txt')
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
    !   write(*,*) 'sample=', j
       do i=1, eqstep
          pn = pn + 0.5*h*fn

          qn = qn + h*pn/m
                                              
          call calForce(fn, qn)

          pn = pn + 0.5*h*fn

          call random_normal(rand)
          pn = a*pn + sqrt((1-a*a)*kT)*sqrt(m)*rand 
          !     write(233, *) qn
          !     write(666, *) pn

       end do

       do i=1, tsstep
          pn = pn + 0.5*h*fn

          qn = qn + h*pn/m
                                              
          call calForce(fn, qn)

          pn = pn + 0.5*h*fn

          call random_normal(rand)
          pn = a*pn + sqrt((1-a*a)*kT)*sqrt(m)*rand 
          !     write(233, *) qn
          !     write(666, *) pn
      !    if (mod(i, tsstep/10+1) .eq. 0) then
      !       write(*,*) real(i)/real(tsstep)*100, '%'
      !       write(*,*) qn, pn
      !    end if
         eptmp = 0.5*qn**2                                                           !needs to be modified
          ektmp = 0.5*pn**2/m
          ettmp = eptmp + ektmp
          if(mod(i-1,ndt) .eq. 0) then
            nt0=nt0+1
            if (nt0 .gt. tt0) exit
      !      write(*,*) i-1
            ep0(nt0)=eptmp
            ek0(nt0)=ektmp
            et0(nt0)=ettmp
            cor_ep(0,j) = cor_ep(0,j)+eptmp**2
            cor_ek(0,j) = cor_ek(0,j)+ektmp**2
            cor_et(0,j) = cor_et(0,j)+ettmp**2
    !        write(*,*) cor_ep(0,j)
  !          write(*,*) nt0,eptmp
          else
            n=mod(i-1,ndt)
            cor_ep(n,j)=cor_ep(n,j)+eptmp*ep0(nt0)
            cor_ek(n,j)=cor_ek(n,j)+ektmp*ek0(nt0)
            cor_et(n,j)=cor_et(n,j)+ettmp*et0(nt0)

          endif

          ep(j)=ep(j)+eptmp
          ek(j)=ek(j)+ektmp
          et(j)=et(j)+ettmp
    !      if (qn>=-bound+width .and. qn<bound+width) then
  !           tmp=(qn+bound)/width
    !         bin(tmp)=bin(tmp)+1
    !      end if
          ep2(j)=ep2(j)+(eptmp)**2
          ek2(j)=ek2(j)+(ektmp)**2
          et2(j)=et2(j)+ettmp**2
          !write(11,'(I8,F16.8,F16.8)') i, 0.5*qn**2-0.1*qn**3+0.1*qn**4, 0.5*pn**2/m
       end do

    !   close(11)
       ep(j) = ep(j)/tsstep
       ek(j) = ek(j)/tsstep
       et(j) = et(j)/tsstep
       ep2(j) = ep2(j)/tsstep
       ek2(j) = ek2(j)/tsstep
    et2(j) = et2(j)/tsstep

       cor_ep(:,j)=(cor_ep(:,j)/tt0-ep(j)**2)/(ep2(j)-ep(j)**2)
       cor_ek(:,j)=(cor_ek(:,j)/tt0-ek(j)**2)/(ek2(j)-ek(j)**2)
       cor_et(:,j)=(cor_et(:,j)/tt0-et(j)**2)/(et2(j)-et(j)**2)
  !     write(100,'(I8,F16.8,F16.8,F16.8,F16.8)') j, ep(j), ek(j), ep2(j), ek2(j)
    write(*,*) cor_ep(0,j)

    end do
    !close(100)
    write(cc,'(F8.2)') h
!    open(11,file="a="//trim(adjustl(bb))//"_dt="//trim(adjustl(cc))//"_cor_energy.dat")
    do i=1,sample
        corep(i)=sum(cor_ep(:,i))
        corek(i)=sum(cor_ek(:,i))
        coret(i)=sum(cor_et(:,i))

!  write(11,"(I5,3(2x,F20.12))") i, i*h, corep(i), corek(i)
     end do
     totalp = 0
     totalk = 0
     totalt = 0
     do i=0,ndt-1
     totalp = totalp + sum(cor_ep(i,:))/sample
     totalk = totalk + sum(cor_ek(i,:))/sample
     totalt = totalt + sum(cor_et(i,:))/sample
        write(233,'(I8,F16.8,F16.8,F16.8)') i,totalk,totalp,totalt !sum(cor_ep(i,:))/sample
        enddo
  !  close(11)
    cortimep = sum(corep)/sample
    cortimek = sum(corek)/sample
    write(*,*) corep(1),cortimep
    cortimep_std = sqrt(sum((corep-cortimep)**2)/(sample-1)/sample)
    cortimek_std = sqrt(sum((corek-cortimek)**2)/(sample-1)/sample)

    ep_ave = sum(ep)/sample
    ep_std = sqrt(sum((ep-ep_ave)**2)/(sample-1)/sample)
    ek_ave = sum(ek)/sample
    ek_std = sqrt(sum((ek-ek_ave)**2)/(sample-1)/sample)
    write(22,'(A8,A8,A16,A16,A16,A16,A16,A16,A16,A16)') 'dt','gamma','ep_ave','ep_std','cortimep','cortimep_std','ek_ave',&
         &'ek_std','cortimek','cortimek_std'
    write(22,'(F8.2,F8.2,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)') h,gamma2, ep_ave,ep_std,cortimep,cortimep_std,&
         &ek_ave,ek_std, cortimek,cortimek_std

!    write(22,111) h,'ep=',ep_ave,ep_std,'ek=',ek_ave,ek_std
!    111 format(F7.3,2x,A5,2x,F16.8,2x,F16.8,2x,A5,2x,F16.8,2x,F16.8)
    deallocate(cor_ep,cor_ek,cor_et,coret,corep,corek)

	close(22)
  enddo

enddo
end subroutine molphys
