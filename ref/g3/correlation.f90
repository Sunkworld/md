module constant
    implicit none

    real*8,parameter :: dt=0.3d0
    integer,parameter :: sample=20
    integer,parameter :: tsstep=6d6/dt
    integer,parameter :: ndt=1.2d2/dt
    integer,parameter :: tt0=tsstep/ndt

end module

program main
    use constant
    implicit none
    integer :: i,j,l
    integer :: nt0=1          !t0 counter
    real*8 :: ep, ek, aveep, aveek, aveek2,aveep2
    real*8 :: ep0(tt0),ek0(tt0)
    real*8 :: cor_ep(0:ndt-1,sample)=0.0d0, cor_ek(0:ndt-1,sample)=0.0d0
    real*8 :: corep(0:ndt-1), corek(0:ndt-1)
    real*8 :: t
    integer :: n
    character(len=30) :: aa

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
          read(99,*) t, ep, ek

          if(mod(j,ndt) .eq. 0) then
             nt0=nt0+1
             ep0(nt0)=ep
             ek0(nt0)=ek
	     cor_ep(0,i)=cor_ep(0,i)+(ep-aveep)**2
             cor_ek(0,i)=cor_ek(0,i)+(ek-aveek)**2
          else

             n=mod(j,ndt)
	     cor_ep(n,i)=cor_ep(n,i)+(ep-aveep)*(ep0(nt0)-aveep)
             cor_ek(n,i)=cor_ek(n,i)+(ek-aveek)*(ek0(nt0)-aveek)
	  end if
       end do

       close(99)
write(*,*)'file closed'
       cor_ep(:,i)=cor_ep(:,i)/tt0/(aveep2-aveep**2)
       cor_ek(:,i)=cor_ek(:,i)/tt0/(aveek2-aveek**2)
    enddo
    close(100)

    open(11,file="cor_energy.dat")

    do j=0,ndt-1
        corep(j)=sum(cor_ep(j,:))/sample
        corek(j)=sum(cor_ek(j,:))/sample
	write(11,"(I5,3(2x,F20.12))") j, j*dt, corep(j), corek(j)
    end do

    close(11)

    stop
end program main
