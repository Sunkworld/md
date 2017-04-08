program main
   implicit none
   integer :: i, j
   integer,parameter :: n = 20000000
   integer,parameter :: nsample = 20
   integer,parameter :: ncor = n/4
   real(8),dimension(n) :: Ek, Ep, Et
   real(8),dimension(0:ncor) :: Ek_acf, Ep_acf, Et_acf
   real(8),dimension(0:ncor) :: Ek_cor_avg, Ep_cor_avg, Et_cor_avg
   real(8),dimension(nsample) :: Ek_avg, Ep_avg, Et_avg
   real(8) :: Ek2_avg, Ep2_avg, Ek2_sum,Ep2_sum, Et2_avg, Et2_sum
   character(len=4) :: jstr
   integer :: tmp
   Ek_cor_avg = 0.0d0
   Ep_cor_avg = 0.0d0
   Ek_avg = 0.0d0
   Ep_avg = 0.0d0
!j=1
   do j = 1, nsample
      write(jstr, '(I0)') j
      open(11, file="traj_"//trim(jstr))
!open(11, file="traj_1")

      do i = 1, n

         read(11, *) tmp, Ep(i), Ek(i)
	 Et(i) = Ep(i) + Ek(i)
         Ek_avg(j) = Ek_avg(j) + Ek(i)
         Ep_avg(j) = Ep_avg(j) + Ep(i)
	 Et_avg(j) = Et_avg(j) + Ek(i) + Ep(i)

      end do
      close(11)

      Ek_avg(j) = Ek_avg(j) /dble(n)
      Ep_avg(j) = Ep_avg(j) /dble(n)
      Et_avg(j) = Et_avg(j) /dble(n)

      Ek = Ek - Ek_avg(j)
      Ep = Ep - Ep_avg(j)
      Et = Et - Et_avg(j)

      call acf(Ek_acf, Ek, n, ncor)
      call acf(Ep_acf, Ep, n, ncor)
      call acf(Et_acf, Et, n, ncor)

      Ek2_avg = 0.0d0
      Ep2_avg = 0.0d0
      Et2_avg = 0.0d0

      do i = 1, n
         Ek2_avg = Ek2_avg + Ek(i)*Ek(i)
         Ep2_avg = Ep2_avg + Ep(i)*Ep(i)
         Et2_avg = Et2_avg + Et(i)*Et(i)
      end do


      Ek2_avg = Ek2_avg / dble(n)
      Ep2_avg = Ep2_avg / dble(n)
      Et2_avg = Et2_avg / dble(n)

      Ek_acf = Ek_acf / Ek2_avg
      Ep_acf = Ep_acf / Ep2_avg
      Et_acf = Et_acf / Et2_avg

      Ek_cor_avg = (Ek_cor_avg*dble(j-1)+Ek_acf)/dble(j)
      Ep_cor_avg = (Ep_cor_avg*dble(j-1)+Ep_acf)/dble(j)
      Et_cor_avg = (Et_cor_avg*dble(j-1)+Et_acf)/dble(j)
   end do

   open(12, file="mdcor")
   do i = 1, ncor
	Ek_cor_avg(i)=Ek_cor_avg(i-1)+Ek_cor_avg(i)
	Ep_cor_avg(i)=Ep_cor_avg(i-1)+Ep_cor_avg(i)
	Et_cor_avg(i)=Et_cor_avg(i-1)+Et_cor_avg(i)
      write(12, '(I8,F16.8,F16.8,F16.8)') i, Ek_cor_avg(i)/Ek_cor_avg(0), Ep_cor_avg(i)/Ep_cor_avg(0), Et_cor_avg(i)/Et_cor_avg(0)
   end do
!write(12,'(F16.8,F16.8)') Ek_cor_avg(ncor)/Ek_cor_avg(0),sum(Ep_cor_avg)/Ep_cor_avg(0)
   close(12)

end program main

!----------------------------------------------------------------------------------
!> calculate auto correlation function of observable A
!!
!! @see [AT90] Allen and Tildesley, "Computer Simulation of Liquids", P.186
!!
!! @param [in]  A(1:iTtotal)  = values of observable A
!! @param [in]  iTtotal       = size of array A
!! @param [in]  iTcorr        = size-1 of correlation function
!! @param [out] acf(0:iTcorr) = normalized auto-correlation function
subroutine autocorr(acf, A, iTtotal, iTcorr)
   implicit none
   integer, intent(in)  ::  iTtotal       ! total number of A
   integer, intent(in)  ::  iTcorr        ! correlation length
   real(8), intent(in)  ::  A(1:iTtotal)  ! array of values of A
   real(8), intent(out) ::  acf(0:iTcorr) ! normalized auto-correlation function

   real(8) ::  nACF(0:iTcorr)             ! number of ACF, for normalization of ACF
   integer :: iT, iT0, iT0max, iTau       ! time index of t0+tau, t0, t0max and tau

   ! initialize all the auto-correlation functions to be calculated
   acf(0:iTcorr)  = 0d0
   nACF(0:iTcorr) = 0d0

   ! for each t0 from 1 to total, calculate auto-correlation function
   do iT0 = 1, iTtotal

      ! the max index of T0 must not exceed the last index
      iT0max = min(iTtotal, iT0+iTcorr)

      ! for all possible t0+tau, calculate a(t0)*a(t0+tau) and add them
      do iT = iT0, iT0max

         ! interval between t and t0
         iTau  = iT - iT0

         ! calculate auto-correlation function for current interval
         acf(iTau) = acf(iTau) + A(iT0) * A(iT)

         ! count the number of auto-correlation function
         nACF(iTau) = nACF(iTau) + 1d0
      end do

   end do

   ! normalize auto-correlation function
   do iT = 0, iTcorr
      acf(iT) = acf(iT) / nACF(iT)
   end do

   return
end subroutine autocorr

!______________________________________________________________________
!
! Use subroutines in FFTW3 library to do fft calculation
!
subroutine fft_r2c(cfout, rfin, n)
   implicit none
   integer :: n
   real(8) :: rfin(n)
   complex(8) :: cfout(n/2+1)

   ! variables for calling fftw
   integer(8) :: plan_forward

   !
   ! Make sure you have this file in local directory, or include directory
   !
   include "fftw3.f"

   !
   !  Set up a plan, and execute the plan to transform the IN data to
   !  the OUT FFT coefficients.
   !
   call dfftw_plan_dft_r2c_1d_ ( plan_forward, n, rfin, cfout, FFTW_ESTIMATE )

   call dfftw_execute_ ( plan_forward )

   !
   !  Discard the information associated with the plans.
   !
   call dfftw_destroy_plan_ ( plan_forward )

end subroutine fft_r2c

!______________________________________________________________________
!
! Use subroutines in FFTW3 library to do fft calculation
!
subroutine fft_c2r(rfout, cfin, n)
   implicit none
   integer :: n
   complex(8) :: cfin(n/2+1)
   real(8) :: rfout(n)

   ! variables for calling fftw
   integer(8) :: plan_backward

   !
   ! Make sure you have this file in local directory, or include directory
   !
   include "fftw3.f"

   !
   !  Set up a plan, and execute the plan to backtransform the
   !  complex FFT coefficients in OUT to real data.
   !
   call dfftw_plan_dft_c2r_1d_ ( plan_backward, n, cfin, rfout, FFTW_ESTIMATE )

   call dfftw_execute_ ( plan_backward )

   !
   !  Discard the information associated with the plans.
   !
   call dfftw_destroy_plan_ ( plan_backward )

end subroutine fft_c2r
!!
!!
subroutine acf(C, A, n, ncor)
   implicit none
   integer :: n, ncor
   real(8) :: A(n), C(0:ncor)
   integer(8) :: plan_forward, plan_backward
   integer :: i, nfft
   real(8), allocatable, dimension(:) :: rfin, rfout
   complex(8),allocatable,dimension(:) :: cfin, cfout

   include "fftw3.f"

   !! size for FFT
   nfft = 2**(ceiling(log(2.0*n-1.0)/log(2.0)))

   allocate(rfin(nfft))
   allocate(rfout(nfft))
   allocate(cfin(nfft/2+1))
   allocate(cfout(nfft/2+1))

   !! zero padding
   rfin = 0.0d0
   rfin(1:n) = A

   !! FFT   
   call fft_r2c(cfout, rfin, nfft)

   !! modulas
   cfin = cfout*conjg(cfout)/dble(nfft)

   !! Inverse FFT
   call fft_c2r(rfout, cfin, nfft)

   !! normalization
   do i = 0, ncor
      c(i) = rfout(i+1) / dble(n-i) 
   end do

   deallocate(rfin)
   deallocate(rfout)
   deallocate(cfin)
   deallocate(cfout)

end subroutine
