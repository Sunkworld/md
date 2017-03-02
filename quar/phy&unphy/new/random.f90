!----------------------------------------------------------------------
!> @brief this module contains PRNG (pseudo random number generator)
!!        and methods for calculating mean value and standard deviation
!!
!!    No need external modules or subroutines
!----------------------------------------------------------------------
module random
contains
   !-------------------------------------------------------------------
   !> @brief generate random number from Normal(Gauss) distribution
   !!        with mean value (default 0) and standard deviation (default 1).
   !!        The Box-Muller method is used, and the cosine and sine operations 
   !!        are avoided using the algorithm in $7.3.4 of [NR3rd]
   !!
   !! @param [out] harvest =  the generated random number
   !!        [in][optional] mean    = mean value
   !!        [in][optional] std_dev = standard deviation
   !!
   !! @see   [NR3rd] W. Press, S. Teukolsky, W. Vetterling, B. Flannery, 
   !!        "Numerical Recipes, The Art of Scientific Computing", 3rd Ed.
   !-------------------------------------------------------------------
   subroutine Random_Normal(harvest, mean, std_dev)
      implicit none
      real(8), intent(out) :: harvest
      real(8), intent(in), optional :: mean, std_dev

      real(8) :: v1, v2, rsq
      real(8) :: fac, s, m
      real(8), save :: stored_val
      logical, save :: stored = .false.      

      if (present(std_dev)) then
         s = std_dev
      else
         s = 1d0   ! default standard deviation value is 1
      end if

      if (present(mean)) then
         m = mean
      else
         m = 0d0    ! default mean value is 0
      end if

      if (stored) then
         harvest = stored_val    ! return the stored value
         stored = .false.        ! now there is no value stored
      else
         rsq = 0d0
         ! check if generated values in the unit circle
         do while ((rsq .le. 0d0).or.(rsq .gt. 1d0)) 

            ! generate two random numbers in uniform distribution
            call my_random_number(v1)  
            call my_random_number(v2)  

            ! extend them from -1 to 1
            v1 = 2d0*v1 - 1d0          
            v2 = 2d0*v2 - 1d0          

            ! calculate square of the radius
            rsq = v1*v1 + v2*v2        
         enddo 
         fac = Sqrt(-2d0*Log(rsq)/rsq) ! Box-Muller transformation
         harvest = v1 * fac * s + m    ! return one value,
         stored_val = v2 * fac * s + m ! store the other
         stored = .true.
      end if

      return
   end subroutine Random_Normal

   !-------------------------------------------------------------------
   !> @brief generate a random number in uniformed distribution (0, 1)
   !!        using intrinsic function 'random_number' with a seed check 
   !!
   !! @param [out] val = random number generated by intrinsic function
   !-------------------------------------------------------------------
   subroutine My_Random_Number(val)
      implicit none
      real(8) :: val

      logical, save :: bSeed = .false.

      if (.not. bSeed) then
         call init_random_seed()
         bSeed = .true.
      end if

      call Random_Number(val)

      return
   end subroutine My_Random_Number

   !---------------------------------------------------------------------
   ! initialize the seed for random number generator
   ! @see [1] http://fortranwiki.org/fortran/show/random_seed 
   !---------------------------------------------------------------------
   subroutine init_random_seed()
      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)

      deallocate(seed)
      return
   end subroutine init_random_seed
end module random
