module init
  implicit none
  real(8), parameter :: kT = 1d0
  real(8), parameter :: gamma = 10d0
  real(8), parameter :: h = 0.75d0
  integer :: Nstep = 100000000
  real(8), parameter :: k = 1d0
  real(8), parameter :: q0 = -1d0
  real(8), parameter :: p0 = -0.1d0
  real(8), parameter :: m = 1d0
end module init

program main

  integer :: d = 2

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
  integer :: i
  write(*,*) kT
  qn = q0
  pn = p0
  call calForce(fn, q0)
  open(2333, file='q_baoab.txt')
  open(6666, file='p_baoab.txt')
  do i=1, Nstep
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
     write(2333,*) qn
     write(6666,*) pn
  end do

end subroutine BAOAB
  
subroutine molphys
  use init
  use random
  implicit none
  real(8) :: rand, qn, pn, qnp1, pnp1, fn, fnp1, a, b
  integer :: i
  qn = q0
  pn = p0
  call calForce(fn, q0)
  open(233, file='q_mp.txt')
  open(666, file='p_mp.txt')

  a = (1-gamma*h/2/m)/(1+gamma*h/2/m)
  b = 1/(1+gamma*h/2/m)
  write(*,*) a,b

  do i=1, Nstep
     call random_normal(rand)
     qnp1 = qn + b*h*pn/m + b*h**2/2/m*fn + b*h/2/m*rand*sqrt(2*gamma*kT*h)
     call calForce(fnp1, qnp1)
     pnp1 = a*pn/m + h/2/m*(a*fn+fnp1) + b/m*rand*sqrt(2*gamma*kT*h)
     qn = qnp1
     pn = pnp1
     fn = fnp1
     write(233, *) qn
     write(666, *) pn

  end do



end subroutine molphys

subroutine vv
  use init
  use random
  implicit none
  real(8) :: rand, qn, pn, qnp1, pnp1, fn, fnp1, c1, c2
  integer :: i
  qn = q0
  pn = p0
  call calForce(fn, q0)
  open(234, file='q_vv.txt')
  open(667, file='p_vv.txt')
  c1 = exp(-gamma*h/2)
  c2 = sqrt((1-c1**2))
  do i=1, Nstep
     call random_normal(rand)
     pn = c1*pn + c2*sqrt(m)*rand
     pn = pn + fn*h/2
     qn = qn + pn/m*h
     call calForce(fn, qn)
     pn = pn + fn*h/2
     call random_normal(rand)
     pn = c1*pn + c2*sqrt(m)*rand
     write(234,*) qn
     write(667,*) pn

  end do
  
end subroutine vv

