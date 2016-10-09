program main
  character(len=30) :: aa
  integer :: n =4
  real,allocatable :: i(:,:)
  allocate(i(n,n))
  i(:,:) =1
  write(*,*) i(2,5)

end program
