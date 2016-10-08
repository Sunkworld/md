program main
  character(len=30) :: aa
  integer :: i

  real(8) :: h = 0.3d0
  write(aa,*) i
  write(*,*)trim(adjustl(aa))//'-e.txt'

    i = 3d6/h
write(*,*) i
end program
