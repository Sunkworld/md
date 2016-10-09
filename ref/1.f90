program main
  character(len=30) :: aa
  real :: i
  open(12,file='222.dat')
  write(12,*) 'sfdsfa'
  write(12,*) 'llll'
  close(12)
  open(12,file='222.dat')
  read(12,*) aa
  write(*,*) aa

end program
