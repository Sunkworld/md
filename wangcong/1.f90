module global
	implicit none
	integer,parameter length=100
	real proportion=0.01
	real datada(length)
	real :: average
	common average
	

subroutine select(nonumb,maxa,ddata)
	implicit none
	use global
	integer :: nonumb
	real :: maxa,ddata
	real :: stda
	integer :: i,j
	
	i=mod(nonumb,length)
	if (nonumb>length) then
		average=average+(ddata-datada(i))/length
		datada(i)=ddata
		stda=sqrt(sum((datada(:)-average)**2)/(length-1)/length)
		if (stda/maxa<proportion) then
			write(nonumb)
		end if
	else
	datada(i)=ddata
	average=average+
	end if	
	
program main
	real :: arra,nonu,cork,corp,cort,maxi,
	real :: aver
	common aver
	integer,parameter :: n
	integer i
	
	open(99,file='mdcor')
	maxi=0.0
	average=0.0
	do i=1,n
		read(99,*) nonu,cork,corp
		if (corp>maxi) then
		maxi=corp
		call select (100,nonu,0.01,maxi,corp)
	
