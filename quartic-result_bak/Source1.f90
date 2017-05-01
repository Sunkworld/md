module global
    implicit none
    integer,parameter :: length=20
    real(8) :: proportion=0.000001d0
	real(8) :: datapot(length),dataham(length)
    integer :: boolepot,booleham!,boolepoh
    common boolepot,booleham!,boolepoh
    end module
    
    subroutine select(nonumb,ddata,boolepoh)!pot,ddataham)
    use global
	implicit none
	integer :: nonumb
	real(8) :: ddata!pot,ddataham
    integer :: boolepoh
	real(8) :: stda
    real(8) :: average
	integer :: i

	i=mod(nonumb,length)+1
    if (boolepoh==0) then
        datapot(i)=ddata
    else
        dataham(i)=ddata
    end if
    
	if (nonumb<length) then
        return
    else
            if (boolepoh==0) then
                average=sum(datapot(:))/length
                stda=sqrt(sum((datapot(:)-average)**2)/(length-1)/length)
                if ((stda/average)<proportion) then
                    boolepot=1
                end if
            else
                average=sum(dataham(:))/length
                stda=sqrt(sum((dataham(:)-average)**2)/(length-1)/length)
                if ((stda/average)<proportion) then
                    booleham=1
                end if
            end if
        return
        end if
    end subroutine
	
program main
    use global
    implicit none
	real(8) :: cork,corp,cort,tstep!,maxi
    real :: gamma1(12)
    character(len=4) :: dt(3)
    integer :: gamma2(5)
    character(len=5) :: gammalist(21)
	integer :: n,nonu
    integer :: m,l,p
	
    data dt /"0.05","0.1","0.3"/
    data gamma1 /0.1,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,4.0/
    data gamma2 /10,20,40,50,100/
    do m=1,12
        write(gammalist(m),'(f3.1)') gamma1(m)
    end do
    do m=13,17
        write(gammalist(m),'(i5)') gamma2(m-12)
    end do
    
    n=500000
    do l=1,3
    open(666,file=dt(l)//'pot.csv',position='append')
    open(333,file=dt(l)//'ham.csv',position='append')
        do p=1,17
            boolepot=0
            booleham=0
	        open(999,file=trim(adjustl(dt(l)))//'-'//trim(adjustl(gammalist(p)))//'-mdcor')
!maxi=0.0
!    if (l>n) then
!        write(*,*) "l<n"
!        stop
!    end if

	        do m=1,n
    		    read(999, '(I8,F16.8,F16.8,F16.8)') nonu,cork,corp,cort
!               write(*,*) nonu,cork,corp,cort
!		        if (corp>maxi) maxi=corp
                read(dt(l),*) tstep
                if (boolepot==0) then
                    call select(m,corp,0)
                    if (boolepot==1) write(666,'(a5,a1,i8,a1,f16.8)') trim(adjustl(gammalist(p))),',',m,',',sum(datapot(:))/length*tstep
                end if
                if (booleham==0) then
                    call select(m,cort,1)
                    if (booleham==1) write(333,'(a5,a1,i8,a1,f16.8)') trim(adjustl(gammalist(p))),',',m,',',sum(dataham(:))/length*tstep
                end if
                if ((boolepot==1) .and. (booleham==1)) exit
            end do
            close(999)
        end do
    close(333)
    close(666)
    end do
    end program
    
	
