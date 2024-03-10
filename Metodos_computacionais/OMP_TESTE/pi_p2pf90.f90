!**********************************************************************
!   pi3f90.f - compute pi by integrating f(x) = 4/(1 + x**2)     
!
!  (C) 2001 by Argonne National Laboratory.
!      See COPYRIGHT in top-level directory.
!     
!   Each node: 
!    1) receives the number of rectangles used in the approximation.
!    2) calculates the areas of it's rectangles.
!    3) Synchronizes for a global summation.
!   Node 0 prints the result.
!
!  Variables:
!
!    pi  the calculated result
!    n   number of points of integration.  
!    x           midpoint of each rectangle's interval
!    f           function to integrate
!    sum,pi      area of rectangles
!    tmp         temporary scratch space for global summation
!    i           do loop index
!****************************************************************************
program main

 use mpi

 implicit none

 double precision  PI25DT
 parameter        (PI25DT = 3.141592653589793238462643d0)

 double precision  mypi, pi, pibuff, h, sum, x, f, a,time
 integer n, myid, numprocs, i, rc,tag,srcid,rcvid, ierr, sizetype, sumtype
 integer :: mpistatus(MPI_STATUS_SIZE)


!                                 function to integrate
 f(a) = 4.d0 / (1.d0 + a*a)
 
 call MPI_INIT( ierr )
 call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
 call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
 print *, 'Process ', myid, ' of ', numprocs, ' is alive'
 
 sizetype   = 1
 sumtype    = 2
 
 !do 
    if ( myid .eq. 0 ) then
       write(6,98)
 98    format('Enter the number of intervals: (0 quits)')
       read(5,99) n
 99    format(i10)
    endif


    time=MPI_WTIME()
      
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!                                 check for quit signal
    if ( n .le. 0 ) stop

!                                 calculate the interval size
    h = 1.0d0/n
 
    sum  = 0.0d0
    do i = myid+1, n, numprocs
       x = h * (dble(i) - 0.5d0)
       sum = sum + f(x)
    enddo
    mypi = h * sum

    print *, myid,mypi

    if ( myid == 0 ) then
        pi=mypi
        do srcid=1,numprocs-1
           tag=srcid
           print *,'Rcvid : ', srcid,pi,pibuff
           call MPI_Recv(pibuff,1,MPI_REAL8, &
                         srcid,tag,MPI_COMM_WORLD,mpistatus,ierr);
           print *, srcid, pibuff
           pi = pi + pibuff
        end do 
        write(6, 97) pi, abs(pi - PI25DT)
 97     format('  pi is approximately: ', F18.16, &
               '  Error is: ', F18.16)
    else
        tag=myid
        rcvid=0
        print *, 'MPI send : ', myid, mypi
        call MPI_Send(mypi,1,MPI_REAL8, &
                  rcvid,tag,MPI_COMM_WORLD,ierr);
    end if
!                                 node 0 prints the answer.
 !enddo

 time=MPI_WTIME()-time
 print *,'MPI time : ',time
 call MPI_FINALIZE(rc)
 stop
end
