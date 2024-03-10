
integer function pivot(n,a,lo,hi)
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: lo
        integer, intent(in) :: hi
        integer, intent(in) :: a(n)
        pivot=(lo+hi)/2
end function


subroutine swap(a,b)
        implicit none
        integer, intent(inout) :: a,b
        integer::tmp
        tmp=b
        b=a
        a=tmp
end subroutine


integer function partitionArray(n,a,lo,hi) result(storeIndex)
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: lo
        integer, intent(in) :: hi
        integer, intent(inout) :: a(n)

        ! locals
        integer :: i,pivotIndex,pivotValue

        interface
        integer function pivot(n,a,lo,hi)
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: lo
        integer, intent(in) :: hi
        integer, intent(in), dimension(n)::a
        end function
        end interface

     
        pivotIndex=pivot(n,a,lo,hi)
        pivotValue=a(pivotIndex)

        call swap(a(hi),a(pivotIndex))

        storeIndex=lo

        do i=lo,hi-1
           if ( a(i) < pivotValue) then
              call swap(a(i),a(storeIndex))
              storeIndex=storeIndex+1
           end if
        end do
        call swap(a(storeIndex),a(hi))
end function


recursive subroutine seq_qsort(n,a,lo,hi)
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: lo
        integer, intent(in) :: hi
        integer, intent(inout) ::a(n)

        integer :: p

        interface
        integer function partitionArray(n,a,lo,hi) result(storeIndex)
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: lo
        integer, intent(in) :: hi
        integer, intent(inout) ::a(n)
        end function
        end interface

        if ( lo < hi ) then
             p=partitionArray(n,a,lo,hi)
             call seq_qsort(n,a,lo,p-1)
             call seq_qsort(n,a,p+1,hi)
        end if
end subroutine
        
recursive subroutine omp_qsort(n,a,lo,hi)
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: lo
        integer, intent(in) :: hi
        integer, intent(inout) :: a(n)

        integer :: p
        integer,parameter :: maxseq=10

        interface
        integer function partitionArray(n,a,lo,hi) result(storeIndex)
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: lo
        integer, intent(in) :: hi
        integer, intent(inout) ::a(n)
        end function
        end interface

        if (  hi-lo+1 < maxseq ) then
            call seq_qsort(n,a,lo,hi)
        else
             p=partitionArray(n,a,lo,hi)
             !$omp task default(none) shared(n,a) firstprivate(lo,p)
             call omp_qsort(n,a,lo,p-1)
             !$omp end task 
             !$omp task default(none) shared(n,a) firstprivate(p,hi)
             call omp_qsort(n,a,p+1,hi)
             !$omp end task
        end if
end subroutine

program qsort_test
use omp_lib
implicit none
integer, parameter   :: N=1000000
integer,dimension(N) :: array
real(8) :: wtime,tbeg,tend
integer :: i,nthr

do i=1,N
   array(i) = 2*N+1-2*i
end do

print *, array(1:10)

wtime=omp_get_wtime()
!$omp parallel shared(array) 
!$omp single
call omp_qsort(N,array,1,N)
!$omp end single nowait
!$omp end parallel
wtime=omp_get_wtime()-wtime

print * 
print *, array(1:10)
print *, 'parallel = ', wtime

do i=1,N
   array(i) = 2*N+1-2*i
end do
call cpu_time(tbeg)
call seq_qsort(N,array,1,N)
call cpu_time(tend)
print *,'serial = ',tend-tbeg

end program

