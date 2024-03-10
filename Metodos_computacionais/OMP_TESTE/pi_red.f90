! This program will numerically compute the integral of
!                   4/(1+x*x) 
! from 0 to 1.  The value of this integral is pi -- which 
! is great since it gives us an easy way to check the answer.

! The is the original sequential program.  It uses the timer
! from the OpenMP runtime library

! History: C code written by Tim Mattson, 11/1999.
! adapted to Fortran code by Helen He, 09/2017. 

 
          PROGRAM MAIN
          USE OMP_LIB
          IMPLICIT NONE

          INTEGER :: i,id,nthreads
          INTEGER, PARAMETER :: num_threads = 4
          INTEGER, PARAMETER :: num_steps = 100000000
          INTEGER, PARAMETER :: l1pad = 8
          REAL(8) :: x, pi, step
          REAL(8) :: start_time, run_time
          REAL(8) :: sumpi

          call omp_set_num_threads(num_threads)
          step = 1.0 / num_steps
          start_time = OMP_GET_WTIME()

          sumpi=0.d0
          !$omp parallel shared(sumpi,step) private(i,x)
          !$omp do reduction(+:sumpi)
          DO i = 1, num_steps
             x = (real(i) - 0.5) * step
             sumpi = sumpi + 4.0 / (1.0 + x * x)
          ENDDO
          !$omp end do
          !$omp end parallel
          pi=sumpi*step
          run_time = OMP_GET_WTIME() - start_time

          WRITE(*,100) num_steps, pi, run_time
100       FORMAT('pi with ',i14,' steps is ',f15.8,' in ',f8.3,' secs')

          END PROGRAM MAIN
