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

          pi=0.0
          !$omp parallel shared(pi,step) private(i,x,id,sumpi,nthreads)
          id=omp_get_thread_num()+1
          nthreads=omp_get_num_threads()
          sumpi = 0.0
          DO i = id, num_steps,nthreads
             x = (real(i) - 0.5) * step
             sumpi = sumpi + 4.0 / (1.0 + x * x)
          ENDDO
          !$omp critical
             pi=pi+sumpi*step
          !$omp end critical
          !$omp end parallel

          run_time = OMP_GET_WTIME() - start_time

          WRITE(*,100) num_steps, pi, run_time
100       FORMAT('pi with ',i14,' steps is ',f15.8,' in ',f8.3,' secs')

          END PROGRAM MAIN
