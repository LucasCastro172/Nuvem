module assertion_mod
!
! Assertion utility for error handling in fortran
!
! reference:
! Rouson, D., Xia, J., Xu, X., 2011. Scientific software design: 
! the object oriented way. Cambridge university press.
!
 use iso_fortran_env, only: error_unit
 implicit none
 private
 public:: error_message, assert, assert_identical

 type error_message
      character(:), allocatable :: string
 end type
contains

        subroutine assert(assertion,text)
        logical, dimension(:), intent(in) :: assertion
        type(error_message), dimension(:), intent(in) :: text
        !
        !locals:
        integer::i
        logical:: any_failures
        call assert_identical([size(assertion),size(text)])
        any_failures=.false.
        do i=1,size(assertion)
           if ( .not. assertion(i) ) then
               any_failures=.true.
               write(error_unit,*) 'Assertion failed with message: '
               if ( allocated(text(i)%string) ) then
                     write(error_unit,*) trim(text(i)%string)
               else
                     write(error_unit,*) '(no error message provided)'
               end if
            end if
        end do
        if ( any_failures ) stop 'Execution halted in failed assertion(s)!'
        end subroutine

        subroutine assert_identical(integers)
        integer, dimension(:), intent(in):: integers
        integer:: i
        logical:: any_mismatches
        any_mismatches=.false.
        do i=1,size(integers)
           if ( integers(i) /= integers(1) ) then
              any_mismatches=.true.
              write(error_unit,*) &
                      'Value ',i,'does not match expected value ', integers(1)
           end if
        end do
        if ( any_mismatches ) stop 'Execution halted in failed assertion!'
        end subroutine
end module assertion_mod

