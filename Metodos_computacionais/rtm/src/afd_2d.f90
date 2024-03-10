program fwi_main
 use number_types,only:fp,sp,dp
 use rtm_mod
 use mpi
implicit none
type(modelagem) :: fwi
integer         :: ierr

call acquisition_init(fwi)

call acquisition_seismic2(fwi)



call MPI_Finalize(ierr)


end program fwi_main
