      module mpi_def_m
c
        implicit none
c
        integer :: master, numpe, myrank
	      common /workfc/ master, numpe, myrank
c
      end module mpi_def_m
