*CMZ :          23/03/95  14.24.24  by  John Apostolakis CERN GP-MIMD 2
*-- Author :    John Apostolakis CERN GP-MIMD 2   23/03/95
#if defined(CERNLIB_PARA)
      subroutine pbarrier
      implicit none
c
c     Any processes making this call will wait, until all processes
c     have arrived at this point.
c
#include "geant321/mpifinc.inc"
#include "geant321/gcunit.inc"
      integer ierr

c     It would make sense when debugging to include the following write:
c     write (lout,*)  ' Entering barrier: if any process cannot reach ',
c     ' this point the program will hang here '
c
      call MPI_barrier( MPI_COMM_WORLD, ierr )
      if ( ierr .lt. 0 ) then
          write (lout,*)  ' There was error ', ierr, ' in pbarrier '
      endif
      return
      end
#else
      SUBROUTINE PBARRIER_DUMMY
      END
#endif
