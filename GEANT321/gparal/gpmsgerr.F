c-----------------------------------------------------------------
c   Author John Apostolakis, GPMIMD-2, CERN-CN,  February 1994
c
c        gpmsgerr(..)   report a error given by MPI (in sending 
c                       or receiving a message).
c-----------------------------------------------------------------
      
      subroutine gpmsgerr( mpi_routine, ierr, srcroutinename )
      implicit none
      character*(*) srcroutinename, mpi_routine
      integer   lunerr, ierr
#include "geant321/multiprox.inc"
#include "geant321/gcunit.inc"

c     As stderr is not defined, use Geant's stdout
      lunerr= lout

      write(LUNERR,*) 'WARNING: MPI routine ', mpi_routine, 
     &      ' returned error code ', ierr , ' on node ', nprank,
     &      ' when called from routine ', srcroutinename

      RETURN
      END
