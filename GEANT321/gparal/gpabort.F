*CMZ :          23/08/95  15.00.00  by  John Apostolakis
*-- Author :    John Apostolakis  CERN/GP-MIMD2   22/07/94
      subroutine gpabort
#if defined(CERNLIB_PARA)
c
c      This makes a "best attempt" to abort all parallel Geant tasks.
c      It should be used if a fatal error is detected in a task that
c      cannot be recovered, or will cause irreperable damage later.
c
      implicit none
      integer  ierror
#include "geant321/mpifinc.inc"
#include "geant321/multiprox.inc"
#include "geant321/gcunit.inc"
      write(chmail,*)
     $' gpabort: calling MPI_abort to terminate program from  node ',
     $ nprank
      call gmail(2,2)

      call MPI_abort( MPI_COMM_WORLD, ierror )
      if (ierror .ne. 0)  then
	 write(chmail,*)
     $     ' ERROR in gpabort: MPI_abort failed to terminate on node ',
     $     nprank, ' error code: ', ierror, ' STOPPING this node.'
         call gmail(2,2)
         stop
      endif
#endif
      return
      end

