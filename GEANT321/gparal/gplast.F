*CMZ :          23/08/94  15.10.00  by  John Apostolakis
*-- Author :    John Apostolakis  CERN/GP-MIMD2   23/08/94
#if defined(CERNLIB_PARA)
      subroutine gplast
c
c       This routine cleans up all the parallel environment on which
c      parallel Geant is layered. No parallel Geant routines should
c      be called after calling gplast.
c       To ensure that compatibility with all current MPI implementations
c      we ask that no code is executed after calling gplast.
c       To implement this choice we call 'STOP' at the end !!
c
      implicit none
      integer  ierror

#include "geant321/multiprox.inc"
* +SEQ, IOUNIT.    to get LUNOUT/LUNERR --- the keep is defined in geanf

      call MPI_finalize( ierror )
      if (ierror .ne. 0)  then
	 write(*,*) ' gplast: Error in MPI_finalize on node ', nprank
*	 write(LUNERR,*) ' gplast: Error in MPI_finalize on node ', nprank
      endif

      write(*,*) ' Calling STOP in MPI_finalize on node ', nprank
      stop
      end
#endif
