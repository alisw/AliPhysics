*CMZ :           1/08/95  16.58.00  by  John Apostolakis
*-- Author :    John Apostolakis  CERN/GP-MIMD2    1/08/95
**DOC
*
*  Routine to allow calling GEANT program to become aware of
*  the extent of the parallel machine it is running on and its place
*  within it.
*
#if defined(CERNLIB_PARA)
      subroutine gprocs( nsize, nrank, nleader )
      implicit none
      integer nsize, nrank, nleader
*  nsize   = Number of nodes or processes
*  nrank   = The id or rank of the this process ( 0 <= rank < size )
*  nleader = The node that does input distibution (if there is one,
*                  and -1 if there is no leader)
**ENDDOC
#include "geant321/multiprox.inc"
      nsize= npsize
      nrank= nprank
      nleader= npleader

      RETURN
      END
#else
      SUBROUTINE GPROCS_DUMMY
      END
#endif
