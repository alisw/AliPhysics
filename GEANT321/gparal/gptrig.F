*CMZ :          02/11/94  15.22.20  by  John Apostolakis CERN GP-MIMD 2
*-- Author :    John Apostolakis CERN GP-MIMD/2  Started: 13/07/94
C------------------------------------------------------------------------
C     gptrig
C           Function: prints some default information about the number of
C             the event simulated, the node and the current RNG seeds.
C
C     called by:  gtrig (if PARALLEL switch is used)
C------------------------------------------------------------------------
#if defined(CERNLIB_PARA)
      subroutine gptrig
      implicit none
#include "geant321/gcflag.inc"
#include "geant321/multiprox.inc"
#include "geant321/multseeds.inc"
      integer  isequence, idbg

      idbg= 0

c     if (idbg) write( *, 10000)  nprank, iseqnc, iseeda, iseedb
10000 FORMAT (' --> Node',I4,' Before grndmq call. ',
     &        ' Sequence expected ', I5, ' Previous Seeds', 2I11 )

      isequence = -1
      call GRNDMQ( iseeda, iseedb, isequence, 'G' )

      if(idbg) write( *, 10005)  nprank, isequence
10005 FORMAT (' --> Node',I4,' After grndmq call. ',
     &        ' Sequence returned ', I5 )

c     call GRNDMQ( iseeda, iseedb, iseqnc, 'G' )

      write( *, 10010)  nprank, ievent, iseeda, iseedb
10010 FORMAT (' --> Node',I4,'  Event',I7,'  Seeds',2I11 )

      return
      end
#endif
