*CMZ :          20/10/94  15.26.50  by  John Apostolakis CERN GP-MIMD 2
*-- Author :    John Apostolakis CERN GP-MIMD 2   13/07/94
#if defined(CERNLIB_PARA)
      subroutine gptrigi
      implicit none
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "geant321/multiprox.inc"
      integer  iseed1, iseed2

      call GRNDMQ( iseed1, iseed2, 0, 'G' )

      WRITE (CHMAIL,10000) IEVENT,nprank, IDEVT, iseed1, iseed2
      CALL GMAIL(0,0)

10000 FORMAT(' **** GTRIGI: IEVENT=',I7,' NODE ', I3 , ' IDEVT=',I7,
     +' Random Seeds = ',I10,2X,I10)

      return
      end
#endif
