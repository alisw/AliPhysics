*CMZ :          26/02/96  11.38.47  by  S.Ravndal
*-- Author :
#if defined(CERNLIB_PARA)
      subroutine checkstr( chREADIN, chExpected, nch, chDescription,
     $ierr )
      implicit none
      character chREADIN*(*), chExpected*(*), chDescription*(*)
      integer   ierr, nch

#include "geant321/gcunit.inc"
      integer   lunerr
      lunerr= lout

*    Example:
*     call checkstr( chSUBSEQ, 'SUBSEQ ', !!                          ! '(a7)'
*    $     ,'Number  of  subsequences ', ichkerr(1) )

      if ( chReadin(1:nch) .ne. chExpected(1:nch) ) then
         write( lunerr, * )
     $        ' Checkstr: Values of string read-in and expected differ '
         write( lunerr, * ) '   Expected>', chExpected , '<'
         write( lunerr, * ) '   Read-in >', chReadin   , '<'
         ierr = 1
      else
         ierr = 0
      endif
      return
      end
#else
      SUBROUTINE CHECKSTR_DUMMY
      END
#endif
