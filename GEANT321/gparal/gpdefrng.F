*CMZ :          26/02/96  11.38.47  by  S.Ravndal
*-- Author :
C------------------------------------------------------------------------
c     gpdefrng( ibase )
c
c     Modified: November 6, 1995  J.A. to provide more flexibility by
c                                       including argument
c
c     Function: Provide seeds for the pseudo-random number
c                 generator different on each node
c               If ibase <= 0 default sequences are used.
c
c     Should be called to provide a default seeding.
c     Only currently available for RANECU/GRNDM.
c
c     ibase=  the RNG sequence for first node.
c             If ibase=0 a system default is used.
c
C     called by:  <USER>, uginit, (if CERNLIB_PARA switch is used)
C------------------------------------------------------------------------
#if defined(CERNLIB_PARA) 
      subroutine gpdefrng( ibase )
      implicit none
      integer ibase

C--------------------------------------------------------------
C     Must be called in uginit only after other routines that do the
C   default RNG seeding for sequential machines
C   have been called (ie after ginit) ... see example uginit
C--------------------------------------------------------------
#include "geant321/mpifinc.inc"
#include "geant321/multiprox.inc"
#include "geant321/multseeds.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
      integer ibasedef, lunerr
      data    ibasedef / 1 /
      integer msgerror, rngerror
      data    msgerror, rngerror / 0, 0 /
c----------------------------------------------------------------------------
c     Random number generator seeding:
c
c     Set the seeds on each node for the start of a different sequence.
c----------------------------------------------------------------------------
c     Done for RANECU
c----------------------------------------------------------------------------

      if ( ibase .lt. 0 ) then
          ibase = 0
          write( chmail, '(a,a,i4,a)' )
     $         ' Warning: gpdefrng: ibase < 0 is invalid. Value =',
     $         ibase , ' Reset to zero '
          call gmail(1,1)
      endif
      if ( ibase .eq. 0 ) then
          ibase = ibasedef
      endif

C     We check to see if the resulting sequence is valid for RANECU

      iseqnc = ibase + nprank
      if( (iseqnc .gt.0 ) .and. (iseqnc .le. 215 )) then
          call GRNDMQ( iseeda, iseedb, iseqnc, 'Q' )
          write (*,*)  ' Node ', nprank, ' init-ed seq ', iseqnc,
     &              ' and got seeds ', iseeda, iseedb
          call GRNDMQ( iseeda, iseedb, iseqnc, 'S' )
      else
          write (lunerr, *) ' Node ', nprank, ': in GRNDMQ exceeds ',
     &        ' number of rng sequences with precalculated seeds.'
	  rngerror = 1

c         --- Other Possibility: ---
c         write (lunerr, *) ' Node ', nprank, ':    wrapping around ... '
c         iseqnc = mod( iseqnc, 215 ) + 1
c         call GRNDMQ( iseeda, iseedb, iseqnc, 'Q' )
      endif
#if defined(CERNLIB_PARA_RANMAR)
c----------------------------------------------------------------------------
c     In the case that RANMAR is used
c
c----------------------------------------------------------------------------
      ibase  = nrndm(1)
      iseeda = ibase * 100 + nprank
      iseqnc = 1
c
c     As changed by Stic/Delphi only the first parameter is used.
c
      call GRNDMQ( iseeda, iseedb, iseqnc, 'S' )
c
#endif
      return
      end
#else
      SUBROUTINE GPDEFRNG_DUMMY
      END
#endif

