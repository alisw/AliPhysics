#include "isajet/pilot.h"
      SUBROUTINE ISAEVT(I,OK,DONE)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods :
C-
C-         Normal operation:
C-    Generate one ISAJET event and return.
C-
C-         "ISALEP" generation:
C-    Generate a TWOJET or DRELLYAN hard scattering. Then make NEVOLVE 
C-    evolutions and NHADRON fragmentations, rejecting events which 
C-    fail the desired cuts using logical functions
C-         REJJET()   tests the QCD evolution stage, e.g. by requiring
C-                    a heavy quark.
C-         REJFRG()   tests the fragmentation stage, e.g. by requiring
C-                    a high-pt lepton.
C-    These functions default to .FALSE.; i.e. they do not reject any 
C-    events. Note that one hard scattering can give more than one 
C-    event. You must choose NEVOLVE and NHADRON carefully.
C-         IEVT   = event number. This is incremented NEVOLVE * NHADRON 
C-                  times for each hard scattering; i.e. it counts the
C-                  number of potential events.
C-         IEVGEN = counter for generated events.
C-         NEVENT = maximum value of hard scatterings. Hence the limit
C-                  for IEVT is NEVENT * NEVOLVE * NHADRON.
C-    The cross section SIGF contains an extra factor of
C-         1 / (NEVOLVE * NHADRON)
C-    to produce the correct final cross section using the weight
C-         SIGF / NEVENT
C-
C-   Input:
C-    I      = number used to control printout
C-   Output:
C-    OK     = logical flag for good event.
C-    DONE   = logical flag for job completion.
C-
C-   Created   3-FEB-1988   Serban D. Protopopescu
C-   Updated  17-APR-1990 Serban D. Protopopescu (add ISALEP option) 
C-   22-JUL-1992: Move PRTEVT and GETTOT statements to end so they
C-                work for TWOJET and DRELLYAN with NOVOLVE. (FEP)
C-
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/idrun.inc"
#include "isajet/keys.inc"
#include "isajet/nodcay.inc"
#include "isajet/primar.inc"
#include "isajet/jetpar.inc"
#include "isajet/partcl.inc"
#include "isajet/jetset.inc"
#include "isajet/isloop.inc"
C
      LOGICAL REJJET,REJFRG,OK,DONE
      INTEGER NPASS,I,NLIMIT
C
      NPASS=0
      OK=.TRUE.
      DONE=.FALSE.
      NLIMIT=NEVENT*NEVOLV*NFRGMN
C
C          Twojet or Drell-Yan events. The evolution and fragmentation
C          loops are done with GO TO statements so that we can exit
C          the loops with a good event and reenter them.
C
      IF(KEYS(1).OR.KEYS(3)) THEN
100     CONTINUE
        IF(IEVOL.EQ.1.AND.IFRG.EQ.1) THEN
          NPASS=NPASS+1
          IF(NPASS.GT.NTRIES) THEN
            WRITE(ITLIS,1001) NTRIES
1001        FORMAT(//' IT IS TAKING MORE THAN',I6,' TRIES TO MAKE',
     $      ' AN EVENT IN ISAEVT.'/
     $      ' CHECK YOUR LIMITS OR OR INCREASE NTRIES.'/
     $      ' CHECK NEVOLVE, NHADRON, AND YOUR REJJET AND REJFRG',
     $      ' FUNCTIONS IF ANY.'/
     $      ' JOB TERMINATED.')
            STOP 99
          ENDIF
          CALL RANFMT
C          Generate appropriate hard scattering
          IF(KEYS(1)) THEN
            CALL TWOJET
          ELSE
            CALL DRLLYN
          ENDIF
        ENDIF
C          QCD evolution
        IF(NOEVOL) THEN
          IEVT=IEVT+NEVOLV*NFRGMN
          GOTO 9999
        ENDIF
C          Continue if in fragmentation loop
        IF(IFRG.NE.1) GO TO 120
C          Begin multiple evolution loop
110       CONTINUE
          NJSET=0
          IEVT=IEVT+1
          CALL EVOLVE
          IEVT=IEVT-1
          IF(NJSET.LT.0) THEN
            IEVT=IEVT+NFRGMN
            GO TO 111
          ENDIF
          IF(REJJET()) THEN
            IEVT=IEVT+NFRGMN
            GO TO 111
          ENDIF
          IF(NOHADR) THEN
            IEVT=IEVT+NFRGMN
            GO TO 9999
          ENDIF
C          Begin multiple fragmentation loop
120         CONTINUE
            NPTCL=0
            CALL FRGMNT
            IEVT=IEVT+1
            IF(REJFRG()) GO TO 121
C          Finish good event
            CALL MBIAS
            IFRG=IFRG+1
            IF(IFRG.GT.NFRGMN) IFRG=1
            IF(IFRG.EQ.1) THEN
              IEVOL=IEVOL+1
              IF(IEVOL.GT.NEVOLV) IEVOL=1
            ENDIF
            GOTO 9999
C          Fragmentation failed - increment counter and loop
121         IFRG=IFRG+1
            IF(IFRG.GT.NFRGMN) THEN
              IFRG=1
            ELSE
              GO TO 120
            ENDIF
C          End of multiple fragmentation loop
C          Evolution failed - increment counter and loop
111       IEVOL=IEVOL+1
          IF(IEVOL.GT.NEVOLV) THEN
            IEVOL=1
            IFRG=1
            GO TO 100
          ELSE
            GO TO 110
          ENDIF
C
C          E+E- events
C
      ELSE IF(KEYS(2)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL ELCTRN
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(.NOT.NOHADR) CALL FRGMNT
        ENDIF
C
C          MINBIAS events
C
      ELSE IF(KEYS(4)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        NPTCL=0
        IF(.NOT.(NOEVOL.OR.NOHADR)) CALL MBIAS
C
C          SUPERSYM events
C
      ELSE IF(KEYS(5)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL TWOJET
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GO TO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          WPAIR events
C
      ELSE IF(KEYS(6)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL TWOJET
        CALL WPAIR
C
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GO TO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          HIGGS events
C
      ELSE IF(KEYS(7)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL DRLLYN
        CALL HIGGS
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GOTO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          PHOTON events
C
      ELSEIF(KEYS(8)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL TWOJET
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GOTO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          TCOLOR events, e.g. techni-rho
C
      ELSEIF(KEYS(9)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL DRLLYN
        CALL HIGGS
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GOTO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          WHIGGS events
C
      ELSE IF(KEYS(10)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL TWOJET
        CALL WHIGGS
C
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GO TO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          EXTRADIM events
C
      ELSE IF(KEYS(11)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL DRLLYN
C 
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GO TO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
C
C          ZJJ events
C
      ELSEIF(KEYS(12)) THEN
        IEVT=IEVT+1
        CALL RANFMT
        CALL ZJJ
C
        IF(.NOT.NOEVOL) THEN
          CALL EVOLVE
          IF(NJSET.LT.0) GO TO 9999
          IF(.NOT.NOHADR) THEN
            CALL FRGMNT
            CALL MBIAS
          ENDIF
        ENDIF
      ENDIF
C
C          Event complete
C
 9999 IEVGEN=IEVGEN+1
      IF(NJSET.LT.0) OK=.FALSE.
      IF(IEVT.GT.NLIMIT) THEN
        OK=.FALSE.
        DONE=.TRUE.
      ELSEIF(IEVT.EQ.NLIMIT) THEN
        DONE=.TRUE.
      ENDIF
      IF(OK) THEN
        CALL PRTEVT(I)
        CALL GETTOT(.FALSE.)
      ENDIF
      RETURN
      END
