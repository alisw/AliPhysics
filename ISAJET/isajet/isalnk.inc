#ifndef CERNLIB_ISAJET_ISALNK_INC
#define CERNLIB_ISAJET_ISALNK_INC
*
*
* isalnk.inc
*
      INTEGER NVD,MQREF,MPQREF
      PARAMETER (NVD=100)
      PARAMETER (MQREF=200)
      PARAMETER (MPQREF=10)
      COMMON/ISALNK/LVD(NVD),QREF(MQREF),PQREF(MPQREF)
      SAVE /ISALNK/
      INTEGER LVD          ! vertex links
      INTEGER QREF         ! initial and final parton links
      INTEGER PQREF        ! links to primary jets

#endif
