#ifdef __APPLE__
#ifndef __INTEL_COMPILER
#define stop CALL EXIT !
#define STOP CALL EXIT !
#endif
#endif
      SUBROUTINE OPENDECAYTABLE(LUN)
      LOGICAL OPENED,EXISTS
      LOGICAL LRMSCH, LRD1O2, LTRASP
      INTEGER LUN
      CHARACTER*1000 FILNAM
      CHARACTER*1000 CHROOT
      CHROOT=' '
      CALL GETENVF('ALICE_ROOT',CHROOT)
      LNROOT = LNBLNK(CHROOT)
      IF(LNROOT.LE.0) THEN
         FILNAM='decaytable.dat'
      ELSE
         FILNAM=CHROOT(1:LNROOT)//'/data/decaytable.dat'
      ENDIF
      INQUIRE(FILE=FILNAM,EXIST=EXISTS)
      IF(.NOT.EXISTS) THEN
         PRINT*,'***********************************'
         PRINT*,'*        OpenDecayTable           *'
         PRINT*,'*        ---------------          *'
         PRINT*,'*   File DECAYTABLE.DAT not found *'
         PRINT*,'*         Program STOP            *'
         PRINT*,'*   Check ALICE_ROOT environment  *'
         PRINT*,'*           variable              *'
         PRINT*,'***********************************'
         STOP
      ENDIF
#ifdef __GFORTRAN__
      OPEN(LUN,FILE=FILNAM,STATUS='OLD',ACTION='READ')
#else
      OPEN(LUN,FILE=FILNAM,STATUS='OLD')
#endif
      RETURN
      END
      SUBROUTINE CLOSEDECAYTABLE(LUN)
      CLOSE(LUN)
      RETURN
      END
