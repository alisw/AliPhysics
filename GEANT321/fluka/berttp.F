*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 22/02/95  12.23.11  by  S.Ravndal
*-- Author :
*=== berttp ===========================================================*
*                                                                      *
      SUBROUTINE BERTTP
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C---------------------------------------------------------------------
C SUBNAME = BERTTP --- READ BERTINI DATA
C---------------------------------------------------------------------
C     --------------------------------- EVAPORATION DATA
#include "geant321/eva0.inc"
#include "geant321/hettp.inc"
#include "geant321/inpflg.inc"
#include "geant321/isotop.inc"
#include "geant321/nucgeo.inc"
#include "geant321/nuclev.inc"
#include "geant321/parevt.inc"
#include "geant321/xsepar.inc"
      LOGICAL OPENED,EXISTS
      LOGICAL LRMSCH, LRD1O2, LTRASP
      CHARACTER*100 FILNAM
#if defined(CERNLIB_UNIX)||defined(CERNLIB_VAX)||defined(CERNLIB_CRAY)
      CHARACTER*100 CHROOT
#endif
C---------------------------------------------------------------------
#if defined(CERNLIB_FDEBUG)
      WRITE( LUNOUT,'(A,I2)')
     & ' -/BERTTP(I): EVAP   DATA READ FROM UNIT ', NBERTP
#endif
      INQUIRE(UNIT=NBERTP, OPENED=OPENED)
      IF(OPENED) THEN
         REWIND NBERTP
      ELSE
#if defined(CERNLIB_UNIX)||defined(CERNLIB_CRAY)
         CHROOT=' '
         CALL GETENVF('CERN_ROOT',CHROOT)
         LNROOT = LNBLNK(CHROOT)
         IF(LNROOT.LE.0) THEN
            FILNAM='flukaaf.dat'
         ELSE
            FILNAM=CHROOT(1:LNROOT)//'/lib/flukaaf.dat'
         ENDIF
         INQUIRE(FILE=FILNAM,EXIST=EXISTS)
         IF(.NOT.EXISTS) THEN
           PRINT*,'**********************************'
           PRINT*,'*        F I F A C E             *'
           PRINT*,'*        -----------             *'
           PRINT*,'*   File FLUKAAF.DAT not found   *'
           PRINT*,'*         Program STOP           *'
           PRINT*,'*   Check CERN_ROOT environment  *'
           PRINT*,'*           variable             *'
           PRINT*,'**********************************'
           STOP
         ENDIF
         OPEN(NBERTP,FILE=FILNAM,STATUS='OLD')
#endif
#if defined(CERNLIB_VAX)
         ISTAT = LIB$SYS_TRNLOG ('CERN_ROOT',NALL,CHROOT,,,%VAL(0))
         IF(ISTAT.NE.1) THEN
            FILNAM='flukaaf.dat'
         ELSE
            FILNAM='CERN_ROOT:[LIB]flukaaf.dat'
         ENDIF
         INQUIRE(FILE=FILNAM,EXIST=EXISTS)
         IF(.NOT.EXISTS) THEN
           PRINT*,'**********************************'
           PRINT*,'*        F I F A C E             *'
           PRINT*,'*        -----------             *'
           PRINT*,'*   File FLUKAAF.DAT not found   *'
           PRINT*,'*         Program STOP           *'
           PRINT*,'*   Check CERN_ROOT environment  *'
           PRINT*,'*           variable             *'
           PRINT*,'**********************************'
           STOP
         ENDIF
         OPEN(NBERTP,FILE=FILNAM,STATUS='OLD',READONLY)
#endif
#if defined(CERNLIB_IBM)
         FILNAM='/FLUKAAF DAT *'
         OPEN(NBERTP,FILE=FILNAM,STATUS='OLD')
#endif
      ENDIF
 
C A. Ferrari: first of all read isotopic data
      READ (NBERTP,2100) ISONDX
      READ (NBERTP,2100) ISOMNM
      READ (NBERTP,2000) ABUISO
      READ (NBERTP,2000) (P0(I),P1(I),P2(I),I=1,1001)
      READ (NBERTP,2100) IA,IZ
      DO 2 I=1,6
         FLA(I)=IA(I)
         FLZ(I)=IZ(I)
    2 CONTINUE
      READ (NBERTP,2000) RHO,OMEGA
      READ (NBERTP,2000) EXMASS
      READ (NBERTP,2000) CAM2
      READ (NBERTP,2000) CAM3
      READ (NBERTP,2000) CAM4
      READ (NBERTP,2000) CAM5
      READ (NBERTP,2000) ((T(I,J),J=1,7),I=1,3)
      DO 3 I=1,7
         T(4,I)=0.D0
    3 CONTINUE
      READ (NBERTP,2000) RMASS
      READ (NBERTP,2000) ALPH
      READ (NBERTP,2000) BET
      READ (NBERTP,2000) WAPS
      READ (NBERTP,2000) APRIME
#if defined(CERNLIB_FDEBUG)
      WRITE( LUNOUT,'(A)' ) ' /DRES(I): USING 1977 WAPS DATA '
#endif
      READ (NBERTP,2200) AHELP , BHELP , LRMSCH, LRD1O2, LTRASP
      IF ( AHELP .NE. ALPHA0 .OR. BHELP .NE. GAMSK0 ) THEN
         WRITE (LUNOUT,*)
     &         ' *** Inconsistent Nuclear Geometry data on file ***'
         STOP
      END IF
      READ (NBERTP,2000) RHOTAB, RHATAB, ALPTAB, RADTAB, SKITAB, HALTAB,
     &              EKATAB, PFATAB, PFRTAB
      READ (NBERTP,2000) AANXSE, BBNXSE, CCNXSE, DDNXSE, EENXSE, ZZNXSE,
     &              EMNXSE, XMNXSE
      READ (NBERTP,2000) AAPXSE, BBPXSE, CCPXSE, DDPXSE, EEPXSE, FFPXSE,
     &              ZZPXSE, EMPXSE, XMPXSE
2000  FORMAT (3(1X,G23.16))
2100  FORMAT (18(1X,I3))
2200  FORMAT (2(1X,G23.16),3(1X,L1))
      CLOSE (UNIT=NBERTP)
      DO 100 JZ = 1, 130
         SHENUC ( JZ, 1 ) = 1.D-03 * ( CAM2 (JZ) + CAM4 (JZ) )
  100 CONTINUE
      DO 200 JA = 1, 200
         SHENUC ( JA, 2 ) = 1.D-03 * ( CAM3 (JA) + CAM5 (JA) )
  200 CONTINUE
      CALL STALIN
      ILVMOD = 1
      IB0 = ILVMOD
#if defined(CERNLIB_FDEBUG)
      WRITE (LUNOUT,*)
      WRITE (LUNOUT,*)' **** Standard EVAP level density used ****'
      WRITE (LUNOUT,*)
     &   ' **** Original Gilbert/Cameron pairing energy used ****'
#endif
      ILVMOD = IB0
      DO 500 JZ = 1, 130
         PAENUC ( JZ, 1 ) = 1.D-03 * CAM4 (JZ)
  500 CONTINUE
      DO 600 JA = 1, 200
         PAENUC ( JA, 2 ) = 1.D-03 * CAM5 (JA)
  600 CONTINUE
      RETURN
      END
