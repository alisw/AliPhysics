*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 15:37:35  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE EVOLIO (NVLIN,NAMIN,NUMIN,NVLOUT,NAMOUT,NUMOUT)
*
C.    ******************************************************************
C.    *                                                                *
C.    *    Called with INWVOL=2 this routine save the volume's tree    *
C.    *                                                                *
C.    *    Called with INWVOL=1 this routine returns                   *
C.    *             - the list of volumes (and number) left            *
C.    *             - the list of volumes (and number) entered         *
C.    *                                                                *
C.    *    ==>Called by : ERTRGO                                       *
C.    *       Author    M.Maire  *********                             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcvolu.inc"
#include "geant321/gctrak.inc"
*
      PARAMETER (NVMX=15)
      COMMON/OLDTRE/NLVOLD,NAMOLD(NVMX),NUMOLD(NVMX),IOFLAG(NVMX),KONLY
*
      DIMENSION NAMIN(NVMX),NUMIN(NVMX),NAMOUT(NVMX),NUMOUT(NVMX)
*
*
      NVLOUT = 0
      NVLIN  = 0
      IF(INWVOL.EQ.0) RETURN
      IF(NLEVEL.EQ.0) RETURN
*
* *** Save the tree when leaving a volume
*
      IF(INWVOL.EQ.2) THEN
         NLVOLD = NLEVEL
         KONLY  = 1
         DO 1 ILV =1,NLEVEL
           NAMOLD(ILV) = NAMES (ILV)
           NUMOLD(ILV) = NUMBER(ILV)
           IOFLAG(ILV) = 2
           IF(GONLY(ILV).NE.1.) KONLY = 0
    1    CONTINUE
         RETURN
      ENDIF
*
* *** Entering in a new volume : compare the trees
*
      IF(INWVOL.NE.1) RETURN
*
* *** Simplest case : old tree is 'ONLY'
*
      IF(KONLY.EQ.1) THEN
        NLMIN = MIN(NLEVEL,NLVOLD)
        NL    = 0
   10   NL    = NL + 1
        IF(NL.GT.NLMIN) GO TO 12
        IF((NAMES (NL).EQ.NAMOLD(NL)).AND.
     *     (NUMBER(NL).EQ.NUMOLD(NL))) GO TO 10
*
   12   CONTINUE
        IF(NL.LE.NLVOLD) THEN
          DO 14 IL = NL,NLVOLD
            NVLOUT = NVLOUT + 1
            NAMOUT(NVLOUT) = NAMOLD(IL)
            NUMOUT(NVLOUT) = NUMOLD(IL)
   14     CONTINUE
        ENDIF
*
        IF(NL.LE.NLEVEL) THEN
          DO 16 IL = NL,NLEVEL
            NVLIN  = NVLIN  + 1
            NAMIN (NVLIN ) = NAMES (IL)
            NUMIN (NVLIN ) = NUMBER(IL)
   16     CONTINUE
        ENDIF
*
        RETURN
      ENDIF
*
* *** Special case for 'NOT ONLY' volumes
*
      IF(KONLY.EQ.0) THEN
         DO 24 IL1 = 1,NLEVEL
           DO 22 IL2 = 1,NLVOLD
              IF((NAMES (IL1).EQ.NAMOLD(IL2)).AND.
     *           (NUMBER(IL1).EQ.NUMOLD(IL2))) THEN
                   IOFLAG(IL2) = 0
                   GO TO 24
              ENDIF
   22      CONTINUE
           NVLIN  = NVLIN  + 1
           NAMIN (NVLIN ) = NAMES (IL1)
           NUMIN (NVLIN ) = NUMBER(IL1)
   24    CONTINUE
*
         DO 26 IL2 = 1,NLVOLD
           IF(IOFLAG(IL2).EQ.0) GO TO 26
           NVLOUT = NVLOUT + 1
           NAMOUT(NVLOUT) = NAMOLD(IL2)
           NUMOUT(NVLOUT) = NUMOLD(IL2)
   26    CONTINUE
*
      ENDIF
*
      END
