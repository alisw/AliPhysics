*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:10  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GHCLOS
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        SUBR. GHCLOS                                            *
C.    *                                                                *
C.    *   For every sensitive detector, calls GGDETV in order to       *
C.    *   complete the JD bank with the list of volume names which     *
C.    *   permit to identify uniquely a given physical detector, the   *
C.    *   list of bit numbers to pack the corresponding volume copy    *
C.    *   numbers,  and the generic path(s) through the JVOLUM tree.   *
C.    *                                                                *
C.    *   Called by : GGCLOS, <USER>                                   *
C.    *   Author    : F.Bruyant                                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
C.
C.    ------------------------------------------------------------------
*
      IF (JSET.NE.0) THEN
*
* ***   Scan the volume structure to retrieve the path through
*       the physical tree for all sensitive detectors
*
        NSET = IQ(JSET-1)
        DO 20 ISET=1,NSET
          JS = LQ(JSET-ISET)
          NDET = IQ(JS-1)
          DO 10 IDET=1,NDET
            CALL GGDETV (ISET, IDET)
   10     CONTINUE
   20   CONTINUE
      ENDIF
*                                                             END GHCLOS
      END
