*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:03:23  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE DVNBKD
C  *******************************************************************
C  *                                                                 *
C  *                                                                 *
C  *     -   D I V O N N E 4  -  PROGRAM FOR MULTIPLE INTEGRATION    *
C  *                                           AND                   *
C  *                               ADAPTIVE IMPORTANCE SAMPLING      *
C  *                                            BY                   *
C  *                                   NESTED PARTITIONING.          *
C  *                                                                 *
C  *                                                                 *
C  *                                                                 *
C  *                                                                 *
C  *     CODED BY:     J. H. FRIEDMAN, SLAC/CERN.                    *
C  *                                  AND                            *
C  *                   M. H. WRIGHT, STANFORD UNIVERSITY             *
C  *     MODIFIED FOR CERN BY T.LINDELOF, JULY 1981                  *
C  *     FORTRAN 77 UPDATE BY T.LINDELOF, MAY 1982                   *
C  *******************************************************************
C
C   INITIALIZATION ROUTINE TO 'SIMULATE' BLOCK DATA
C     BLOCK DATA
#include "d151dt.inc"
      COMMON /PRINT/ IPRINT
      COMMON /ISTRGE/ MXRGNS , ISTOR(12000)
      COMMON /RSTRGE/ RSTSZE,RSTOR(18001)
      INTEGER RSTSZE
      COMMON /QUADRE/ IDEG
      COMMON /START/ ISTART
      COMMON /EXFILE/ NFILE
      COMMON /DISPOS/ IDISP
      COMMON /DEPTHS/ ISTDPH , INCDPH
      COMMON /SAMPLE/ NPOINT
      COMMON /CUTOLS/ BNDTOL, FRACT, REGNTL, FNLTOL
      COMMON /BNDLMT/ FLOBD,FUPBD
      COMMON /PRSTOP/ NSTOP
      COMMON /ZEETRM/ ITRMF
C
      DATA INITL/0/
      SAVE INITL
      IF(INITL.NE.0) RETURN
      INITL=1
C
      IDATE='17/12/80'
      IPRINT=1
      MXRGNS=3000
      RSTSZE=18001
      IDEG=0
      ISTART=1
      NFILE=1
      IDISP=0
      ISTDPH=3
      INCDPH=5
      NPOINT=50
      BNDTOL=.05
      FRACT=.5
      REGNTL=.9
      FNLTOL=.1
#if defined(CERNLIB_IBM)||defined(CERNLIB_SINGLE)
      FLOBD=-9.9E71
      FUPBD=9.9E71
#endif
#if (defined(CERNLIB_DOUBLE))&&(!defined(CERNLIB_IBM))
      FLOBD=-9.9E37
      FUPBD= 9.9E37
#endif
      NSTOP=5
      ITRMF=0
C
      CALL DVNOPT
C
      RETURN
      END
