*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:17  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GRFILE(LUN,CHFILE,CHOPTT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to open a GEANT/RZ data base.                    *
C.    *                                                                *
C.    *       LUN logical unit number associated to the file           *
C.    *                                                                *
C.    *       CHFILE RZ file name                                      *
C.    *                                                                *
C.    *       CHOPT is a character string which may be                 *
C.    *          'N' To create a new file                              *
C.    *          'U' to open an existing file for update               *
C.    *          ' ' to open an existing file for read only            *
C.    *          'Q' The initial allocation (default 1000 records)     *
C.    *              is given in IQUEST(10)                            *
C.    *          'X' Open the file in exchange format                  *
C.    *          'I' Read all data structures from file to memory      *
C.    *          'O' Write all data structures from memory to file     *
C.    *                                                                *
C.    *    Note:                                                       *
C.    *         If options 'I' or 'O' all data structures are read or  *
C.    *            written from/to file and the file is closed.        *
C.    *         See routine GRMDIR to create subdirectories            *
C.    *         See routines GROUT,GRIN to write,read objects          *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      COMMON/QUEST/IQUEST(100)
      CHARACTER*(*) CHFILE,CHOPTT
      CHARACTER*8 CHOPT,CHTOP,CHTAGS(2)
      CHARACTER*10 CHOPT1,CHOPT2
      SAVE CHTAGS,CHTOP
      DATA CHTAGS,CHTOP/'Object','Version','        '/
C.
C.    ------------------------------------------------------------------
C.
*
      CHOPT=CHOPTT
      CALL CLTOU(CHOPT)
      IOPTN=INDEX(CHOPT,'N')
      IOPTU=INDEX(CHOPT,'U')
      IOPTQ=INDEX(CHOPT,'Q')
      IOPTI=INDEX(CHOPT,'I')
      IOPTO=INDEX(CHOPT,'O')
      IF(IOPTN.EQ.0.AND.IOPTU.EQ.0)IOPTO=0
      IF(IOPTN.NE.0)IOPTI=0
*
*          Open file
*
      CHOPT2 = CHOPT//'XW'
      CHOPT1 = CHOPT//'X'
      CALL RZOPEN(LUN,CHTOP,CHFILE,CHOPT2,1024,ISTAT)
      IF(ISTAT.NE.0)GO TO 99
*
      IF(IOPTN.NE.0)THEN
         IF(IOPTQ.NE.0)THEN
            NQUOT=IQUEST(10)
            IF(NQUOT.LT.100.OR.NQUOT.GT.100000)NQUOT=1000
         ELSE
            NQUOT=1000
         ENDIF
         CALL RZMAKE(LUN,CHTOP,2,'HI',CHTAGS,NQUOT,CHOPT1)
      ELSE
         CALL RZFILE(LUN,CHTOP,CHOPT1)
      ENDIF
*
*           Read/write initialisation structures
*
      IF(IOPTI.NE.0)THEN
         CALL GRIN('INIT',0,' ')
      ENDIF
      IF(IOPTO.NE.0)THEN
         CALL GROUT('INIT',1,' ')
      ENDIF
      IF(IOPTI.NE.0.OR.IOPTO.NE.0)THEN
         CALL RZEND(CHTOP)
         CLOSE(LUN)
      ENDIF
*
  99  END
