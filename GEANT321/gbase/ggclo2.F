*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_OLD)
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GGCLOS
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Closes off the geometry setting.                            *
C.    *    Initializes the search list for the contents of each        *
C.    *    volume following the order they have been positioned, and   *
C.    *    inserting the content '0' when a call to GSNEXT (-1) has    *
C.    *    been required by the user.                                  *
C.    *    Performs the development of the JVOLUM structure for all    *
C.    *    volumes with variable parameters, by calling GGDVLP.        *
C.    *    Interprets the user calls to GSORD, through GGORD.          *
C.    *    Computes and stores in a bank (next to JVOLUM mother bank)  *
C.    *    the number of levels in the geometrical tree and the        *
C.    *    maximum number of contents per level, by calling GGNLEV.    *
C.    *    Sets status bit for CONCAVE volumes, through GGCAVE.        *
C.    *    Completes the JSET structure with the list of volume names  *
C.    *    which identify uniquely a given physical detector, the      *
C.    *    list of bit numbers to pack the corresponding volume copy   *
C.    *    numbers, and the generic path(s) in the JVOLUM tree,        *
C.    *    through the routine GHCLOS.                                 *
C.    *                                                                *
C.    *    Called by : <USER>                                          *
C.    *    Authors   : R.Brun, F.Bruyant  *********                    *
C.    *                                                                *
C.    *    Modified by S.Egli at 15.9.90: automatic sorting of volumes *
C     *    done by calling GGORDQ for each volume                      *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gclist.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcopti.inc"
      CHARACTER*4 NAME
      LOGICAL BTEST
C.
C.    ------------------------------------------------------------------
C.
*
* *** Stop the run in case of serious anomaly during initialization
*
      IF (IEORUN.NE.0) THEN
         WRITE (CHMAIL, 1001)
         CALL GMAIL (0, 0)
         STOP
      ENDIF
*
      IF (NVOLUM.LE.0) THEN
         WRITE (CHMAIL, 1002) NVOLUM
         CALL GMAIL (0, 0)
         GO TO 999
      ENDIF
*
      NPUSH = NVOLUM -IQ(JVOLUM-2)
      CALL MZPUSH (IXCONS, JVOLUM, NPUSH, NPUSH,'I')
*
* *** Loop over volumes, create default JNear banks as relevant,
*      and release unused bank space
*
      IDO = 0
      DO 80 IVO = 1,NVOLUM
         JVO = LQ(JVOLUM-IVO)
*
* *** Check if Tracking medium has been defined
*
         NMED=Q(JVO+4)
         IF(NMED.LE.0.OR.NMED.GT.IQ(JTMED-2))THEN
            WRITE(CHMAIL,1003)IQ(JVOLUM+IVO)
            CALL GMAIL (0, 0)
         ELSE
            IF(LQ(JTMED-NMED).EQ.0)THEN
               WRITE(CHMAIL,1003)IQ(JVOLUM+IVO)
               CALL GMAIL (0, 0)
            ENDIF
         ENDIF
         IF (BTEST(IQ(JVO),0)) GO TO 80
         IDO = 1
         IQ(JVO) = IBSET(IQ(JVO),0)
         NINL  = IQ(JVO-2)
         NIN   = Q(JVO+3)
         NUSED = IABS(NIN)
         IF (NIN.GT.0) THEN
*           reserve enough additional space for sorted volumes
            IF(NIN.LE.1.OR.NIN.GT.500.OR.IOPTIM.LT.0)THEN
              NUSED=NUSED+1
            ELSE
              NUSED=NUSED+2
            ENDIF
         ENDIF
*
         NPUSH = NUSED -NINL
         DO 90 IN=NINL,NUSED+1,-1
            JIN = LQ(JVO-IN)
            IF(JIN.GT.0) THEN
               CALL MZDROP(IXCONS,JIN,'L')
            ENDIF
  90     CONTINUE
         CALL MZPUSH (IXCONS, JVO, NPUSH, 0, 'I')
         IF (NIN.LE.0) GO TO 80
*
         IF(BTEST(IQ(JVO),3)) THEN
            IZERO=1
         ELSE
            IZERO=0
         ENDIF
         NEL = NIN +IZERO
         JN = LQ(JVO-NIN-1)
         IF(JN.EQ.0) THEN
            CALL MZBOOK (IXCONS,JN,JVO,-NIN-1,'VONE',0,0,NEL+1,2,0)
         ENDIF
         IQ(JN-5) = IVO
         IQ(JN+1) = NEL
         JN = JN +1
         DO 29 I = 1,NIN
            IQ(JN+IZERO+I) = I
   29    CONTINUE
         IF (IZERO.NE.0) IQ(JN+1) = 0
*
   80 CONTINUE
*
      IF (IDO.NE.0) THEN
*
* ***    Perform development of JVOLUM structure where necessary
*
         CALL GGDVLP
*
* ***    Fill GSORD ordering banks if required
*
* Modified by S.Egli to allow GGORDQ to find the optimum sorting for
* all volumes
*
         IF(IOPTIM.GE.1)THEN
            WRITE(6,'(A)')' GGCLOS: Start automatic volume ordering:'
         ENDIF
         DO 91 IVO = 1,NVOLUM
            JVO = LQ(JVOLUM-IVO)
            NIN = Q(JVO+3)
            ISEARC=Q(JVO+1)
            IF(ISEARC.GT.0) GO TO 91
*           check if sorting not possible or not wanted
            IF(NIN.LE.1.OR.NIN.GT.500.OR.IOPTIM.LT.0)THEN
               Q(JVO+1)=0.
               IF(NIN.GT.500.AND.IOPTIM.GE.1)THEN
                 CALL UHTOC(IQ(JVOLUM+IVO),4,NAME,4)
                 WRITE (CHMAIL,1004) NAME,NIN
                 CALL  GMAIL (0, 0)
               ENDIF
            ELSEIF(IOPTIM.EQ.0)THEN
               IF(ISEARC.LT.0)CALL GGORD (IVO)
            ELSEIF(IOPTIM.EQ.1)THEN
               IF(ISEARC.EQ.0) THEN
                  CALL GGORDQ(IVO)
               ELSE
                  CALL GGORD (IVO)
               END IF
            ELSE
               CALL GGORDQ(IVO)
            ENDIF
   91    CONTINUE
*
* ***    Set status bit for concave volumes
*
         CALL GGCAVE
*
* ***    Compute maximum number of levels and of contents per level
*
         CALL GGNLEV
*
      ENDIF
*
* *** Scan the volume structure to retrieve the path through
*      the physical tree for all sensitive detectors
*
       CALL GHCLOS
*
* *** Books STAT banks if data card STAT is submitted
*
      IF (NSTAT.GT.0)  CALL GBSTAT
*
      CALL MZGARB (IXCONS, 0)
*
 1001 FORMAT (' Severe diagnostic in initialization phase. STOP')
 1002 FORMAT (' GGCLOS : NVOLUM =',I5,' *****')
 1003 FORMAT (' Illegal tracking medium number in volume : ',A4)
 1004 FORMAT (' GGORDQ : Volume ',A4,' has more than 500 (',
     +        I3,') daughters ; volume sorting not possible !')
*                                                             END GGCLOS
  999 END
 
#endif
