*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:20  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCGHI(IMOD,NOBJ,NWWS,LSTEP,ISG,NTVOL,ISHAPE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       This routine allows :                                    *
C.    *                                                                *
C.    *          - Slicing by a plane if some cut has been set,        *
C.    *          - Clipping by any kind of volume if one or more       *
C.    *            clipping volumes have been defined                  *
C.    *          - to Insert the CG Objects into the Hidden Structure  *
C.    *          - to convert the CG Objects in Wire Objects           *
C.    *                                                                *
C.    *       Input Parameters:                                        *
C.    *                                                                *
C.    *          IMOD:  Number for indicating the task to be perfor-   *
C.    *                 med.  See GDCGOB subroutine                    *
C.    *          NOBJ:  Counter for CG objects                         *
C.    *          NWWS:  Size of Wire Structure                         *
C.    *          LSTEP:  N. of bodies forming each CG volume           *
C.    *          ISG  :  i-th body forming a Pcon or Pgon              *
C.    *          NTVOL: N. of bodies forming each Pcon or Pgon         *
C.    *          ISHAPE:Shape of CG object                             *
C.    *                                                                *
C.    *                                                                *
C.    *    ==>Called by : GDCGOB                                       *
C.    *                                                                *
C.    *       Authors :  J.Salt, S.Giani                               *
C.    ******************************************************************
C.
*
#include "geant321/gcbank.inc"
#include "geant321/gcgobj.inc"
#include "geant321/gcunit.inc"
#include "geant321/gchiln.inc"
#include "geant321/gcspee.inc"
#include "geant321/gcmutr.inc"
      SAVE NFIRST, NSHA, IOLDNS
      DATA NFIRST/1/
      DATA NSHA/0/
      DATA IOLDNS/0/
*
*   Slicing or clipping Actions
*
      LLEP=ABS(LEP)
      IVFUN=1
      LENGHT=4000
      ICPOIN=JCGOBJ+1
      IZPOIN=JCGOBJ+20000
      IF(ISHAPE.EQ.11.OR.ISHAPE.EQ.12)ICGP=ICGP+1
*
      IF(IMOD.EQ.3)THEN
         IST=2
      ELSE IF(IMOD.EQ.4)THEN
         IST=3
         CALL CGSLIC(Q(ICPOIN),ABCD,LENGHT,Q(IZPOIN))
**SG
      ELSE IF(IMOD.EQ.5)THEN
         IST=3
*  Clip object
         CALL GDCGCL(ISHAPE)
**SG
      ELSE
         WRITE(CHMAIL,10000)IMOD
         CALL GMAIL(0,0)
      ENDIF
*
**      IF(IVFUN.EQ.0)GOTO 999
*
***SG
*      Inserting CG object in the Hidden Structure
*
      IHPOIN=JCG+1
      ICPOIN=JCGOBJ+1
      ISPOIN=JCGOBJ+8000
      IZPOIN=JCGOBJ+20000
      IF(ISHAPE.EQ.29.OR.ISHAPE.EQ.99)ICPOIN=ISPOIN
      IF(IMOD.EQ.4.OR.IMOD.EQ.5)ICPOIN=IZPOIN
*
*   Evaluate size of Hide Structure after all volumes
*   If it's too big, go back to GDRAW
*   If it's ok, build the Wire Structure
*
      CALL CGHINS(Q(ICPOIN),Q(IHPOIN),ISHAPE)
      IF(KCGST.EQ.-8)THEN
         IF(NFIRST.EQ.1)THEN
           NFIRST=0
           WRITE(CHMAIL,10002)
           CALL GMAIL(0,0)
         ENDIF
         NNN=NOBJ+1
         IF(ISHAPE.EQ.11.OR.ISHAPE.EQ.12)THEN
           IF(ISG.EQ.NTVOL)NNN=NNN+LSTEP-1
         ENDIF
         IF(NNN.EQ.NCLAS2)THEN
           KCGST=-9
           NFIRST=1
         ENDIF
         GOTO 999
      ENDIF
*
*      Building Up the Wiring structure.
*
      ILCG=ILCG+1
      IST=ILCG
      IF(IST.GE.IQ(JCOUNT-1))THEN
        WRITE(CHMAIL,10001)
        CALL GMAIL(0,0)
        GOTO 999
      ENDIF
      IF(LLEP.NE.1)THEN
       IQ(IMCOUN+IST)=IQ(IMCOUN+IST-1)+NSHA
      ENDIF
      IQ(JCOUNT+IST)=IQ(JCOUNT+IST-1)+NTCUR
      IF(IFCG.EQ.4.AND.ILCG.EQ.4)THEN
        IQ(JCOUNT+4)=IQ(JCOUNT+4)-8000+1
        IF(LLEP.NE.1)THEN
         IQ(IMCOUN+4)=IQ(IMCOUN+4)-8000+1
        ENDIF
      ENDIF
      IF(IFCG.EQ.3.AND.ILCG.EQ.3)THEN
        IQ(JCOUNT+3)=IQ(JCOUNT+3)-4000+1
        IF(LLEP.NE.1)THEN
         IQ(IMCOUN+3)=IQ(IMCOUN+3)-4000+1
        ENDIF
      ENDIF
      IF(LLEP.NE.1)THEN
       MMPOIN=IMPOIN+IQ(IMCOUN+IST)
      ENDIF
      IWPOIN=JCG+IQ(JCOUNT+IST)
      ICPOIN=JCGOBJ+1
      IF(ISHAPE.EQ.29.OR.ISHAPE.EQ.99) ICPOIN=JCGOBJ+8000
      IF(IMOD.EQ.4.OR.IMOD.EQ.5) ICPOIN=JCGOBJ+20000
*
*    Evaluate size of Wire structure after all volumes
*    If it's too big, go back to GDRAW
*    If it's ok, go on.
*
      IF(LLEP.NE.1)THEN
      CALL CGWIRE(Q(ICPOIN),NTRCG,-1,NWWS,Q(IWPOIN),ISHAPE,IQ(MMPOIN))
      ELSE
      CALL CGWIRE(Q(ICPOIN),NTRCG,-1,NWWS,Q(IWPOIN),ISHAPE,0)
      ENDIF
      IF(KCGST.EQ.-8)THEN
         IF(NFIRST.EQ.1)THEN
           NFIRST=0
           WRITE(CHMAIL,10002)
           CALL GMAIL(0,0)
         ENDIF
        NNN=NOBJ+1
        IF(ISHAPE.EQ.11.OR.ISHAPE.EQ.12)THEN
           IF(ISG.EQ.NTVOL)NNN=NNN+LSTEP-1
        ENDIF
        IF(NNN.EQ.NCLAS2)THEN
          KCGST=-10
          NFIRST=1
        ENDIF
        NTCUR=0
        GOTO 999
      ENDIF
      ITSTCU=Q(IWPOIN)
      IF(LLEP.NE.1)NSHACU=IQ(MMPOIN)
      IF(ITSTCU.NE.0)THEN
       IOLDCU=ITSTCU
       IF(LLEP.NE.1)IOLDNS=NSHACU
      ENDIF
      NTCUR=Q(IWPOIN)
      IF(LLEP.NE.1)NSHA=IQ(MMPOIN)
      NFILT=NFILT+NTCUR
      IF(LLEP.NE.1)NTNEX=NTNEX+NSHA
      IF(NTCUR.EQ.0)THEN
        NTCUR=IOLDCU
        IF(LLEP.NE.1)NSHA=IOLDNS
        IF(IFCG.EQ.4.AND.ILCG.EQ.4)THEN
          IQ(JCOUNT+4)=IQ(JCOUNT+4)-8000+1
          IF(LLEP.NE.1)THEN
           IQ(IMCOUN+4)=IQ(IMCOUN+4)-8000+1
          ENDIF
        ENDIF
        IF(IFCG.EQ.3.AND.ILCG.EQ.3)THEN
          IQ(JCOUNT+3)=IQ(JCOUNT+3)-4000+1
          IF(LLEP.NE.1)THEN
           IQ(IMCOUN+3)=IQ(IMCOUN+3)-4000+1
          ENDIF
        ENDIF
        ILCG=ILCG-1
        ICGP=ICGP-1
      ENDIF
***SG
10000    FORMAT('  IMOD  = ',I5,' is not defined  ')
10001    FORMAT(' *** GDCGHI *** : Please, increase size of',
     +          ' the Zebra store; the drawing',
     +          ' will not be completed.')
10002    FORMAT(' *** GDCGHI *** : The memory size is not enough;',
     +          ' the program is going on evaluating the number of',
     +          ' words missing.')
*
  999 END
