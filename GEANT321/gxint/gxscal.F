*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GXSCAL(ID,NAME,CHOPT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      To generate a LEGO plot of the scan geometry              *
C.    *  Generates and plot a table of physics quantities such as      *
C.    *  the total number of radiation lengths or interaction lengths  *
C.    *  in function of the SCAN parameters TETA,PHI.                  *
C.    *    CHOPT='O' table is generated at Exit  of volume NAME.       *
C.    *    CHOPT='I' table is generated at Entry of volume NAME        *
C.    *    CHOPT='X' radiation lengths                                 *
C.    *    CHOPT='L' Interaction lengths                               *
C.    *    CHOPT='P' Plot the table                                    *
C.    *  If VOLUME='XXXX' Mother volume is used.                       *
C.    *                                                                *
C.    *       Author:    R.Brun      **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcscan.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcscal.inc"
      CHARACTER*(*) NAME,CHOPT
      CHARACTER*80 CHTITL
      CHARACTER*11 CHCASE(2)
      CHARACTER*6 CHLOC(2)
      DIMENSION IOPT(5)
      EQUIVALENCE (IOPTX,IOPT(1)),(IOPTL,IOPT(2))
      EQUIVALENCE (IOPTI,IOPT(3)),(IOPTO,IOPT(4))
      DATA CHCASE/'Radiation','Absorption'/
      DATA CHLOC/'before','after'/
C.
C.    ------------------------------------------------------------------
C.
      IF(LSCAN.LE.0)THEN
         ID=0
         GO TO 99
      ENDIF
      NTETA=Q(LSCAN+1)
      NPHI     = IQ(LSCAN+1)
      NTETA    = IQ(LSCAN+2)
      MODTET   = IQ(LSCAN+3)
      NSLIST   = IQ(LSCAN+4)
      NSLMAX   = IQ(LSCAN+5)
      PHIMIN   = Q (LSCAN+11)
      PHIMAX   = Q (LSCAN+12)
      TETMIN   = Q (LSCAN+13)
      TETMAX   = Q (LSCAN+14)
      VSCAN(1) = Q (LSCAN+15)
      VSCAN(2) = Q (LSCAN+16)
      VSCAN(3) = Q (LSCAN+17)
      FACTX0   = Q (LSCAN+18)
      FACTL    = Q (LSCAN+19)
      FACTR    = Q (LSCAN+20)
      DO 10 I=1,NSLIST
         ISLIST(I)=IQ(LSCAN+20+I)
  10  CONTINUE
*
      CHTITL=' '
      CALL UOPTC(CHOPT,'XLIO',IOPT)
      IF(IOPTL.EQ.0)IOPTX=1
      IF(IOPTI.EQ.0)IOPTO=1
      ICASE=1+IOPTL
      ILOC =1+IOPTO
      IF(NAME.EQ.'XXXX')THEN
         IHNAME=IQ(JVOLUM+1)
         WRITE(CHTITL,1000)CHCASE(ICASE),IHNAME
 1000    FORMAT('Number of ',A,' lengths in ',A4)
      ELSE
         CALL UCTOH(NAME,IHNAME,4,4)
         ISL=IUCOMP(ISLIST,NSLIST,IHNAME)
         IF(ISL.LE.0)THEN
            PRINT *,' Unknown SCAN name: ',NAME
            ID=0
            GO TO 99
         ENDIF
         WRITE(CHTITL,2000)CHCASE(ICASE),CHLOC(ILOC),IHNAME
 2000    FORMAT('Number of ',A,' lengths ',A,1X,A)
      ENDIF
      CALL HBOOK2(ID,CHTITL,NPHI,PHIMIN,PHIMAX,NTETA,TETMIN,TETMAX,0.)
*
      DPHI=(PHIMAX-PHIMIN)/NPHI
      DTETA=(TETMAX-TETMIN)/NTETA
      DO 50 IPHI=1,NPHI
         LPHI=LQ(LSCAN-IPHI)
         IF(LPHI.LE.0)GO TO 50
         PHI=PHIMIN+DPHI*(IPHI-0.5)
         DO 20 ITETA=1,NTETA
            TETA=TETMIN+DTETA*(ITETA-0.5)
            IDES=IQ(LPHI+ITETA)
            NW=IBITS(IDES,16,6)
            ISCUR=IBITS(IDES,0,16)
            IW1=IQ(LPHI+ISCUR+2*NW-2)
            IW2=IQ(LPHI+ISCUR+2*NW-1)
            SX0 =IBITS(IW2,17,15)/FACTX0
            SABS=IBITS(IW1,22,10)/FACTL
            IF(IOPTX.NE.0)THEN
               W=SX0
            ELSE
               W=SABS
            ENDIF
            CALL HFILL(ID,PHI,TETA,W)
  20     CONTINUE
  50  CONTINUE
*
  99  END
 
