*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDLENS(ISEL1,KLSA)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Zoom by the factors ZOOMFU,ZOOMFV                        *
C.    *       around screen coordinates UUZ0,VVZ0                      *
C.    *       and places the resulting picture at UU0,VV0              *
C.    *                                                                *
C.    *       After GDZOOM has been CALLed, any drawing operation      *
C.    *       will actually draw accordingly to last GDZOOM parameters *
C.    *                                                                *
C.    *       If ZOOMFU or ZOOMFV is 0 (or negative) the zoom          *
C.    *       parameters (GZUA,...,GZVC in common GCDRAW) are reset    *
C.    *                                                                *
C.    *       It is now possible to show in different windows the      *
C.    *       `zoomed-views' obtained just clicking in different       *
C.    *       points of the main picture or of the zoomed-pictures     *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Author: S. Giani 1992                                    *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
C.
C.    ------------------------------------------------------------------
C.
      DIMENSION RVAL(2)
      CHARACTER*4 KLSA
*
      IPX1=1
      IPXSAV=0
      INGR=0
      ISTR=0
   10 GZUA=1
      GZVA=1
      GZUB=0
      GZVB=0
      GZUC=0
      GZVC=0
      MO=0
      CALL IGLOC2(1,NT,UZ1,VZ1,UZ2,VZ2,ISTAT,'R')
      IF (ISTAT.EQ.0)THEN
***              CALL IGCHWK(9,0,0,600,600)
         GZUA=1
         GZVA=1
         GZUB=0
         GZVB=0
         GZUC=0
         GZVC=0
         GO TO 999
      ENDIF
      PLTX=20.
      PLTY=20.
      IF (UZ2-UZ1.EQ.0.) UZ2=UZ1+PLTX*.005
      IF (VZ2-VZ1.EQ.0.) VZ2=VZ1+PLTY*.005
      ZOOMFU=PLTX/ABS(UZ2-UZ1)
      ZOOMFV=PLTY/ABS(VZ2-VZ1)
      ZZFU=ZOOMFU
      ZZFV=ZOOMFV
      UUZ0=(UZ2+UZ1)*.5
      VVZ0=(VZ2+VZ1)*.5
      UU0=(ABS(UZ2-UZ1))*.5
      VV0=20.-(ABS(VZ2-VZ1))*.5
***         CALL IDAWK(1)
      GZUA=GZUA*ZOOMFU
      GZVA=GZVA*ZOOMFV
      GZUB=(GZUB+GZUC)*ZOOMFU
      GZVB=(GZVB+GZVC)*ZOOMFV
      GZUC=UU0-UUZ0*ZOOMFU
      GZVC=VV0-VVZ0*ZOOMFV
      CALL IGQWK(1,'MXDS',RVAL)
      ISU1=UZ1*RVAL(1)*.05
      ISV1=(20.-VZ1)*RVAL(2)*.05
      ISU21=(UZ2-UZ1)*RVAL(1)*.05
      ISV21=(VZ1-VZ2)*RVAL(2)*.05
      IPXSAV=2+IPXSAV
      CALL IGPXMP(IPXSAV,ISU21,ISV21,'OC')
      CALL IGPXMP(IPXSAV,ISU1,ISV1,'S')
      IPX1=2+IPX1
*                             CALL GDSHOW(ISEL1)
      CALL IGPXMP(IPX1,ISU21,ISV21,'O')
***         CALL IGCHWK(9,ISU1,ISV1,ISU21,ISV21)
***         CALL IOPWK(2,1,9)
***         CALL IACWK(2)
      CALL GDSHOW(ISEL1)
      CALL IGPXMP(IPX1,ISU1,ISV1,'CDR')
   20 CONTINUE
      CALL IRQLC(1,MO,JISTA,NT,PX,PY)
      IF(MO.NE.-2.AND.MO.NE.2)THEN
         IF(JISTA.EQ.0)THEN
***           CALL IDAWK(2)
***           CALL ICLWK(2)
***           CALL IACWK(1)
            IF(INGR.EQ.1)THEN
               CALL IGPXMP(IPXSAV,ISU1,ISV1,'DR')
            ELSE
               CALL IGPXMP(IPXSAV,IPXX,IPYY,'DR')
            ENDIF
            GOTO 10
         ELSEIF(JISTA.EQ.1)THEN
            INGR=1
            ISTR=0
            IF(KLSA.EQ.'ON')THEN
               MO=-2
            ELSE
               MO=2
            ENDIF
            GZUA=1
            GZVA=1
            GZUB=0
            GZVB=0
            GZUC=0
            GZVC=0
            ZOOMFU=(.5*PY*PLTX)/ABS(UZ2-UZ1)
            ZOOMFV=(.5*PY*PLTY)/ABS(VZ2-VZ1)
            ZZFU=ZOOMFU
            ZZFV=ZOOMFV
            GZUA=GZUA*ZOOMFU
            GZVA=GZVA*ZOOMFV
            GZUB=(GZUB+GZUC)*ZOOMFU
            GZVB=(GZVB+GZVC)*ZOOMFV
            GZUC=UU0-UUZ0*ZOOMFU
            GZVC=VV0-VVZ0*ZOOMFV
***           CALL ISFAIS(1)
***           CALL ISFACI(0)
***           CALL IGBOX(0.,20.,20.,0.)
            IPX1=IPX1+2
            CALL IGPXMP(IPX1,ISU21,ISV21,'O')
            CALL GDSHOW(ISEL1)
            CALL IGPXMP(IPX1,ISU1,ISV1,'CDR')
            GOTO 20
         ELSEIF(JISTA.EQ.2)THEN
            ISTR=1
            INGR=0
            IF(KLSA.EQ.'ON')THEN
               MO=-2
            ELSE
               MO=2
            ENDIF
            GZUA=1
            GZVA=1
            GZUB=0
            GZVB=0
            GZUC=0
            GZVC=0
            UUZ0=PX
            VVZ0=PY
            GZUA=GZUA*ZOOMFU
            GZVA=GZVA*ZOOMFV
            GZUB=(GZUB+GZUC)*ZOOMFU
            GZVB=(GZVB+GZVC)*ZOOMFV
            GZUC=UU0-UUZ0*ZOOMFU
            GZVC=VV0-VVZ0*ZOOMFV
            IPX1=IPX1+2
            CALL IGPXMP(IPX1,ISU21,ISV21,'O')
            CALL GDSHOW(ISEL1)
            IPXX=PX*RVAL(1)*.05-ISU21*.5
            IPYY=(20.-PY)*RVAL(2)*.05-ISV21*.5
            CALL IGPXMP(IPX1,IPXX,IPYY,'C')
            CALL IGPXMP(IPXSAV,ISU1,ISV1,'DR')
            IPXSAV=IPXSAV+2
            CALL IGPXMP(IPXSAV,ISU21,ISV21,'OC')
            CALL IGPXMP(IPXSAV,IPXX,IPYY,'S')
            CALL IGPXMP(IPX1,IPXX,IPYY,'DR')
            GOTO 20
         ENDIF
      ELSE
         IF(JISTA.EQ.0)THEN
            IF(INGR.EQ.1)THEN
*            CALL IGPXMP(IPXSAV,ISU1,ISV1,'DR')
               MO=0
               GOTO 20
            ELSE
*            CALL IGPXMP(IPXSAV,IPXX,IPYY,'DR')
               MO=0
               ISU1=IPXX
               ISV1=IPYY
               GOTO 20
            ENDIF
*           GOTO 15
         ELSEIF((JISTA.EQ.-1.AND.INGR.EQ.1).OR. (JISTA.GT.0.AND.INGR.EQ
     +   .1.AND.MO.EQ.2))THEN
            GZUA=1
            GZVA=1
            GZUB=0
            GZVB=0
            GZUC=0
            GZVC=0
            ZOOMFU=(.5*PY*PLTX)/ABS(UZ2-UZ1)
            ZOOMFV=(.5*PY*PLTY)/ABS(VZ2-VZ1)
            ZZFU=ZOOMFU
            ZZFV=ZOOMFV
            GZUA=GZUA*ZOOMFU
            GZVA=GZVA*ZOOMFV
            GZUB=(GZUB+GZUC)*ZOOMFU
            GZVB=(GZVB+GZVC)*ZOOMFV
            GZUC=UU0-UUZ0*ZOOMFU
            GZVC=VV0-VVZ0*ZOOMFV
*          CALL ICLRWK(2,1)
***           CALL ISFAIS(1)
***           CALL ISFACI(0)
***           CALL IGBOX(0.,20.,20.,0.)
            IPX1=IPX1+2
            CALL IGPXMP(IPX1,ISU21,ISV21,'O')
            CALL GDSHOW(ISEL1)
            CALL IGPXMP(IPX1,ISU1,ISV1,'CDR')
            GOTO 20
         ELSEIF((JISTA.EQ.-1.AND.ISTR.EQ.1).OR. (JISTA.GT.0.AND.ISTR.EQ
     +   .1.AND.MO.EQ.2))THEN
            GZUA=1
            GZVA=1
            GZUB=0
            GZVB=0
            GZUC=0
            GZVC=0
            UUZ0=PX
            VVZ0=PY
            GZUA=GZUA*ZOOMFU
            GZVA=GZVA*ZOOMFV
            GZUB=(GZUB+GZUC)*ZOOMFU
            GZVB=(GZVB+GZVC)*ZOOMFV
            GZUC=UU0-UUZ0*ZOOMFU
            GZVC=VV0-VVZ0*ZOOMFV
            IPX1=IPX1+2
            CALL IGPXMP(IPX1,ISU21,ISV21,'O')
            CALL GDSHOW(ISEL1)
            CALL IGPXMP(IPX1,IPXX,IPYY,'C')
            CALL IGPXMP(IPXSAV,IPXX,IPYY,'DR')
            IPXX=PX*RVAL(1)*.05-ISU21*.5
            IPYY=(20.-PY)*RVAL(2)*.05-ISV21*.5
            IPXSAV=IPXSAV+2
            CALL IGPXMP(IPXSAV,ISU21,ISV21,'OC')
            CALL IGPXMP(IPXSAV,IPXX,IPYY,'S')
            CALL IGPXMP(IPX1,IPXX,IPYY,'DR')
            GOTO 20
         ENDIF
      ENDIF
C
  999 END
