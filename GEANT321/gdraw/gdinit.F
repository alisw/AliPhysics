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
      SUBROUTINE GDINIT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Initialize the drawing package                           *
C.    *                                                                *
C.    *    ==>Called by : <USER>  UGINIT                               *
C.    *       Author : P.Zanarini ; S.Giani   *********                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcgobj.inc"
#include "geant321/gcmutr.inc"
#include "geant321/gcspee.inc"
      COMMON/INIFIR/NFIRST
      DIMENSION GR0(10),GT0(3)
#if defined(CERNLIB_CG)
      DIMENSION     V(3),T(4,3)
#endif
C
      SAVE GR0,GT0
      DATA GR0/1.,0.,0.,0.,1.,0.,0.,0.,1.,0./
      DATA GT0/0.,0.,0./
C.
C.    ------------------------------------------------------------------
C.
C             Initialize general drawing parameters:
C
C             NGVIEW = Flag for GDFR3D, set to 0 if one wants GDFR3D to
C                      compute new transformation matrix
C             NKVIEW = Total number of view banks
C             IDVIEW = Actual view bank (0 means screen)
C             GTHETA,GPHI,...,GSCV = Initial values for GDRAW calling parameters
C             GZUA,GZVA,GZUB,GZVB,GZUC,GZVC = Zoom parameters
C             DPERS  = Distance from the center of projection
C                      (i.e. viewing position) at Z=-DPERS
C                      and the projection plane at Z=0
C             LINBUF = Flag to determine if GDRAWV has to perform
C                      line buffering logic (LINBUF=1) or not
C             IPIONS = Flag for GDFR3D, set to 1 if PIONS is used
C             IGMR   = Flag for GDFR3D, set to 1 if APOLLO-GMR is used
C             IPKHIT = Flag for GPHITS (normally 0); if >0 then print only
C                      hit nr. IPKHIT
C             IDRNUM = Flag for GDRAW, set to 1 when called by GDRVOL
C             PLTRNX = Screen and plotter X range, PLTRNX x PLTRNY cm.
C             PLTRNY = Screen and plotter Y range, PLTRNX x PLTRNY cm.
***           RNGU1,RNGU2,RNGV1,RNGV2 = Limits of user coordinates (2D window)
C
      IF(IFINIT(14).NE.0)RETURN
      IFINIT(14)=1
      NFIRST=0
      NGVIEW=0
      NKVIEW = 0
      IDVIEW= 0
      GTHETA=45.
      GPHI=135.
      GPSI=0.
      GU0=10.
      GV0=10.
      GSCU=0.015
      GSCV=0.015
      SINPSI=SIN(GPSI*DEGRAD)
      COSPSI=COS(GPSI*DEGRAD)
      GZUA=1
      GZVA=1
      GZUB=0
      GZVB=0
      GZUC=0
      GZVC=0
      DPERS=1000
      LINBUF=0
      IPIONS=0
      IGMR=0
      IPKHIT=0
      IDRNUM=0
      PLTRNX=20.
      PLTRNY=20.
C
C             LINATT = Current line attributes
C             LINATP = Permanent line attributes
C             ITXATT = Current text attributes
C
      LINCOL=1
      LINWID=1
      LINSTY=1
      LINFIL=0
      LEP=1
C
      LINATT=0
      CALL MVBITS(LINCOL,0,8,LINATT,16)
      CALL MVBITS(LINWID,0,3,LINATT,7)
      CALL MVBITS(LINSTY,0,3,LINATT,10)
      CALL MVBITS(LINFIL,0,3,LINATT,13)
C
      LINATP=LINATT
C
      ITXWID=1
      ITXCOL=1
C
      ITXATT=0
      CALL MVBITS(ITXWID,0,3,ITXATT,0)
      CALL MVBITS(ITXCOL,0,4,ITXATT,3)
C
C             Initialize drawing options:
C
C             ITR3D  = Flag for GDFR3D : ITR3D=0 for standard projection,
C                      ITR3D <> 0 for rotation + projection i.e. R-Z projection
C                      (refer to the header comments of GDFR3D
C                      to have a detailed explanation of the steps
C                      from 3D to 2D)
C
C             ITHRZ  = Flag for GDHITS/GDCXYZ/GDXYZ : if not set to 'OFF ' all t
C                      hits and tracks drawing routines will apply an R-Z
C                      projection in case of SIDE or TOP view
C
C             IOBJ   = Type of object being currently drawn :
C
C                         1 : VOLUME
C                         2 : TRACK
C                         3 : HIT
C                         0 : None or any other object (axes,header,text,etc.)
C
C             IPRJ   = Flag for GDFR3D :
C
C                      'PARA' : Ortographic parallel projection (default)
C                      'PERS' : Perspective projection
C
C             ITRKOP = Flag for GDXYZ :
C
C                      'LINE' : Lines joining track points are drawn
C                      'POIN' : Only the track points are drawn
C
      ITR3D=0
      CALL UCTOH('OFF ',ITHRZ,4,4)
      IOBJ=0
      CALL UCTOH('PARA',IPRJ,4,4)
      CALL UCTOH('LINE',ITRKOP,4,4)
      CALL UCTOH('OFF ',IHIDEN,4,4)
C
C             Initialize Color look up table
C
      DO 10 I=1,16
  10  LOOKTB(I)=I
C
C             Define tree parameters:
C
C             NWCUT = Max words allocated to be used by the cut routines
C
      NWCUT=5000
C
C             Initialize to unitary transformation
C
      DO 20 I=1,3
   20 GTRAN0(I)=GT0(I)
      DO 30 I=1,10
   30 GRMAT0(I)=GR0(I)
C
C             Set display and user coordinate systems
C
      CALL ISWN(1,0.,PLTRNX,0.,PLTRNY)
      CALL ISVP(1,0.,1.,0.,1.)
      CALL IGSET('TXFP',-61.)
      CALL IGRNG(PLTRNX,PLTRNY)
      CALL IGSA(0)
*JS
*
*  Setting the Viewing parameters to THETA =0., PHI=0., PSI=0.
*
#if defined(CERNLIB_CG)
      NCVOLS=0
      NSHIFT=0
      KSHIFT=0
      IHOLE=0
      CGXMIN=-10000.
      CGXMAX=-9999.
      CGYMIN=-10000.
      CGYMAX=-9999.
      CGZMIN=-10000.
      CGZMAX=-9999.
      V(1)=0.
      V(2)=0.
      V(3)=0.
      CALL GDCGVW(V,T)
      CALL CGTSET(NTRCG,T,IREP)
#endif
*JS
C
C             Book JDRAW structure for view banks
C
      IF (JDRAW.LE.0) CALL MZBOOK(IXCONS,JDRAW,JDRAW,1,'DRAW',0,0,0,3,0)
C
      END
 
