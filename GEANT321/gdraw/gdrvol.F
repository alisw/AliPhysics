*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:27  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDRVOL(N,NNAME,NNUMB,NRS,UTHET,UPHI,UPSI,UU0,UV0,SU,SV)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draws the object NNAME(N),NNUMB(N) with its contents.    *
C.    *       The parameters UTHET,UPHI,...,SV are the same of GDRAW.  *
C.    *       NNAME(1),...,NNAME(N) contain the volume names and       *
C.    *       NNUMB(1),...,NNUMB(N) the volume numbers defining the    *
C.    *       path to go from the top volume to the one to be drawn.   *
C.    *                                                                *
C.    *       This routine differs from GDRAW in the following :       *
C.    *       - the object to be drawn is identified by a full path,   *
C.    *         giving so the possibility to draw a particular copy    *
C.    *         or division of a volume, or even a volume that has     *
C.    *         more than one mother in the geometry tree              *
C.    *       - the object can be drawn either with respect to the     *
C.    *         MAster Reference System (NRS=0) or with respect to its *
C.    *         Daugther Reference System; in the first case it is     *
C.    *         drawn where it stands in the real world, while in      *
C.    *         the second one it is drawn like GDRAW would do         *
C.    *       - in this latter case, track and hit points will be      *
C.    *         drawn with respect to the DRS of the volume last       *
C.    *         drawn by this routine, and not with respect to the     *
C.    *         MARS as it done normally; to reset to the normal       *
C.    *         case a call with NRS=0 or N=0 is required              *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Author : P.Zanarini    **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcvolu.inc"
#include "geant321/gcdraw.inc"
      DIMENSION NNAME(15),NNUMB(15),GR0(10),GT0(3)
      CHARACTER*4 NAME
      SAVE GR0,GT0
C
      DATA GR0/1.,0.,0.,0.,1.,0.,0.,0.,1.,0./
      DATA GT0/0.,0.,0./
C.
C.    ------------------------------------------------------------------
C.
      IF (N.EQ.0) THEN
C
C             Restore unitary transformation in /GCDRAW/
C
         DO 10 I=1,3
   10    GTRAN0(I)=GT0(I)
         DO 20 I=1,10
   20    GRMAT0(I)=GR0(I)
         GO TO 999
C
      ENDIF
C
      IF (N.GT.15.OR.N.LT.1) GO TO 999
C
C             Save /GCVOLU/ if necessary
C
      IFCVOL=0
      IF (NLEVEL.GT.0) THEN
         CALL GSCVOL
         IFCVOL=1
      ENDIF
C
C             Force GLVOLU to recompute the path
C
      NLEVEL=0
      CALL GLVOLU(N,NNAME,NNUMB,IER)
      IF (IER.NE.0) GO TO 90
C
      IF (NRS.EQ.0) THEN
C
C             Store unitary transformation in /GCDRAW/
C
         DO 30 I=1,3
   30    GTRAN0(I)=GT0(I)
         DO 40 I=1,10
   40    GRMAT0(I)=GR0(I)
C
      ELSE
C
C             Store N-th level transformation in /GCDRAW/
C
         DO 50 I=1,3
   50    GTRAN0(I)=GTRAN(I,N)
         DO 60 I=1,10
   60    GRMAT0(I)=GRMAT(I,N)
C
C             Replace N-th level transformation with the unitary one
C
         DO 70 I=1,3
   70    GTRAN(I,N)=GT0(I)
         DO 80 I=1,10
   80    GRMAT(I,N)=GR0(I)
C
      ENDIF
C
C             NLEVEL < 0 to inform GDRAW that GSCVOL has already been ca
C
      IF (IFCVOL.EQ.1) NLEVEL=-NLEVEL
C
      IDRNUM=1
      CALL UHTOC(NNAME(N),4,NAME,4)
      CALL GDRAW(NAME,UTHET,UPHI,UPSI,UU0,UV0,SU,SV)
      IDRNUM=0
C
   90 CONTINUE
      IF (IFCVOL.EQ.1) THEN
         CALL GFCVOL
      ELSE
         NLEVEL=0
      ENDIF
  999 RETURN
      END
