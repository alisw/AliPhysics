*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE INDEX2.FOR
*COPY INDEX2
*
*=== index2 ===========================================================*
*
      SUBROUTINE INDEX2(KA,KB,IND)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
      KP=KA*KB
      KS=KA+KB
      IF(KP.EQ.1)IND=1
      IF(KP.EQ.2)IND=2
      IF(KP.EQ.3)IND=3
      IF(KP.EQ.4.AND.KS.EQ.5)IND=4
      IF(KP.EQ.5)IND=5
      IF(KP.EQ.6.AND.KS.EQ.7)IND=6
      IF(KP.EQ.4.AND.KS.EQ.4)IND=7
      IF(KP.EQ.6.AND.KS.EQ.5)IND=8
      IF(KP.EQ.8)IND=9
      IF(KP.EQ.10)IND=10
      IF(KP.EQ.12.AND.KS.EQ.8)IND=11
      IF(KP.EQ.9)IND=12
      IF(KP.EQ.12.AND.KS.EQ.7)IND=13
      IF(KP.EQ.15)IND=14
      IF(KP.EQ.18)IND=15
      IF(KP.EQ.16)IND=16
      IF(KP.EQ.20)IND=17
      IF(KP.EQ.24)IND=18
      IF(KP.EQ.25)IND=19
      IF(KP.EQ.30)IND=20
      IF(KP.EQ.36)IND=21
      RETURN
      END
