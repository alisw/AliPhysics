*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:28:58  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:25  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE PDFSTA
C
C **********************************************************************
C *                                                                    *
C *   Routine to print status of parameters in the COMMON/W50511/,     *
C *                                                COMMON/W50512/,     *
C *                                                COMMON/W50513/      *
C *                                                COMMON/W50514/      *
C *                                                COMMON/W505110/     *
C *                                            and COMMON/W505122/     *
C *                                                                    *
C *   Author:  H. Plothow-Besch                                        *
C *            CERN-PPE, CH - 1211 Geneva 23, Switzerland              *
C *                                                                    *
C **********************************************************************
C
#include "pdf/impdp.inc"
C
#include "pdf/w5051p2.inc"
#include "pdf/w5051p7.inc"
#include "pdf/w50510.inc"
#include "pdf/w50511.inc"
#include "pdf/w50512.inc"
#include "pdf/w50513.inc"
#include "pdf/w50514.inc"
#include "pdf/w50515.inc"
#include "pdf/w50517.inc"
#include "pdf/w505110.inc"
#include "pdf/w505122.inc"
      DATA ZEROD/0.D0/
C
      IFLSTA = 1
C.
C print values
      IF(N6.LE.0) N6 = L6
      WRITE(N6,2000)
      WRITE(N6,2001) NPTYPE,NGROUP,NSET,SFNAME(NPTYPE,NGROUP,NSET),
     +               MODECR(NPTYPE,NGROUP,NSET)
      WRITE(N6,2002) NFL,LO,TMAS
      WRITE(N6,2003) QCDL4,QCDL5
      WRITE(N6,2004) XMIN,XMAX,Q2MIN,Q2MAX
 2000 FORMAT(/,'  PDFLIB : Summary from PDFSTA')
 2001 FORMAT(1H ,' Nptype = ',I1,'  Ngroup = ',I1,'  Nset = ',I2,
     +         '  Name  = "',A8,'"  CrMode = ',I3)
 2002 FORMAT(1H ,' Nfl    = ',I2,',  LO = ',I1,
     +           ',  Tmas = ',F7.2,' GeV/c**2')
 2003 FORMAT(1H ,' QCDL4  = ',F7.4,' GeV,  QCDL5 = ',F7.4,' GeV')
 2004 FORMAT(1H ,' Xmin   = ',1X,E8.2,',  Xmax = ',E11.5,
     +          ',  Q2min = ',F7.3,' (GeV/c)**2,  Q2max = ',E8.2,
     +          ' (GeV/c)**2')
C.
      WRITE(N6,3000)
      IF(WXMIN .GT. ZEROD) THEN
         PERC = 100.D0*WXMIN/WTXMIN
         WRITE(N6,3001) WXMIN,WTXMIN,PERC,XMIN
      ENDIF
      IF(WXMAX .GT. ZEROD) THEN
         PERC = 100.D0*WXMAX/WTXMAX
         WRITE(N6,3002) WXMAX,WTXMAX,PERC,XMAX
      ENDIF
      IF(WQ2MIN .GT. ZEROD) THEN
         PERC = 100.D0*WQ2MIN/WTQ2MIN
         WRITE(N6,3003) WQ2MIN,WTQ2MIN,PERC,Q2MIN
      ENDIF
      IF(WQ2MAX .GT. ZEROD) THEN
         PERC = 100.D0*WQ2MAX/WTQ2MAX
         WRITE(N6,3004) WQ2MAX,WTQ2MAX,PERC,Q2MAX
      ENDIF
      IF((WXMIN+WXMAX+WQ2MIN+WQ2MAX) .EQ. ZEROD) WRITE(N6,3005)
 3000 FORMAT(' ')
 3001 FORMAT(1H ,' PDFSTA:    X values ',E9.3,' times out of ',E9.3,
     +         ' ( = ',E9.3,' %) smaller allowed minimum of ',E9.3)
 3002 FORMAT(1H ,' PDFSTA:    X values ',E9.3,' times out of ',E9.3,
     +         ' ( = ',E9.3,' %) bigger  allowed maximum of ',E11.5)
 3003 FORMAT(1H ,' PDFSTA: Q**2 values ',E9.3,' times out of ',E9.3,
     +         ' ( = ',E9.3,' %) smaller allowed minimum of ',E9.3,
     +         ' (GeV/c)**2')
 3004 FORMAT(1H ,' PDFSTA: Q**2 values ',E9.3,' times out of ',E9.3,
     +         ' ( = ',E9.3,' %) bigger  allowed maximum of ',E9.3,
     +         ' (GeV/c)**2')
 3005 FORMAT(1H ,' PDFSTA:    NO errors occured')
C.
      RETURN
      END
