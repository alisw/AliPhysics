*
* $Id$
*
* $Log$
* Revision 1.3  1996/10/30 11:05:22  cernlib
* Conflicts resolved
*
* Revision 1.1.1.2  1996/10/30 08:32:07  cernlib
* Version 7.04
*
* Revision 1.2  1996/05/31 13:57:48  pdflib
* replace # by +
*
* Revision 1.1.1.1  1996/04/12 15:29:01  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C-----------------------------------------------------------------------
#if defined(CERNLIB_SINGLE)
      FUNCTION ALPHAS2(SCALE)
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION FUNCTION ALPHAS2(SCALE)

#endif
C
C **********************************************************************
C *                                                                    *
C *   Program to calculate alfa strong to second order                 *
C *   as a function of lambda of 5 flavours                            *
C *   and the desired number of flavours (NFL)                         *
C *   for the selected set of structure functions which fixes lambda.  *
C *   NFL is limited to max. of 6 flavours.                            *
C *   The value of alfa is matched at the thresholds q = mq.           *
C *   When invoked with NFL < 0, it chooses NFL                        *
C *   as the number of flavours with mass less then q.                 *
C *   (mc = 1.5 GeV, mb = 4.75 GeV, mt = 100 GeV)                      *
C *                                                                    *
C *   Input:   SCALE   = QCD scale in GeV                              *
C *                                                                    *
C *   Output:  ALPHAS2 = alpha strong to second order,                 *
C *                      if not LO = 1                                 *
C *                                                                    *
C *   The variables NPTYPE,NGROUP,NSET,NFL,LO,TMAS,QCDL4,QCDL5         *
C *   should be provided by the user via a call to the                 *
C *   subroutine PDFSET at the initialization phase.                   *
C *                                                                    *
C *   Author:   H. Plothow-Besch (CERN-PPE)                            *
C *                                                                    *
C **********************************************************************
C
#include "pdf/impdp.inc"
C
#include "pdf/w5051p1.inc"
#include "pdf/w5051p2.inc"
#include "pdf/w5051p7.inc"
#include "pdf/w50510.inc"
#include "pdf/w50511.inc"
#include "pdf/w50512.inc"
#include "pdf/w50513.inc"
#include "pdf/w50515.inc"
#include "pdf/w50516.inc"
#include "pdf/w50517.inc"
#include "pdf/w50519.inc"
#include "pdf/w505120.inc"
#include "pdf/w505121.inc"
C
#include "pdf/expdp.inc"
     +    NF,KNF
      CHARACTER*20 PARM(NCHDIM)
#include "pdf/expdp.inc"
     +       VAL(NCHDIM)
      DATA XMC/1.5D0/,XMB/4.75D0/,XMT/100.D0/
      DATA ZEROD/0.D0/,PONED/0.001D0/,ONED/1.D0/,TWOD/2.D0/
C.
C.    initialization of the parameters for the coupling constant
C.
C.    QCDL4  : QCD scale, which is matched to the scale used in the
C.             parton distribution function chosen above for 4 flavours.
C.
C.    QCDL5  : QCD scale, which is matched to the scale used in the
C.              parton distribution function chosen above.
C.
C. be aware, this lambda is given for 5 flavours (mb = 4.75 GeV/c**2) !!
C.
C
C   set the default values
C
      ALPHAS2 = ZEROD
C
      IF(N6.LE.0) N6 = L6
C.
#include "pdf/w50511c.inc"
C.
C User wants new version (4.0 or bigger) of PDFLIB format
      IF (IFLSET.NE.1) THEN
         IF(FIRST) THEN
         WRITE(N6,*) ' Warning : NO initialisation via PDFSET made !!'
         WRITE(N6,*)
     +   '           ALL Parameters set to default (Nucleon PDFs) !!'
         ENDIF
         NPTYPE = LPTYPE
         NGROUP = LGROUP
         NSET   = LNSET
C
         PARM(1) = 'NPTYPE'
         VAL(1)  = NPTYPE
         PARM(2) = 'NGROUP'
         VAL(2)  = NGROUP
         PARM(3) = 'NSET'
         VAL(3)  = NSET
C
         CALL PDFSET(PARM,VAL)
      ENDIF
C.
C User wants old version (3.0 or less) of PDFLIB format
      IF(.NOT.NEWVER) THEN
        IF (MODE .GE.0 .AND. MODE .LE.MODEMX) THEN
           IF(MODE.EQ.0) THEN
             NPTYPE = LPTYPO
             NGROUP = LGROPO
             NSET   = LNSETO
           ELSE
             NPTYPE = NPTYCR(MODE)
             NGROUP = NGROCR(MODE)
             NSET   = NSETCR(MODE)
C   Check on validity of parameter values
             IF (NPTYPE.LT.0 .OR. NGROUP.LT.0 .OR. NSET.LT.0) THEN
                IF(FIRST) THEN
                WRITE(N6,*) ' PDFLIB  : MODE value INCORRECT, MODE = ',
     +          MODE
               WRITE(N6,*) ' Warning : MODE value set to OLD default !!'
                ENDIF
                NPTYPE = LPTYPO
                NGROUP = LGROPO
                NSET   = LNSETO
             ENDIF
           ENDIF
        ENDIF
      ENDIF
C.
         IF (QCDL5 .LE. PONED) THEN
            IF(FIRST) THEN
            WRITE (N6,*) ' PDFLIB: QCD-Lambda5 is not defined !! '
            ENDIF
            ALPHAS2 = ZEROD
            RETURN
         ENDIF
C.
        PI=4.0D0*ATAN(ONED)
        B6  = (33.D0-2.D0*6.D0)/PI/12.D0
        BP6 = (153.D0 - 19.D0*6.D0) / PI / TWOD / (33.D0 - 2.D0*6.D0)
        B5  = (33.D0-2.D0*5.D0)/PI/12.D0
        BP5 = (153.D0 - 19.D0*5.D0) / PI / TWOD / (33.D0 - 2.D0*5.D0)
        B4  = (33.D0-2.D0*4.D0)/PI/12.D0
        BP4 = (153.D0 - 19.D0*4.D0) / PI / TWOD / (33.D0 - 2.D0*4.D0)
        B3  = (33.D0-2.D0*3.D0)/PI/12.D0
        BP3 = (153.D0 - 19.D0*3.D0) / PI / TWOD / (33.D0 - 2.D0*3.D0)
        XLC = TWOD * LOG( XMC/QCDL5)
        XLB = TWOD * LOG( XMB/QCDL5)
        XLT = TWOD * LOG( XMT/QCDL5 * TMAS/XMT)
        XLLC = LOG( XLC)
        XLLB = LOG( XLB)
        XLLT = LOG( XLT)
        C65  =  ONED/( ONED/(B5 * XLT) - XLLT*BP5/(B5 * XLT)**2 )
     +        - ONED/( ONED/(B6 * XLT) - XLLT*BP6/(B6 * XLT)**2 )
        C45  =  ONED/( ONED/(B5 * XLB) - XLLB*BP5/(B5 * XLB)**2 )
     +        - ONED/( ONED/(B4 * XLB) - XLLB*BP4/(B4 * XLB)**2 )
        C35  =  ONED/( ONED/(B4 * XLC) - XLLC*BP4/(B4 * XLC)**2 )
     +        - ONED/( ONED/(B3 * XLC) - XLLC*BP3/(B3 * XLC)**2 ) + C45
C
      Q   = SCALE
      XLQ = TWOD *  LOG( Q/QCDL5 )
      XLLQ =  LOG( XLQ )
      KNF = NFL
      NF = KNF
      IF  ( NF .LT. ZEROD) THEN
        IF      ( Q .GT. XMT * TMAS/XMT) THEN
          NF = 6.D0
        ELSEIF  ( Q .GT. XMB ) THEN
          NF = 5.D0
        ELSEIF  ( Q .GT. XMC ) THEN
          NF = 4.D0
        ELSE
          NF = 3.D0
        ENDIF
      ENDIF
      IF(NF .GT. 6.D0) NF = 6.D0
      IF      ( NF .EQ. 6.D0 ) THEN
       ALF = ONED/(ONED/(ONED/(B6*XLQ)- BP6/(B6*XLQ)**2*XLLQ) + C65)
       IF (LO.EQ.1) ALF = ONED/B6/XLQ
      ELSEIF  ( NF .EQ. 5.D0 ) THEN
       ALF = ONED/(B5 * XLQ) -  BP5/(B5 * XLQ)**2 * XLLQ
       IF (LO.EQ.1) ALF = ONED/B5/XLQ
      ELSEIF  ( NF .EQ. 4.D0 ) THEN
       ALF = ONED/(ONED/(ONED/(B4*XLQ)- BP4/(B4*XLQ)**2*XLLQ) + C45)
       IF (LO.EQ.1) ALF = ONED/B4/XLQ
      ELSEIF  ( NF .EQ. 3.D0 ) THEN
       ALF = ONED/(ONED/(ONED/(B3*XLQ)- BP3/(B3*XLQ)**2*XLLQ) + C35)
       IF (LO.EQ.1) ALF = ONED/B3/XLQ
      ELSE
       WRITE(N6,*)'Error in Alphas2: unimplemented # of light flavours',
     +             NFL
       STOP
      ENDIF
      ALPHAS2 = ALF
      RETURN
      END
