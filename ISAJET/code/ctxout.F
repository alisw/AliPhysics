#include "isajet/pilot.h"
      SUBROUTINE CTXOUT(NVC,VC,MXVC)
C-----------------------------------------------------------------------
C  Purpose:
C          Save the context for an ISAJET job:
C          Save in NVC words of VC all common blocks NOT associated only
C          with a single event. Call this and CTXIN to generate mixed
C          events.
C          PARAMETER (MXVC=20000)
C          REAL    VC(MXVC)
C          ...
C          CALL CTXIN(NVC,VC,MXVC)
C
C          Note that the MSSM common blocks are not saved, so different
C          SUSY runs cannot be mixed.
C
C          Ver. 7.02: Equivalenced dummy variables to avoid mixed 
C                     arguments in MOVLEV or multiple EQUIVALENCEd
C                     arguments to CTXIN/CTXOUT.
C
C  Author:
C          F.E. Paige, April 1992     
C-----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/dkytab.inc"
#include "isajet/dylim.inc"
#include "isajet/dypar.inc"
#include "isajet/eepar.inc"
#include "isajet/final.inc"
#include "isajet/force.inc"
#include "isajet/frgpar.inc"
#include "isajet/hcon.inc"
#include "isajet/idrun.inc"
#include "isajet/isloop.inc"
#include "isajet/itapes.inc"
#include "isajet/jetlim.inc"
#include "isajet/keys.inc"
#include "isajet/limevl.inc"
#include "isajet/lstprt.inc"
#include "isajet/mbgen.inc"
#include "isajet/mbpar.inc"
#include "isajet/nodcay.inc"
#include "isajet/primar.inc"
#include "isajet/prtout.inc"
#include "isajet/ptpar.inc"
#include "isajet/q1q2.inc"
#include "isajet/qcdpar.inc"
#include "isajet/qlmass.inc"
#include "isajet/tcpar.inc"
#include "isajet/times.inc"
#include "isajet/totals.inc"
#include "isajet/types.inc"
#include "isajet/wcon.inc"
C
      INTEGER NVC,MXVC,NC,NN,I
      REAL VC(MXVC)
      CHARACTER*8 CLIST(290)
      EQUIVALENCE (CLIST(1),PARTYP(1))
C
C          Dummy real variables for integers
      REAL VLOOK(MXLOOK+6*MXDKY)
      EQUIVALENCE (VLOOK(1),LOOK(1))
      REAL VNKINF(5)
      EQUIVALENCE (VNKINF(1),NKINF)
      REAL VFORCE(9*MXFORC+1)
      EQUIVALENCE (VFORCE(1),NFORCE)
      REAL VIDVER(5)
      EQUIVALENCE (VIDVER(1),IDVER)
      REAL VEVOLV(4)
      EQUIVALENCE (VEVOLV(1),NEVOLV)
      REAL VITDKY(4)
      EQUIVALENCE (VITDKY(1),ITDKY)
      REAL VIKEYS(12)
      EQUIVALENCE (VIKEYS(1),IKEYS)
      REAL VSTPRT
      EQUIVALENCE (VSTPRT,LSTPRT)
      REAL VNJET(9)
      EQUIVALENCE (VNJET(1),NJET)
      REAL VEVPRT(2)
      EQUIVALENCE (VEVPRT(1),NEVPRT)
      REAL VKINPT(5)
      EQUIVALENCE (VKINPT(1),NKINPT)
      REAL VLOC(100)
      EQUIVALENCE (VLOC(1),LOC(1))
C          Dummy real variables for logicals
      REAL VFLW(13)
      EQUIVALENCE (VFLW(1),FLW)
      REAL VNODCY(6)
      EQUIVALENCE (VNODCY(1),NODCAY)
      REAL VGOQ(3*MXGOQ+135)
      EQUIVALENCE (VGOQ(1),GOQ(1,1))
C
      NC=0
C          DKYTAB
      NN=MXLOOK+6*MXDKY
      CALL MOVLEV(VLOOK(1),VC(NC+1),NN)
      NC=NC+NN
C          DYLIM
      CALL MOVLEV(QMIN,VC(NC+1),24)
      NC=NC+24
C          DYPAR
      CALL MOVLEV(VFLW(1),VC(NC+1),13)
      NC=NC+13
C          EEPAR
      CALL MOVLEV(SGMXEE,VC(NC+1),1)
      NC=NC+1
C          FINAL
      CALL MOVLEV(VNKINF(1),VC(NC+1),5)
      NC=NC+5
C          FORCE
      NN=9*MXFORC+1
      CALL MOVLEV(VFORCE(1),VC(NC+1),NN)
      NC=NC+NN
C          FRGPAR
      CALL MOVLEV(PUD,VC(NC+1),41)
      NC=NC+41
C          HCON
      CALL MOVLEV(HMASS,VC(NC+1),69)
      NC=NC+69
C          IDRUN
      CALL MOVLEV(VIDVER(1),VC(NC+1),5)
      NC=NC+5
C          ISLOOP
      CALL MOVLEV(VEVOLV(1),VC(NC+1),4)
      NC=NC+4
C          ITAPES
      CALL MOVLEV(VITDKY(1),VC(NC+1),4)
      NC=NC+4
C          JETLIM
      CALL MOVLEV(PMIN(1),VC(NC+1),72)
      NC=NC+72
C          KEYS
      CALL MOVLEV(VIKEYS(1),VC(NC+1),12)
      NC=NC+12
      CALL CTXC2I(REAC,VC(NC+1),8)
      NC=NC+8
C          LIMEVL
      CALL MOVLEV(ETTHRS,VC(NC+1),3)
      NC=NC+3
C          LSTPRT
      CALL MOVLEV(VSTPRT,VC(NC+1),1)
      NC=NC+1
C          MBGEN
      NN=4*LIMPOM+8
      CALL MOVLEV(POMWT(1),VC(NC+1),NN)
      NC=NC+NN
C          MBPAR
      CALL MOVLEV(PUD0,VC(NC+1),19)
      NC=NC+19
C          NODCAY
      CALL MOVLEV(VNODCY(1),VC(NC+1),6)
      NC=NC+6
C          PRIMAR
      CALL MOVLEV(VNJET(1),VC(NC+1),9)
      NC=NC+9
C          PRTOUT
      CALL MOVLEV(VEVPRT(1),VC(NC+1),2)
      NC=NC+2
C          PTPAR
      CALL MOVLEV(PTFUN1,VC(NC+1),6)
      NC=NC+6
C          Q1Q2
      CALL MOVLEV(VGOQ(1),VC(NC+1),3*MXGOQ+135)
      NC=NC+3*MXGOQ+135
C          QCDPAR
      CALL MOVLEV(ALAM,VC(NC+1),4)
      NC=NC+4
C          QLMASS
      CALL MOVLEV(AMLEP(1),VC(NC+1),55)
      NC=NC+55
C          TCPAR
      CALL MOVLEV(TCMRHO,VC(NC+1),2)
      NC=NC+2
C          TIMES
      CALL MOVLEV(TIME1,VC(NC+1),2)
      NC=NC+2
C          TOTALS
      CALL MOVLEV(VKINPT(1),VC(NC+1),5)
      NC=NC+5
C          TYPES
      CALL MOVLEV(VLOC(1),VC(NC+1),100)
      NC=NC+100
      DO 100 I=1,290
        CALL CTXC2I(CLIST(I),VC(NC+1),8)
        NC=NC+8
100   CONTINUE
C          WCON
#if defined(CERNLIB_SINGLE)
      NN=514
#endif
#if defined(CERNLIB_DOUBLE)
      NN=514+97
#endif
      CALL MOVLEV(SIN2W,VC(NC+1),NN)
      NC=NC+NN
C
      IF(NC.LE.MXVC) THEN
        NVC=NC
        RETURN
      ELSE
        WRITE(ITLIS,9000) NC
9000    FORMAT(//' ERROR IN CTXOUT, NC = ',I5)
        STOP99
      ENDIF
      END
