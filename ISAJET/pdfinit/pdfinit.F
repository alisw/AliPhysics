      SUBROUTINE PDFINIT(PDFPAR, PDFVAL)
#include "isajet/types.inc" 
#if defined(CERNLIB_PDFLIB)
#include "isajet/w50510.inc" 
#include "isajet/w50517.inc"
      CHARACTER*20 PDFPAR(20)

#if defined(CERNLIB_SINGLE)
      REAL PDFVAL(20)
      REAL DX,DSCALE,DXPDF(-6:6)
#endif
#if defined(CERNLIB_DOUBLE)
      DOUBLE PRECISION PDFVAL(20)
      DOUBLE PRECISION DX,DSCALE,DXPDF(-6:6)
#endif

        IF(ISTRUC.EQ.-999) THEN
          WRITE(ITLIS,1200)
1200      FORMAT(//
     $    '1********************************'/
     $    ' *                              *'/
     $    ' * INITIALIZE PDFLIB FOR ISAJET *'/
     $    ' *                              *'/
     $    ' ********************************'/)
          N6=ITLIS
          IFLPRT=2
          CALL PDFSET(PDFPAR,PDFVAL)
          CALL PFTOPDG(0.5D0,1.0D2,DXPDF)
          IFLPRT=0
        ENDIF

        DO 541 I=1,20
          PDFPAR(I)='                    '
          PDFVAL(I)=0
541     CONTINUE
        READ(ITCOM,*) (PDFPAR(I),PDFVAL(I),I=1,20)
        DO 542 I=1,20
          IF(PDFPAR(I).NE.'                    ') THEN
            WRITE(ITLIS,*) PDFPAR(I),PDFVAL(I)
          ENDIF
542     CONTINUE
        ISTRUC=-999
        LOC(54)=NSEL
#endif
        
        END




