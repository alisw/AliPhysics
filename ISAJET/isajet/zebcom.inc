#ifndef CERNLIB_ISAJET_ZEBCOM_INC
#define CERNLIB_ISAJET_ZEBCOM_INC
*
*
* zebcom.inc
*
C
C  ZEBCOM is the main zebra common block for event data storage
C
      INTEGER NNQ,NREF
      PARAMETER (NNQ=200000)
      PARAMETER (NREF=9)
      COMMON/ZEBCOM/IXCOM,IXMAIN,IXDVR,FENCE,LISAE,LISAB,LREF,
     $  ZSTOR,ENDZS
      SAVE /ZEBCOM/
      INTEGER IXCOM    ! store number
     $       ,IXMAIN   ! event division number
     $       ,IXDVR    ! run division number
      INTEGER FENCE(8),LREF(NREF),ZSTOR(NNQ),ENDZS
      INTEGER LISAE     ! pointer to ISAE bank
      INTEGER LISAB    ! pointer to ISAB bank
      REAL Q(NNQ)
      INTEGER IQ(NNQ),LQ(NNQ)
      EQUIVALENCE (LISAE,LQ(1)),(LQ(9),IQ(1),Q(1))
C

#endif
