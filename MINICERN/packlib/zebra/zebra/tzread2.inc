*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:56  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_TZREAD2_INC
#define CERNLIB_ZEBRA_TZREAD2_INC
*
*
* tzread2.inc
*
#if !defined(CERNLIB_CRAZYERR)
      READ (LUNTQ,LINE(JFMTC1:JFMTC2),END=96) (LQ(J),J=LPUTA,LPUTE-1)
#else
      IF (NAME(5).NE.3)  THEN
          READ (LUNTQ,LINE(JFMTC1:JFMTC2),ERR=64,END=96)
     +                (LQ(J),J=LPUTA,LPUTE-1)
        ELSE
          READ (LUNTQ,LINE(JFMTC1:JFMTC2),ERR=64,END=96)
     +                (Q(J),J=LPUTA-8,LPUTE-9)
        ENDIF
   64 CONTINUE
#endif

#endif
