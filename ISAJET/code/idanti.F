#include "isajet/pilot.h"
      INTEGER FUNCTION IDANTI(ID)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods :
C-     Return value of antiparticle id
C-
C-   Inputs  :
C-   ID = particle id
C-
C-   Created   1-JUN-1988   Serban D. Protopopescu
C-   3-Jan-1993: Expand self-conjugate list for MSSM and simplify
C                structure.  FEP
C    17-Mar-1997: Correctly handle mesons with IDENT>10000
C-
C----------------------------------------------------------------------
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
      INTEGER ID,IFL1,IFL2,IFL3,IDABS
      INTEGER NSELF,I
      PARAMETER (NSELF=14)
      INTEGER IDSELF(NSELF)
      SAVE IDSELF
      DATA IDSELF/9,10,20,29,30,40,50,60,81,82,83,84,90,91/
C----------------------------------------------------------------------
      IDABS=IABS(ID)
      IFL1=MOD(IDABS/1000,10)
C
C          Baryons and diquarks
C
      IF(IFL1.NE.0) THEN
        IDANTI=-ID
        RETURN
      ENDIF
C
C          Mesons
C
      IF(IDABS.GT.100.AND.IFL1.EQ.0) THEN
        IFL2=MOD(IDABS/100,10)
        IFL3=MOD(IDABS/10,10)
        IF(IFL2.EQ.IFL3) THEN
          IDANTI=+ID
        ELSE
          IDANTI=-ID
        ENDIF
        RETURN
      ENDIF
C
C          Other particles
C
      DO 100 I=1,NSELF
        IF(IDABS.EQ.IDSELF(I)) THEN
          IDANTI=+ID
          RETURN
        ENDIF
100   CONTINUE
      IDANTI=-ID
      RETURN
      END
