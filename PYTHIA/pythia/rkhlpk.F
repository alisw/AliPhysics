 
*==================================================================
 
      SUBROUTINE RKHLPK(NUM,LGL1,LGL2,LLV,LL1,LL2)
      IMPLICIT INTEGER(A-Z)
      SAVE
      DIMENSION CONFIG(32,6)
      DATA INIT/0/
      IF(INIT.EQ.0) THEN
        INIT=1
        MUM=0
        DO 140 GL1=1,-1,-2
          DO 130 GL2=1,-1,-2
            DO 120 LV=1,-1,-2
              DO 110 L1=1,-1,-2
                DO 100 L2=1,-1,-2
                  MUM=MUM+1
                  CONFIG(MUM,1)=GL1
                  CONFIG(MUM,2)=GL2
                  CONFIG(MUM,3)=LV
                  CONFIG(MUM,4)=L1
                  CONFIG(MUM,5)=L2
  100           CONTINUE
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE
  140   CONTINUE
      ENDIF
      LGL1=CONFIG(NUM,1)
      LGL2=CONFIG(NUM,2)
      LLV =CONFIG(NUM,3)
      LL1 =CONFIG(NUM,4)
      LL2 =CONFIG(NUM,5)
      END
