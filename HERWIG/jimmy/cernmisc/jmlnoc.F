      FUNCTION JMLNOC (CHV)
C
C CERN PROGLIB# M507    LENOCC          .VERSION KERNFOR  4.21  890323
C ORIG. March 85, A.Petrilli, re-write 21/02/89, JZ
C
C-    Find last non-blank character in CHV
 
      CHARACTER    CHV*(*)
 
      N = LEN(CHV)
 
      DO 17  JJ= N,1,-1
      IF (CHV(JJ:JJ).NE.' ') GO TO 99
 17    CONTINUE
      JJ = 0
 
 99    JMLNOC = JJ
      RETURN
      END
