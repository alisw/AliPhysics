*CMZ :          17/07/98  15.49.05  by  Federico Carminati
*-- Author :
      SUBROUTINE SHEVUT(JTYPE,JEDIT)
c	==============================

c	Event information output; formatted if JTYPE=1, unformatted if
c	JTYPE=0
c	All particles in output if JEDIT=0, final state only if JEDIT=1

*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,SHGENE.
      COMMON /SHGENE/ IEVT,NPI0,NETA,NPIC,NPRO,NKAC,NKA0,NRHO,NOME,
     +                  NPHI,NPSI,NDRY
*KEEP,SHRUNP.
      COMMON /SHRUNP/ VMAJ,IMIN,NRUN,NEVTOT
*KEND.

      IF (JEDIT.EQ.1) CALL LUEDIT(1)

      IF (JTYPE.EQ.0) THEN
        WRITE(2)NRUN,IEVT,N
 	  DO 100 J=1,N
          WRITE(2)K(J,1),K(J,2),K(J,3),K(J,4),K(J,5),
     +              P(J,1),P(J,2),P(J,3),P(J,4),P(J,5)
100	  CONTINUE
      ENDIF

      IF (JTYPE.EQ.1) THEN
        WRITE(2,1001)NRUN,IEVT,N
 	  DO 200 J=1,N
          WRITE(2,1002)K(J,1),K(J,2),K(J,3),K(J,4),K(J,5),
     +                   P(J,1),P(J,2),P(J,3),P(J,4),P(J,5)
200	  CONTINUE
      ENDIF

      RETURN
1001	FORMAT(I6,I6,I6)
1002	FORMAT(5I6,5F12.6)
      END
