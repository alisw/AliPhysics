*CMZ :          17/07/98  15.49.56  by  Federico Carminati
*-- Author :
      SUBROUTINE SH2GLW
C	=================

C	Output in GLHID format, weights are also output

*KEEP,SHRUNP.
      COMMON /SHRUNP/ VMAJ,IMIN,NRUN,NEVTOT
*KEEP,SHGENE.
      COMMON /SHGENE/ IEVT,NPI0,NETA,NPIC,NPRO,NKAC,NKA0,NRHO,NOME,
     +                  NPHI,NPSI,NDRY
*KEEP,SHWATE.
      COMMON /SHWATE/ WEI(200000)

*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEND.

      CHARACTER   CODEP*16

      FTHET = 1.E6
      FPHI  = 1.E6
      FP    = 1.E5
      FE    = 1.E5
      FW    = 1.E8

      CALL LUEDIT(1)

      KDAT  = 0
      KTIM  = 0
      KRUN  = NRUN
      KEVT  = IEVT
      KPART = N

      KPROJ = 0
      KTARG = 0
      KZPRO = 0
      KZTAR = 0
      KDECA = 0

      KMPAC = 0
      KPREF = 0
      KECM  = 0
      KDB   = 0
      KDM   = 0

      WRITE(7) KDAT,KTIM,KRUN,KEVT,KPART
      WRITE(7) KPROJ,KTARG,KZPRO,KZTAR,KDECA
      WRITE(7) KMPAC,KPREF,KECM,KDB,KDM


      DO 10 JPART=1,N
        IF (K(JPART,2).EQ.22) THEN
          IPART = 1
        ELSE IF (K(JPART,2).EQ.-11) THEN
          IPART = 2
        ELSE IF (K(JPART,2).EQ.11) THEN
          IPART = 3
        ELSE IF (ABS(K(JPART,2)).EQ.12) THEN
            IPART = 4
        ELSE IF (ABS(K(JPART,2)).EQ.14) THEN
          IPART = 4
        ELSE IF (ABS(K(JPART,2)).EQ.16) THEN
          IPART = 4
        ELSE IF (K(JPART,2).EQ.-13) THEN
          IPART = 5
        ELSE IF (K(JPART,2).EQ.13) THEN
          IPART = 6
        ELSE IF (K(JPART,2).EQ.111) THEN
          IPART = 7
        ELSE IF (K(JPART,2).EQ.211) THEN
          IPART = 8
        ELSE IF (K(JPART,2).EQ.-211) THEN
          IPART = 9
        ELSE IF (K(JPART,2).EQ.130) THEN
          IPART = 10
        ELSE IF (K(JPART,2).EQ.321) THEN
          IPART = 11
        ELSE IF (K(JPART,2).EQ.-321) THEN
          IPART = 12
        ELSE IF (K(JPART,2).EQ.2112) THEN
          IPART = 13
        ELSE IF (K(JPART,2).EQ.2212) THEN
          IPART = 14
        ELSE IF (K(JPART,2).EQ.-2212) THEN
          IPART = 15
        ELSE IF (K(JPART,2).EQ.310) THEN
          IPART = 16
        ELSE
          CALL LUNAME(KFA,CODEP)
          WRITE(MSTU(11),*)'ERROR:'
          WRITE(MSTU(11),*)CODEP,'NOT generated with JWEI=1'
          WRITE(MSTU(11),*)'EXECUTION STOPPED!'
        ENDIF

        THETA = PLU(JPART,14)
        PHI   = PLU(JPART,16)
        PP    = PLU(JPART,8)
        EE    = P(JPART,4)
        WW    = WEI(JPART)

        KTHET = NINT(THETA*FTHET)
        KPHI  = NINT(PHI*FPHI)
        KP    = NINT(PP*FP)
        KE    = NINT(EE*FE)
        KW    = NINT(WW*FW)

        WRITE(7) IPART,KTHET,KPHI,KP,KE,KW

10	CONTINUE
      RETURN
      END
