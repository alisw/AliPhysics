C...  ALICE interface to PDFLIB with possibility to select nuclear structure 
C...  functions. 
C...  
C...  The MSTP array in the PYPARS common block is used to enable and 
C...  select the nuclear structure functions. 
C...  MSTP(52)  : (D=1) choice of proton and nuclear structure-function library
C...          =1: internal PYTHIA acording to MSTP(51) 
C...          =2: PDFLIB proton  s.f., with MSTP(51)  = 1000xNGROUP+NSET
C...              MSTP( 51)  = 1000xNPGROUP+NPSET
C...              MSTP(151)  = 1000xNAGROUP+NASET
C...  MSTP(192) : Mass number of nucleus side 1
C...  MSTP(193) : Mass number of nucleus side 2
C...
C...
C...  MINT(124) : side (1 or 2)


      SUBROUTINE PDFSET_ALICE(PARM, VALUE)
C...
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Interface to PDFLIB.
      COMMON/LW50512/QCDL4,QCDL5
      SAVE /LW50512/
      DOUBLE PRECISION QCDL4,QCDL5
      COMMON/LW50513/XMIN,XMAX,Q2MIN,Q2MAX
      SAVE /LW50513/
      DOUBLE PRECISION XMIN,XMAX,Q2MIN,Q2MAX
C...
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)  
      DOUBLE PRECISION VALUE(20)
      CHARACTER*20 PARM(20)
      write(6,*) MSTP(52)
      write(6,*) PARM
      write(6,*) VALUE

      IF (MSTP(192) .GT. 0 .AND. MSTP(193) .GT. 0) THEN
         PARM(5)='NATYPE'
         VALUE(5)=4
         PARM(6)='NAGROUP'
         VALUE(6)=MSTP(191)/1000
         PARM(7)='NASET'
         VALUE(7)=MOD(MSTP(191),1000)
         CALL PDFSET(PARM,VALUE,
     >        MSTU(11),MSTP(51),MSTP(53),MSTP(55),
     >        QCDL4,QCDL5,
     >        XMIN,XMAX,Q2MIN,Q2MAX)

         IF (MSTP(194) .EQ. 0) THEN 
            CALL SETLHAPARM("EKS98")
         ELSE IF (MSTP(194) .EQ.  9) THEN
            CALL SETLHAPARM("EPS09LO")
         ELSE IF (MSTP(194) .EQ. 19) THEN
            CALL SETLHAPARM("EPS09NLO")
         ELSE IF (MSTP(194) .EQ.  8) THEN
            CALL SETLHAPARM("EPS08")
         ELSE
            CALL SETLHAPARM("EPS09LO")
         ENDIF
       ELSE 
         write(6,*) "-> pdfset"
         CALL PDFSET(PARM,VALUE,
     >        MSTU(11),MSTP(51),MSTP(53),MSTP(55),
     >        QCDL4,QCDL5,
     >        XMIN,XMAX,Q2MIN,Q2MAX)
      ENDIF
      write(6,*) "done"
      END



      SUBROUTINE STRUCTM_ALICE
     +     (XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
C...
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
C      write(6,*) "structm_alice->"
      IF (MSTP(192) .GT. 0 .AND. MSTP(193) .GT. 0) THEN
         A=MSTP(191+MINT(124))
C         write(6,*) mint(124), "-> structa ", A
          CALL STRUCTA(XX,QQ,A,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
      ELSE
C         write(6,*) mint(124), "-> structm "
         CALL STRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)
      ENDIF
      END


