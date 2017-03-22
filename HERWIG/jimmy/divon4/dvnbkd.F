*CMZ :          23/08/93  13.30.10  by  Jonathan Butterworth
*-- Author :
      SUBROUTINE DVNBKD
C  *******************************************************************
C  *                                                                 *
C  *                                                                 *
C  *     -   D I V O N N E 4  -  PROGRAM FOR MULTIPLE INTEGRATION    *
C  *                                           AND                   *
C  *                               ADAPTIVE IMPORTANCE SAMPLING      *
C  *                                            BY                   *
C  *                                   NESTED PARTITIONING.          *
C  *                                                                 *
C  *                                                                 *
C  *                                                                 *
C  *                                                                 *
C  *     CODED BY:     J. H. FRIEDMAN, SLAC/CERN.                    *
C  *                                  AND                            *
C  *                   M. H. WRIGHT, STANFORD UNIVERSITY             *
C  *     MODIFIED FOR CERN BY T.LINDELOF, JULY 1981                  *
C  *     FORTRAN 77 UPDATE BY T.LINDELOF, MAY 1982                   *
C  *******************************************************************
C
C   INITIALIZATION ROUTINE TO 'SIMULATE' BLOCK DATA
C     BLOCK DATA
      COMMON /D151DT/ IDATE
      REAL*8 IDATE
      COMMON /PRINT/ IPRINT
      COMMON /ISTRGE/ MXRGNS , ISTOR(12000)
      COMMON /RSTRGE/ RSTSZE,RSTOR(18001)
      INTEGER RSTSZE
      COMMON /QUADRE/ IDEG
      COMMON /START/ ISTART
      COMMON /EXFILE/ NFILE
      COMMON /DISPOS/ IDISP
      COMMON /DEPTHS/ ISTDPH , INCDPH
      COMMON /SAMPLE/ NPOINT
      COMMON /CUTOLS/ BNDTOL, FRACT, REGNTL, FNLTOL
      COMMON /BNDLMT/ FLOBD,FUPBD
      COMMON /PRSTOP/ NSTOP
      COMMON /ZEETRM/ ITRMF
C
      DATA INITL/0/
CPH      SAVE INITL
      SAVE
      IF(INITL.NE.0) RETURN
      INITL=1
C
      CALL JMUCPY(8H17/12/80,IDATE,2)
      IPRINT=1
      MXRGNS=3000
      RSTSZE=18001
      IDEG=0
      ISTART=1
      NFILE=1
      IDISP=0
      ISTDPH=3
      INCDPH=5
      NPOINT=50
      BNDTOL=.05
      FRACT=.5
      REGNTL=.9
      FNLTOL=.1
      FLOBD=-9.9E37
      FUPBD= 9.9E37
      NSTOP=5
      ITRMF=0
C
      CALL DVNOPT
C
      RETURN
      END
