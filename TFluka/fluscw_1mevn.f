*$ CREATE FLUSCW.FOR
*COPY FLUSCW
*                                                                      *
*=== fluscw ===========================================================*
*                                                                      *
      DOUBLE PRECISION FUNCTION FLUSCW
     #(IJ,PLA,TXX,TYY,TZZ,WEE,XX,YY,ZZ,NREG,IOLREG,LLO,NSURF)
 
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      SAVE
 
      INCLUDE '(PAPROP)'
      INCLUDE '(USRBIN)'
      INCLUDE '(USRBDX)'
      INCLUDE '(USRTRC)'
      INCLUDE '(SCOHLP)'
*
*----------------------------------------------------------------------*
*                                                                      *
* This functions returns neutron, proton and pion  displacement damage *
* weight factors.                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Input variables:                                                 *
*                                                                      *
*           Ij = (generalized) particle code                           *
*          Pla = particle momentum (if > 0), or kinetic energy (if <0 )*
*    Txx,yy,zz = particle direction cosines                            *
*          Wee = particle weight                                       *
*     Xx,Yy,Zz = position                                              *
*         Nreg = (new) region number                                   *
*       Iolreg = (old) region number                                   *
*          Llo = particle generation                                   *
*        Nsurf = transport flag (ignore!)                              *
*                                                                      *
*     Output variables:                                                *
*                                                                      *
*       Fluscw = factor the scored amount will be multiplied by        *
*       Lsczer = logical flag, if true no amount will be scored        *
*                regardless of Fluscw                                  *
*                                                                      *
*     Useful variables (common SCOHLP):                                *
*                                                                      *
*     Flux like binnings/estimators (Fluscw):                          *
*          ISCRNG = 1 --> Boundary crossing estimator                  *
*          ISCRNG = 2 --> Track  length     binning                    *
*          ISCRNG = 3 --> Track  length     estimator                  *
*          ISCRNG = 4 --> Collision density estimator                  *
*          ISCRNG = 5 --> Yield             estimator                  *
*          JSCRNG = # of the binning/estimator                         *
*                                                                      *
*----------------------------------------------------------------------*
*
      LOGICAL LFIRST


      DATA LFIRST /.TRUE./
*
* the default is unit weight
      FLUSCW = ONEONE
      LSCZER = .FALSE.
*
* skip all scorings other than tracklength
      IF (ISCRNG.NE.2) RETURN
* skip all particles other than protons, pions and neutrons
      IF ((IJ.NE.1).AND.(IJ.NE.13).AND.(IJ.NE.14).AND.(IJ.NE.8)) RETURN
*
* read displacement damage factors
      IF (LFIRST) THEN
         WRITE(LUNOUT,1000)
 1000    FORMAT(1X,'FLUSCW:  direct conversion of neutron fluence to',
     &       ' 1 MeV neutron equivalent requested')
         WRITE(LUNOUT,*) ' neutrons'
         OPEN(19,FILE='disdam_n.dat',STATUS='UNKNOWN')
         CALL SKIP(19,7)
         CALL RDXSC(19,LUNOUT,IDNDAM,-3)
         CLOSE(19)
         WRITE(LUNOUT,*) ' protons '
         OPEN(19,FILE='disdam_p.dat',STATUS='UNKNOWN')
         CALL SKIP(19,7)
         CALL RDXSC(19,LUNOUT,IDPDAM,-2)
         CLOSE(19)
         WRITE(LUNOUT,*) ' pions   '
         OPEN(19,FILE='disdam_pi.dat',STATUS='UNKNOWN')
         CALL SKIP(19,7)
         CALL RDXSC(19,LUNOUT,IDODAM,-2)
         CLOSE(19)
         LFIRST = .FALSE.
      ENDIF
*
* should be always called with ekin ( pla < 0 ), 
* but we leave the check for the moment..
      IF (PLA.LT.0.0D0) THEN
         EKIN = ABS(PLA)
      ELSEIF (PLA.GT.0.0D0) THEN
         EKIN = SQRT(PLA**2+AM(IJ)**2)-AM(IJ)
      ELSE
         RETURN
      ENDIF
*
* calculate the weight
      IF (IJ.EQ.1) THEN
         FLUSCW = XSECT(IDPDAM,EKIN,0)
      ELSEIF ((IJ.EQ.13).OR.(IJ.EQ.14)) THEN
         FLUSCW = XSECT(IDODAM,EKIN,0)
      ELSEIF (IJ.EQ.8) THEN
         FLUSCW = XSECT(IDNDAM,EKIN,0)
      ELSE
         STOP ' FLUSCW: inconsistent IJ '
      ENDIF

      RETURN
      END
*
*===rdxsc==============================================================*
*
      SUBROUTINE RDXSC(LINP,LCHK,IDXSEC,MODE)

************************************************************************
*   LINP           logical input unit                                  *
*   LCHK    > 0    cross section data are plotted to unit LCHK         *
*                   in equidistant log. binning using double-log.      *
*                   interpolation                                      *
*                   (otherwise they are not plotted)                   *
*   IDXSEC         index of stored cross section in common             *
*   |MODE|  = 1    cross section data given as histogram               *
*                    i.e. E_lo(i)  sigma(i)                            *
*                         E_hi(i)  sigma(i)                            *
*                   with increasing energy                             *
*           = 2    cross section data given for bin averages with      *
*                   increasing energy                                  *
*           = 3    cross section data given for bin averages with      *
*                   decreasing energy                                  *
*   if MODE < 0 the energy unit is assumed to be MeV (otherwise GeV)   *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      PARAMETER (MAXSDT = 10,
     &           MAXSBI = 1400)
      COMMON /CSECT/ XSEAV(MAXSDT,MAXSBI),XS(MAXSDT,MAXSBI),
     &               NXSBIN(MAXSDT),NXSDT
      DIMENSION TMPXS(2,MAXSBI)

      LOGICAL LFIRST
      DATA LFIRST /.TRUE./

      IF (LFIRST) THEN
         NXSDT = 0
         DO 1 IDT=1,MAXSDT
            NXSBIN(IDT) = 0
            DO 2 IBIN=1,MAXSBI
               XSEAV(IDT,IBIN) = 0.0D0
               XS(IDT,IBIN)    = 0.0D0
    2       CONTINUE
    1    CONTINUE
         LFIRST = .FALSE.
      ENDIF

      NXSDT = NXSDT+1

      NEBIN = 0
   10 CONTINUE
         NEBIN = NEBIN+1
         IF (NEBIN.GT.MAXSBI) STOP ' nebin > maxsbi !!!'
         IF (IABS(MODE).EQ.1) THEN
            READ(LINP,*,END=9000) ELOW,XS(NXSDT,NEBIN)
            READ(LINP,*) ELHI,XS(NXSDT,NEBIN)
            XSEAV(NXSDT,NEBIN) = SQRT(ELOW*ELHI)
            IF (MODE.LT.0) XSEAV(NXSDT,NEBIN) = 1.D-3*XSEAV(NXSDT,NEBIN)
         ELSEIF (IABS(MODE).EQ.2) THEN
            READ(LINP,*,END=9000) XSEAV(NXSDT,NEBIN),XS(NXSDT,NEBIN)
            IF (MODE.LT.0) XSEAV(NXSDT,NEBIN) = 1.D-3*XSEAV(NXSDT,NEBIN)
         ELSEIF (IABS(MODE).EQ.3) THEN
            READ(LINP,*,END=9000) TMPXS(1,NEBIN),TMPXS(2,NEBIN)
            IF (MODE.LT.0) TMPXS(1,NEBIN) = 1.D-3*TMPXS(1,NEBIN)
         ELSE
            STOP ' RDXSC: unsupported mode ! '
         ENDIF
         GOTO 10
 9000 CONTINUE

      IF (IABS(MODE).EQ.3) THEN
         DO 3 I=1,NEBIN-1
            IDX = NEBIN-I
            XSEAV(NXSDT,IDX) = TMPXS(1,I)
            XS(NXSDT,IDX) = TMPXS(2,I)
    3    CONTINUE
      ENDIF

      NXSBIN(NXSDT) = NEBIN-1
      IDXSEC = NXSDT

      IF (LCHK.GT.0) THEN
         AELO = LOG10(XSEAV(IDXSEC,1))
         AEHI = LOG10(XSEAV(IDXSEC,NXSBIN(IDXSEC)))
         DAE  = (AEHI-AELO)/DBLE(NXSBIN(IDXSEC))
         DO 4 I=1,NXSBIN(IDXSEC)+1
            AE = AELO+DBLE(I-1)*DAE
            E  = 10.0D0**AE
            WRITE(LCHK,'(1X,2E15.5)') E,XSECT(IDXSEC,E,0)
    4    CONTINUE
      ENDIF

      RETURN
      END
*
*===xsect==============================================================*
*
      DOUBLE PRECISION FUNCTION XSECT(IDXSEC,E,MODE)

************************************************************************
*   IDXSEC         index of stored cross section in common             *
*   E              energy                                              *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      PARAMETER (MAXSDT = 10,
     &           MAXSBI = 1400)
      COMMON /CSECT/ XSEAV(MAXSDT,MAXSBI),XS(MAXSDT,MAXSBI),
     &               NXSBIN(MAXSDT),NXSDT

      IF (E.LE.XSEAV(IDXSEC,1)) THEN
         CS = 0.0D0
      ELSEIF (E.GT.XSEAV(IDXSEC,NXSBIN(IDXSEC))) THEN
         N  = NXSBIN(IDXSEC)-1
         CS = XSCINT(IDXSEC,E,N)
      ELSE
         DO 1 J=1,NXSBIN(IDXSEC)-1
            IF ((E.GT.XSEAV(IDXSEC,J)).AND.(E.LE.XSEAV(IDXSEC,J+1))) 
     &                                                          THEN
               CS = XSCINT(IDXSEC,E,J)
               GOTO 2
            ENDIF
   1     CONTINUE
         STOP ' xsection value not found '
   2     CONTINUE
      ENDIF
      XSECT = CS

      RETURN
      END
*
*===xscint=============================================================*
*
      DOUBLE PRECISION FUNCTION XSCINT(IDXSEC,E,J)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      PARAMETER (MAXSDT = 10,
     &           MAXSBI = 1400)
      COMMON /CSECT/ XSEAV(MAXSDT,MAXSBI),XS(MAXSDT,MAXSBI),
     &               NXSBIN(MAXSDT),NXSDT

      FAC = (LOG10(E)                -LOG10(XSEAV(IDXSEC,J)))/
     &      (LOG10(XSEAV(IDXSEC,J+1))-LOG10(XSEAV(IDXSEC,J)))
      XSCINT = LOG10(XS(IDXSEC,J))
     &         +(LOG10(XS(IDXSEC,J+1))-LOG10(XS(IDXSEC,J)))*FAC
      XSCINT = 10**(XSCINT)

      RETURN
      END
*
*===skip===============================================================*
*
      SUBROUTINE SKIP(LUNIT,NSKIP)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER (LOUT=6,LLOOK=9)

      CHARACTER*132 ACARD

      IF (NSKIP.EQ.0) RETURN

      DO 1 K=1,NSKIP
         READ(LUNIT,'(A132)') ACARD
    1 CONTINUE

      RETURN
      END
