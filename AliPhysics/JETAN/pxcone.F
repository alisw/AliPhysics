      SUBROUTINE PXCONE(MODE,NTRAK,ITKDM,PTRAK,CONER,EPSLON,OVLIM,
     +                   MXJET,NJET,PJET,IPASS,IJMUL,IERR)
*.*********************************************************
*. ------
*. PXCONE
*. ------
*.
*.********** Pre Release Version 26.2.93
*.
*. Driver for the Cone  Jet finding algorithm of L.A. del Pozo.
*. Based on algorithm from D.E. Soper.
*. Finds jets inside cone of half angle CONER with energy > EPSLON.
*. Jets which receive more than a fraction OVLIM of their energy from
*. overlaps with other jets are excluded.
*. Output jets are ordered in energy.
*. If MODE.EQ.2 momenta are stored as (eta,phi,<empty>,pt)
*. Usage     :
*.
*.      INTEGER  ITKDM,MXTRK
*.      PARAMETER  (ITKDM=4.or.more,MXTRK=1.or.more)
*.      INTEGER  MXJET, MXTRAK, MXPROT
*.      PARAMETER  (MXJET=10,MXTRAK=500,MXPROT=500)
*.      INTEGER  IPASS (MXTRAK),IJMUL (MXJET)
*.      INTEGER  NTRAK,NJET,IERR,MODE
*.      DOUBLE PRECISION  PTRAK (ITKDM,MXTRK),PJET (5,MXJET)
*.      DOUBLE PRECISION  CONER, EPSLON, OVLIM
*.      NTRAK = 1.to.MXTRAK
*.      CONER   = ...
*.      EPSLON  = ...
*.      OVLIM   = ...
*.      CALL PXCONE (MODE,NTRAK,ITKDM,PTRAK,CONER,EPSLON,OVLIM,MXJET,
*.     +             NJET,PJET,IPASS,IJMUL,IERR)
*.
*. INPUT     :  MODE      1=>e+e-, 2=>hadron-hadron
*. INPUT     :  NTRAK     Number of particles
*. INPUT     :  ITKDM     First dimension of PTRAK array
*. INPUT     :  PTRAK     Array of particle 4-momenta (Px,Py,Pz,E)
*. INPUT     :  CONER     Cone size (half angle) in radians
*. INPUT     :  EPSLON    Minimum Jet energy (GeV)
*. INPUT     :  OVLIM     Maximum fraction of overlap energy in a jet
*. INPUT     :  MXJET     Maximum possible number of jets
*. OUTPUT    :  NJET      Number of jets found
*. OUTPUT    :  PJET      5-vectors of jets
*. OUTPUT    :  IPASS(k)  Particle k belongs to jet number IPASS(k)
*.                        IPASS = -1 if not assosciated to a jet
*. OUTPUT    :  IJMUL(i)  Jet i contains IJMUL(i) particles
*. OUTPUT    :  IERR      = 0 if all is OK ;   = -1 otherwise
*.
*. CALLS     : PXSEAR, PXSAME, PXNEW, PXTRY, PXORD, PXUVEC, PXOLAP
*. CALLED    : User
*.
*. AUTHOR    :  L.A. del Pozo
*. CREATED   :  26-Feb-93
*. LAST MOD  :   2-Mar-93
*.
*. Modification Log.
*. 2-Jan-97: M Wobisch    - fix bug concerning COS2R in eta phi mode
*. 4-Apr-93: M H Seymour  - Change 2d arrays to 1d in PXTRY & PXNEW
*. 2-Apr-93: M H Seymour  - Major changes to add boost-invariant mode
*. 1-Apr-93: M H Seymour  - Increase all array sizes
*. 30-Mar-93: M H Seymour - Change all REAL variables to DOUBLE PRECISION
*. 30-Mar-93: M H Seymour - Change OVLIM into an input parameter
*. 2-Mar-93: L A del Pozo - Fix Bugs in PXOLAP
*. 1-Mar-93: L A del Pozo - Remove Cern library routine calls
*. 1-Mar-93: L A del Pozo - Add Print out of welcome and R and Epsilon
*.
*.*********************************************************
C+SEQ,DECLARE.
*** External Arrays
      INTEGER  ITKDM,MXJET,NTRAK,NJET,IERR,MODE
      INTEGER  IPASS (NTRAK),IJMUL (MXJET)
      DOUBLE PRECISION  PTRAK (ITKDM,NTRAK),PJET (5,MXJET), 
     +     CONER, EPSLON, OVLIM
*** Internal Arrays
      INTEGER MXPROT, MXTRAK
      PARAMETER (MXPROT=5000, MXTRAK=5000)
      DOUBLE PRECISION PP(4,MXTRAK), PU(3,MXTRAK), PJ(4,MXPROT)
      LOGICAL JETLIS(MXPROT,MXTRAK)
*** Used in the routine.
      DOUBLE PRECISION COSR,COS2R, VSEED(3), VEC1(3), VEC2(3),PTSQ,PPSQ,
     +     COSVAL,PXMDPI
cMWobisch
      DOUBLE PRECISION RSEP
cMWobisch
      LOGICAL UNSTBL
      INTEGER I,J,N,MU,N1,N2, ITERR
      INTEGER NCALL, NPRINT
      DOUBLE PRECISION ROLD, EPSOLD, OVOLD
      SAVE NCALL,NPRINT,ROLD, EPSOLD, OVOLD
      DATA NCALL,NPRINT /0,0/
      DATA ROLD,EPSOLD,OVOLD/0.,0.,0./

cMWobisch
c***************************************
      RSEP  = 2D0
c***************************************
cMWobisch
      IERR=0
*
*** INITIALIZE
      IF(NCALL.LE.0)  THEN
         ROLD = 0.
         EPSOLD = 0.
         OVOLD = 0.
      ENDIF
      NCALL = NCALL + 1
*
*** Print welcome and Jetfinder parameters
      IF((CONER.NE.ROLD .OR. EPSLON.NE.EPSOLD .OR. OVLIM.NE.OVOLD)
     +     .AND. NPRINT.LE.10) THEN
         WRITE (6,*)
         WRITE (6,*) ' *********** PXCONE: Cone Jet-finder ***********'
         WRITE (6,*) '    Written by Luis Del Pozo of OPAL'
         WRITE (6,*) '    Modified for eta-phi by Mike Seymour'
         WRITE(6,1000)'   Cone Size R = ',CONER,' Radians'
         WRITE(6,1001)'   Min Jet energy Epsilon = ',EPSLON,' GeV'
         WRITE(6,1002)'   Overlap fraction parameter = ',OVLIM
         WRITE (6,*) ' ***********************************************'
cMWobisch
         IF (RSEP .lt. 1.999) THEN
            WRITE(6,*) ' '
            WRITE (6,*) ' ******************************************'
            WRITE (6,*) ' ******************************************'
            WRITE(6,*) ' M Wobisch: private change !!!!!!!!!!!! '
            WRITE(6,*) '      Rsep is set to ',RSEP
            WRITE(6,*) ' this is ONLY meaningful in a NLO calculation'
            WRITE(6,*) '      ------------------------  '
            WRITE(6,*) '  please check what you''re doing!!'
            WRITE(6,*) '   or ask:  Markus.Wobisch@desy.de --'
            WRITE (6,*) ' ******************************************'
            WRITE (6,*) ' ******************************************'
            WRITE (6,*) ' ******************************************'
            WRITE(6,*) ' '
            WRITE(6,*) ' '
         ENDIF
cMWobisch

         WRITE (6,*)
1000     FORMAT(A18,F5.2,A10)
1001     FORMAT(A29,F5.2,A5)
1002     FORMAT(A33,F5.2)
         NPRINT = NPRINT + 1
         ROLD=CONER
         EPSOLD=EPSLON
         OVOLD=OVLIM
      ENDIF
*
*** Copy calling array PTRAK  to internal array PP(4,NTRAK)
*
      IF (NTRAK .GT. MXTRAK) THEN
         WRITE (6,*) ' PXCONE: Ntrak too large: ',NTRAK
         IERR=-1
         RETURN
      ENDIF
      IF (MODE.NE.2) THEN
         DO  100 I=1, NTRAK
            DO  101 J=1,4
               PP(J,I)=PTRAK(J,I)
101         CONTINUE
100      CONTINUE
      ELSE
*** Converting to eta,phi,pt if necessary
         DO  104 I=1,NTRAK
            PTSQ=PTRAK(1,I)**2+PTRAK(2,I)**2
            PPSQ=(SQRT(PTSQ+PTRAK(3,I)**2)+ABS(PTRAK(3,I)))**2
            IF (PTSQ.LE.4.25E-18*PPSQ) THEN
               PP(1,I)=20
            ELSE
               PP(1,I)=0.5*LOG(PPSQ/PTSQ)
            ENDIF
            PP(1,I)=SIGN(PP(1,I),PTRAK(3,I))
            IF (PTSQ.EQ.0) THEN
               PP(2,I)=0
            ELSE
               PP(2,I)=ATAN2(PTRAK(2,I),PTRAK(1,I))
            ENDIF
            PP(3,I)=0
            PP(4,I)=SQRT(PTSQ)
            PU(1,I)=PP(1,I)
            PU(2,I)=PP(2,I)
            PU(3,I)=PP(3,I)
104      CONTINUE
      ENDIF
*
*** Zero output variables
*
      NJET=0
      DO 102 I = 1, NTRAK
         DO 103 J = 1, MXPROT
           JETLIS(J,I) = .FALSE.
103      CONTINUE
102   CONTINUE
      CALL PXZERV(4*MXPROT,PJ)
      CALL PXZERI(MXJET,IJMUL)
*
      IF (MODE.NE.2) THEN
         COSR = COS(CONER)
         COS2R = COS(CONER)
      ELSE
*** Purely for convenience, work in terms of 1-R**2
         COSR = 1-CONER**2
cMW -- select Rsep: 1-(Rsep*CONER)**2
         COS2R =  1-(RSEP*CONER)**2
cORIGINAL         COS2R =  1-(2*CONER)**2
      ENDIF
      UNSTBL = .FALSE.
      IF (MODE.NE.2) THEN
         CALL PXUVEC(NTRAK,PP,PU,IERR)
         IF (IERR .NE. 0) RETURN
      ENDIF
*** Look for jets using particle diretions as seed axes
*
      DO 110 N = 1,NTRAK
        DO 120 MU = 1,3
          VSEED(MU) = PU(MU,N)
120     CONTINUE
        CALL PXSEAR(MODE,COSR,NTRAK,PU,PP,VSEED,
     &                   NJET,JETLIS,PJ,UNSTBL,IERR)
         IF (IERR .NE. 0) RETURN
110   CONTINUE

cMW - for Rsep=1 goto 145
c      GOTO 145

*** Now look between all pairs of jets as seed axes.
      DO 140 N1 = 1,NJET-1
         VEC1(1)=PJ(1,N1)
         VEC1(2)=PJ(2,N1)
         VEC1(3)=PJ(3,N1)
         IF (MODE.NE.2) CALL PXNORV(3,VEC1,VEC1,ITERR)
         DO 150 N2 = N1+1,NJET
            VEC2(1)=PJ(1,N2)
            VEC2(2)=PJ(2,N2)
            VEC2(3)=PJ(3,N2)
            IF (MODE.NE.2) CALL PXNORV(3,VEC2,VEC2,ITERR)
            CALL PXADDV(3,VEC1,VEC2,VSEED,ITERR)
            IF (MODE.NE.2) THEN
               CALL PXNORV(3,VSEED,VSEED,ITERR)
            ELSE
               VSEED(1)=VSEED(1)/2
               VSEED(2)=VSEED(2)/2
            ENDIF
C---ONLY BOTHER IF THEY ARE BETWEEN 1 AND 2 CONE RADII APART
            IF (MODE.NE.2) THEN
              COSVAL=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
            ELSE
               IF (ABS(VEC1(1)).GE.20.OR.ABS(VEC2(1)).GE.20) THEN
                  COSVAL=-1000
               ELSE
                  COSVAL=1-
     +              ((VEC1(1)-VEC2(1))**2+PXMDPI(VEC1(2)-VEC2(2))**2)
               ENDIF
            ENDIF

            IF (COSVAL.LE.COSR.AND.COSVAL.GE.COS2R)
     +           CALL PXSEAR(MODE,COSR,NTRAK,PU,PP,VSEED,NJET,
     +           JETLIS,PJ,UNSTBL,IERR)
c            CALL PXSEAR(MODE,COSR,NTRAK,PU,PP,VSEED,NJET,
c     +           JETLIS,PJ,UNSTBL,IERR)
            IF (IERR .NE. 0) RETURN
150      CONTINUE
140   CONTINUE
      IF (UNSTBL) THEN
        IERR=-1
        WRITE (6,*) ' PXCONE: Too many iterations to find a proto-jet'
        RETURN
      ENDIF

 145  CONTINUE
*** Now put the jet list into order by jet energy, eliminating jets
*** with energy less than EPSLON.
       CALL PXORD(EPSLON,NJET,NTRAK,JETLIS,PJ)
*
*** Take care of jet overlaps
       CALL PXOLAP(MODE,NJET,NTRAK,JETLIS,PJ,PP,OVLIM)
*
*** Order jets again as some have been eliminated, or lost energy.
       CALL PXORD(EPSLON,NJET,NTRAK,JETLIS,PJ)
*
*** All done!, Copy output into output arrays
      IF (NJET .GT. MXJET) THEN
         WRITE (6,*) ' PXCONE:  Found more than MXJET jets'
         IERR=-1
         GOTO 99
      ENDIF
      IF (MODE.NE.2) THEN
         DO 300 I=1, NJET
            DO 310 J=1,4
               PJET(J,I)=PJ(J,I)
310         CONTINUE
300      CONTINUE
      ELSE
         DO 315 I=1, NJET
            PJET(1,I)=PJ(4,I)*COS(PJ(2,I))
            PJET(2,I)=PJ(4,I)*SIN(PJ(2,I))
            PJET(3,I)=PJ(4,I)*SINH(PJ(1,I))
            PJET(4,I)=PJ(4,I)*COSH(PJ(1,I))
 315     CONTINUE
      ENDIF
      DO 320 I=1, NTRAK
         IPASS(I)=-1
         DO 330 J=1, NJET
            IF (JETLIS(J,I)) THEN
               IJMUL(J)=IJMUL(J)+1
               IPASS(I)=J
            ENDIF
330      CONTINUE
320   CONTINUE
99    RETURN
      END
*CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper
*-- Author :
C-----------------------------------------------------------------------
      SUBROUTINE PXNORV(N,A,B,ITERR)
      INTEGER I,N,ITERR
      DOUBLE PRECISION A(N),B(N),C
      C=0
      DO 10 I=1,N
        C=C+A(I)**2
 10   CONTINUE
      IF (C.LE.0) RETURN
      C=1/SQRT(C)
      DO 20 I=1,N
        B(I)=A(I)*C
 20   CONTINUE
      END
*CMZ :  2.00/00 10/01/95  10.17.57  by  P. Schleper
*CMZ :  1.06/00 15/03/94  12.17.46  by  P. Schleper
*-- Author :
*
C+DECK,PXOLAP.
      SUBROUTINE PXOLAP(MODE,NJET,NTRAK,JETLIS,PJ,PP,OVLIM)
*
*** Looks for particles assigned to more than 1 jet, and reassigns them
*** If more than a fraction OVLIM of a jet's energy is contained in
*** higher energy jets, that jet is neglected.
*** Particles assigned to the jet closest in angle (a la CDF, Snowmass).
C+SEQ,DECLARE.
      INTEGER MXTRAK, MXPROT
      PARAMETER (MXTRAK=5000,MXPROT=5000)
      INTEGER NJET, NTRAK, MODE
      LOGICAL JETLIS(MXPROT,MXTRAK)
      DOUBLE PRECISION PJ(4,MXPROT),PP(4,MXTRAK),PXMDPI
      INTEGER I,J,N,MU
      LOGICAL OVELAP
      DOUBLE PRECISION EOVER
      DOUBLE PRECISION OVLIM
      INTEGER ITERR, IJMIN, IJET(MXPROT), NJ
      DOUBLE PRECISION VEC1(3), VEC2(3), COST, THET, THMIN
      DATA IJMIN/0/
*
      IF (NJET.LE.1) RETURN
*** Look for jets with large overlaps with higher energy jets.
      DO 100 I = 2,NJET
*** Find overlap energy between jets I and all higher energy jets.
       EOVER = 0.0
       DO 110 N = 1,NTRAK
         OVELAP = .FALSE.
         DO 120 J= 1,I-1
           IF (JETLIS(I,N).AND.JETLIS(J,N)) THEN
            OVELAP = .TRUE.
           ENDIF
120      CONTINUE
         IF (OVELAP) THEN
           EOVER = EOVER + PP(4,N)
         ENDIF
110     CONTINUE
*** Is the fraction of energy shared larger than OVLIM?
        IF (EOVER.GT.OVLIM*PJ(4,I)) THEN
*** De-assign all particles from Jet I
            DO 130 N = 1,NTRAK
              JETLIS(I,N) = .FALSE.
130         CONTINUE
         ENDIF
100   CONTINUE
*** Now there are no big overlaps, assign every particle in
*** more than 1 jet to the closet jet.
*** Any particles now in more than 1 jet are assigned to the CLOSET
*** jet (in angle).
      DO 140 I=1,NTRAK
         NJ=0
         DO 150 J=1, NJET
         IF(JETLIS(J,I)) THEN
            NJ=NJ+1
            IJET(NJ)=J
         ENDIF
150      CONTINUE
         IF (NJ .GT. 1) THEN
*** Particle in > 1 jet - calc angles...
            VEC1(1)=PP(1,I)
            VEC1(2)=PP(2,I)
            VEC1(3)=PP(3,I)
            THMIN=0.
            DO 160 J=1,NJ
               VEC2(1)=PJ(1,IJET(J))
               VEC2(2)=PJ(2,IJET(J))
               VEC2(3)=PJ(3,IJET(J))
               IF (MODE.NE.2) THEN
                  CALL PXANG3(VEC1,VEC2,COST,THET,ITERR)
               ELSE
                  THET=(VEC1(1)-VEC2(1))**2+PXMDPI(VEC1(2)-VEC2(2))**2
               ENDIF
               IF (J .EQ. 1) THEN
                  THMIN=THET
                  IJMIN=IJET(J)
               ELSEIF (THET .LT. THMIN) THEN
                  THMIN=THET
                  IJMIN=IJET(J)
               ENDIF
160         CONTINUE
*** Assign track to IJMIN
            DO 170 J=1,NJET
               JETLIS(J,I) = .FALSE.
170         CONTINUE
            JETLIS(IJMIN,I)=.TRUE.
         ENDIF
140   CONTINUE
*** Recompute PJ
      DO 200 I = 1,NJET
        DO 210 MU = 1,4
          PJ(MU,I) = 0.0
210     CONTINUE
        DO 220 N = 1,NTRAK
          IF( JETLIS(I,N) ) THEN
             IF (MODE.NE.2) THEN
                DO 230 MU = 1,4
                   PJ(MU,I) = PJ(MU,I) + PP(MU,N)
230             CONTINUE
             ELSE
                PJ(1,I)=PJ(1,I)
     +               + PP(4,N)/(PP(4,N)+PJ(4,I))*(PP(1,N)-PJ(1,I))
                PJ(2,I)=PJ(2,I)
     +               + PP(4,N)/(PP(4,N)+PJ(4,I))*PXMDPI(PP(2,N)-PJ(2,I))
                PJ(4,I)=PJ(4,I)+PP(4,N)
             ENDIF
          ENDIF
220     CONTINUE
200   CONTINUE
      RETURN
      END
*CMZ :  2.00/00 10/01/95  10.17.57  by  P. Schleper
*CMZ :  1.06/00 14/03/94  15.37.45  by  P. Schleper
*-- Author :
*
C+DECK,PXORD.
       SUBROUTINE PXORD(EPSLON,NJET,NTRAK,JETLIS,PJ)
*
*** Routine to put jets into order and eliminate tose less than EPSLON
C+SEQ,DECLARE.
      INTEGER MXTRAK,MXPROT
      PARAMETER (MXTRAK=5000,MXPROT=5000)
       INTEGER I, J, INDEX(MXPROT)
       DOUBLE PRECISION PTEMP(4,MXPROT), ELIST(MXPROT)
       INTEGER NJET,NTRAK
       LOGICAL JETLIS(MXPROT,MXTRAK)
       LOGICAL LOGTMP(MXPROT,MXTRAK)
       DOUBLE PRECISION EPSLON,PJ(4,MXPROT)
*** Puts jets in order of energy: 1 = highest energy etc.
*** Then Eliminate jets with energy below EPSLON
*
*** Copy input arrays.
      DO 100 I=1,NJET
         DO 110 J=1,4
            PTEMP(J,I)=PJ(J,I)
110      CONTINUE
         DO 120 J=1,NTRAK
            LOGTMP(I,J)=JETLIS(I,J)
120      CONTINUE
100   CONTINUE
      DO 150 I=1,NJET
         ELIST(I)=PJ(4,I)
150   CONTINUE
*** Sort the energies...
      CALL PXSORV(NJET,ELIST,INDEX,'I')
*** Fill PJ and JETLIS according to sort ( sort is in ascending order!!)
      DO 200 I=1, NJET
         DO 210 J=1,4
            PJ(J,I)=PTEMP(J,INDEX(NJET+1-I))
210      CONTINUE
         DO 220 J=1,NTRAK
            JETLIS(I,J)=LOGTMP(INDEX(NJET+1-I),J)
220      CONTINUE
200   CONTINUE
** Jets are now in order
*** Now eliminate jets with less than Epsilon energy
      DO 300, I=1, NJET
         IF (PJ(4,I) .LT. EPSLON) THEN
            NJET=NJET-1
            PJ(4,I)=0.
         ENDIF
300   CONTINUE
      RETURN
      END

********************************************************************
*CMZ :  2.00/00 10/01/95  10.17.57  by  P. Schleper
*CMZ :  1.06/00 14/03/94  15.37.44  by  P. Schleper
*-- Author :
C+DECK,PXSEAR.
      SUBROUTINE PXSEAR(MODE,COSR,NTRAK,PU,PP,VSEED,NJET,
     +                JETLIS,PJ,UNSTBL,IERR)
*
C+SEQ,DECLARE.
      INTEGER MXTRAK, MXPROT
      PARAMETER (MXTRAK=5000,MXPROT=5000)
      INTEGER NTRAK, IERR, MODE
      DOUBLE PRECISION COSR,PU(3,MXTRAK),PP(4,MXTRAK),VSEED(3)
      LOGICAL UNSTBL
      LOGICAL JETLIS(MXPROT,MXTRAK)
      INTEGER NJET
      DOUBLE PRECISION  PJ(4,MXPROT)
*** Using VSEED as a trial axis , look for a stable jet.
*** Check stable jets against those already found and add to PJ.
*** Will try up to MXITER iterations to get a stable set of particles
*** in the cone.
      INTEGER MU,N,ITER
      LOGICAL PXSAME,PXNEW,OK
      LOGICAL NEWLIS(MXTRAK),OLDLIS(MXTRAK)
      DOUBLE PRECISION OAXIS(3),NAXIS(3),PNEW(4)
      INTEGER MXITER
      PARAMETER(MXITER = 30)
*
      DO 100 MU=1,3
        OAXIS(MU) = VSEED(MU)
100   CONTINUE
      DO 110 N = 1,NTRAK
        OLDLIS(N) = .FALSE.
110   CONTINUE
      DO 120 ITER = 1,MXITER
        CALL PXTRY(MODE,COSR,NTRAK,PU,PP,OAXIS,NAXIS,PNEW,NEWLIS,OK)
*** Return immediately if there were no particles in the cone.
       IF (.NOT.OK) THEN
         RETURN
       ENDIF
       IF(PXSAME(NEWLIS,OLDLIS,NTRAK)) THEN
*** We have a stable jet.
             IF (PXNEW(NEWLIS,JETLIS,NTRAK,NJET)) THEN
*** And the jet is a new one. So add it to our arrays.
*** Check arrays are big anough...
             IF (NJET .EQ. MXPROT) THEN
             WRITE (6,*) ' PXCONE:  Found more than MXPROT proto-jets'
                IERR = -1
                RETURN
             ENDIF
               NJET = NJET + 1
               DO 130 N = 1,NTRAK
                 JETLIS(NJET,N) = NEWLIS(N)
130            CONTINUE
               DO 140 MU=1,4
                 PJ(MU,NJET)=PNEW(MU)
140          CONTINUE
             ENDIF
             RETURN
       ENDIF
*** The jet was not stable, so we iterate again
       DO 150 N=1,NTRAK
         OLDLIS(N)=NEWLIS(N)
150    CONTINUE
       DO 160 MU=1,3
         OAXIS(MU)=NAXIS(MU)
160    CONTINUE
120   CONTINUE
      UNSTBL = .TRUE.
      RETURN
      END
*CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper
*-- Author :
C-----------------------------------------------------------------------
      SUBROUTINE PXSORV(N,A,K,OPT)
C     Sort A(N) into ascending order
C     OPT = 'I' : return index array K only
C     OTHERWISE : return sorted A and index array K
C-----------------------------------------------------------------------
      INTEGER NMAX
      PARAMETER (NMAX=5000)
*
*      INTEGER N,I,J,K(N),IL(NMAX),IR(NMAX)
*LUND
      INTEGER N,I,J,K(NMAX),IL(NMAX),IR(NMAX)
      CHARACTER OPT
*
*      DOUBLE PRECISION A(N),B(NMAX)
      DOUBLE PRECISION A(NMAX),B(NMAX)
*LUND
      IF (N.GT.NMAX) STOP 'Sorry, not enough room in Mike''s PXSORV'
      IL(1)=0
      IR(1)=0
      DO 10 I=2,N
      IL(I)=0
      IR(I)=0
      J=1
   2  IF(A(I).GT.A(J)) GO TO 5
   3  IF(IL(J).EQ.0) GO TO 4
      J=IL(J)
      GO TO 2
   4  IR(I)=-J
      IL(J)=I
      GO TO 10
   5  IF(IR(J).LE.0) GO TO 6
      J=IR(J)
      GO TO 2
   6  IR(I)=IR(J)
      IR(J)=I
  10  CONTINUE
      I=1
      J=1
      GO TO 8
  20  J=IL(J)
   8  IF(IL(J).GT.0) GO TO 20
   9  K(I)=J
      B(I)=A(J)
      I=I+1
      IF(IR(J)) 12,30,13
  13  J=IR(J)
      GO TO 8
  12  J=-IR(J)
      GO TO 9
  30  IF(OPT.EQ.'I') RETURN
      DO 31 I=1,N
  31  A(I)=B(I)
 999  END

*********************************************************************
*CMZ :  2.00/00 10/01/95  10.17.57  by  P. Schleper
*CMZ :  1.06/00 14/03/94  15.37.44  by  P. Schleper
*-- Author :
*
C+DECK,PXTRY.
       SUBROUTINE PXTRY(MODE,COSR,NTRAK,PU,PP,OAXIS,NAXIS,
     +                  PNEW,NEWLIS,OK)
*
C+SEQ,DECLARE.
      INTEGER MXTRAK
      PARAMETER (MXTRAK=5000)
       INTEGER NTRAK,MODE
*** Note that although PU and PP are assumed to be 2d arrays, they
*** are used as 1d in this routine for efficiency
       DOUBLE PRECISION COSR,PU(3*MXTRAK),PP(4*MXTRAK),OAXIS(3),PXMDPI
       LOGICAL OK
       LOGICAL NEWLIS(MXTRAK)
       DOUBLE PRECISION NAXIS(3),PNEW(4)
*** Finds all particles in cone of size COSR about OAXIS direction.
*** Calculates 4-momentum sum of all particles in cone (PNEW) , and
*** returns this as new jet axis NAXIS (Both unit Vectors)
       INTEGER N,MU,NPU,NPP
       DOUBLE PRECISION COSVAL,NORMSQ,NORM
*
       OK = .FALSE.
       DO 100 MU=1,4
          PNEW(MU)=0.0
100    CONTINUE
       NPU=-3
       NPP=-4
       DO 110 N=1,NTRAK
          NPU=NPU+3
          NPP=NPP+4
          IF (MODE.NE.2) THEN
             COSVAL=0.0
             DO 120 MU=1,3
                COSVAL=COSVAL+OAXIS(MU)*PU(MU+NPU)
120          CONTINUE
          ELSE
             IF (ABS(PU(1+NPU)).GE.20.OR.ABS(OAXIS(1)).GE.20) THEN
                COSVAL=-1000
             ELSE
                COSVAL=1-
     +           ((OAXIS(1)-PU(1+NPU))**2+PXMDPI(OAXIS(2)-PU(2+NPU))**2)
             ENDIF
          ENDIF
          IF (COSVAL.GE.COSR)THEN
             NEWLIS(N) = .TRUE.
             OK = .TRUE.
             IF (MODE.NE.2) THEN
                DO 130 MU=1,4
                   PNEW(MU) = PNEW(MU) + PP(MU+NPP)
130             CONTINUE
             ELSE
                PNEW(1)=PNEW(1)
     +              + PP(4+NPP)/(PP(4+NPP)+PNEW(4))*(PP(1+NPP)-PNEW(1))
                PNEW(2)=PNEW(2)
     +              + PP(4+NPP)/(PP(4+NPP)+PNEW(4))
     +               *PXMDPI(PP(2+NPP)-PNEW(2))
                PNEW(4)=PNEW(4)+PP(4+NPP)
             ENDIF
          ELSE
             NEWLIS(N)=.FALSE.
          ENDIF
110   CONTINUE
*** If there are particles in the cone, calc new jet axis
       IF (OK) THEN
          IF (MODE.NE.2) THEN
             NORMSQ = 0.0
             DO 140 MU = 1,3
                NORMSQ = NORMSQ + PNEW(MU)**2
140          CONTINUE
             NORM = SQRT(NORMSQ)
          ELSE
             NORM = 1
          ENDIF
          DO 150 MU=1,3
             NAXIS(MU) = PNEW(MU)/NORM
150       CONTINUE
       ENDIF
       RETURN
       END

*********************************************************************
*CMZ :  2.00/00 10/01/95  10.17.57  by  P. Schleper
*CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper
*-- Author :
C+DECK,PXUVEC.
*
       SUBROUTINE PXUVEC(NTRAK,PP,PU,IERR)
*
*** Routine to calculate unit vectors PU of all particles PP
C+SEQ,DECLARE.
      INTEGER MXTRAK
      PARAMETER (MXTRAK=5000)
      INTEGER NTRAK, IERR
      DOUBLE PRECISION PP(4,MXTRAK)
      DOUBLE PRECISION PU(3,MXTRAK)
      INTEGER N,MU
      DOUBLE PRECISION MAG
       DO 100 N=1,NTRAK
          MAG=0.0
          DO 110 MU=1,3
             MAG=MAG+PP(MU,N)**2
110       CONTINUE
          MAG=SQRT(MAG)
          IF (MAG.EQ.0.0) THEN
             WRITE(6,*)' PXCONE: An input particle has zero mod(p)'
             IERR=-1
             RETURN
          ENDIF
          DO 120 MU=1,3
           PU(MU,N)=PP(MU,N)/MAG
120       CONTINUE
100    CONTINUE
       RETURN
       END
*CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper
*-- Author :
C-----------------------------------------------------------------------
      SUBROUTINE PXZERI(N,A)
      INTEGER I,N,A(N)
      DO 10 I=1,N
        A(I)=0
 10   CONTINUE
      END
*CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper
*-- Author :
C-----------------------------------------------------------------------
C     This is a set of routines written by Mike Seymour to provide the
C     services presumably normally provided by standard OPAL routines
C     PXZERV zeroes a vector
C     PXZERI zeroes a vector of integers
C     PXNORV normalizes a vector
C     PXADDV adds two vectors
C     PXSORV sorts a vector (copied from HERWIG)
C     PXANG3 finds the angle (and its cosine) between two vectors
C     PXMDPI moves its argument onto the range [-pi,pi)
C-----------------------------------------------------------------------
      SUBROUTINE PXZERV(N,A)
      INTEGER I,N
      DOUBLE PRECISION A(N)
      DO 10 I=1,N
        A(I)=0
 10   CONTINUE
      END
*-- Author :
C-----------------------------------------------------------------------
      SUBROUTINE PXADDV(N,A,B,C,ITERR)
      INTEGER I,N,ITERR
      DOUBLE PRECISION A(N),B(N),C(N)
      DO 10 I=1,N
        C(I)=A(I)+B(I)
 10   CONTINUE
      END
*CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper
*-- Author :
C-----------------------------------------------------------------------
      SUBROUTINE PXANG3(A,B,COST,THET,ITERR)
      INTEGER ITERR
      DOUBLE PRECISION A(3),B(3),C,COST,THET
      C=(A(1)**2+A(2)**2+A(3)**2)*(B(1)**2+B(2)**2+B(3)**2)
      IF (C.LE.0) RETURN
      C=1/SQRT(C)
      COST=(A(1)*B(1)+A(2)*B(2)+A(3)*B(3))*C
      THET=ACOS(COST)
      END
*CMZ :  1.06/00 14/03/94  15.41.57  by  P. Schleper
*-- Author :    P. Schleper   28/02/94
       LOGICAL FUNCTION PXNEW(TSTLIS,JETLIS,NTRAK,NJET)
*
      INTEGER MXTRAK,MXPROT
      PARAMETER (MXTRAK=5000,MXPROT=5000)
       INTEGER NTRAK,NJET
*** Note that although JETLIS is assumed to be a 2d array, it
*** it is used as 1d in this routine for efficiency
       LOGICAL TSTLIS(MXTRAK),JETLIS(MXPROT*MXTRAK)
*** Checks to see if TSTLIS entries correspond to a jet already found
*** and entered in JETLIS
       INTEGER N, I, IN
       LOGICAL MATCH
*
       PXNEW = .TRUE.
       DO 100 I = 1,NJET
          MATCH = .TRUE.
          IN=I-MXPROT
          DO 110 N = 1,NTRAK
            IN=IN+MXPROT
            IF(TSTLIS(N).NEQV.JETLIS(IN)) THEN
             MATCH = .FALSE.
             GO TO 100
            ENDIF
110       CONTINUE
          IF (MATCH) THEN
           PXNEW = .FALSE.
           RETURN
          ENDIF
100    CONTINUE
       RETURN
       END
*CMZ :  1.06/00 14/03/94  15.41.57  by  P. Schleper
*-- Author :    P. Schleper   28/02/94
       LOGICAL FUNCTION PXSAME(LIST1,LIST2,N)
*
       LOGICAL LIST1(*),LIST2(*)
       INTEGER N
*** Returns T if the first N elements of LIST1 are the same as the
*** first N elements of LIST2.
       INTEGER I
*
       PXSAME = .TRUE.
       DO 100 I = 1,N
        IF ( LIST1(I).NEQV.LIST2(I) ) THEN
          PXSAME = .FALSE.
          RETURN
        ENDIF
100    CONTINUE
       RETURN
       END
*CMZ :  1.06/00 28/02/94  15.44.44  by  P. Schleper
*-- Author :
C-----------------------------------------------------------------------
      FUNCTION PXMDPI(PHI)
      IMPLICIT NONE
C---RETURNS PHI, MOVED ONTO THE RANGE [-PI,PI)
      DOUBLE PRECISION PXMDPI,PHI,PI,TWOPI,THRPI,EPS
      PARAMETER (PI=3.141592654,TWOPI=6.283185307,THRPI=9.424777961)
      PARAMETER (EPS=1E-15)
      PXMDPI=PHI
      IF (PXMDPI.LE.PI) THEN
        IF (PXMDPI.GT.-PI) THEN
          GOTO 100
        ELSEIF (PXMDPI.GT.-THRPI) THEN
          PXMDPI=PXMDPI+TWOPI
        ELSE
          PXMDPI=-MOD(PI-PXMDPI,TWOPI)+PI
        ENDIF
      ELSEIF (PXMDPI.LE.THRPI) THEN
        PXMDPI=PXMDPI-TWOPI
      ELSE
        PXMDPI=MOD(PI+PXMDPI,TWOPI)-PI
      ENDIF
 100  IF (ABS(PXMDPI).LT.EPS) PXMDPI=0
      END
