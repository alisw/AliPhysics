*$ CREATE COMSCW.FOR
*COPY COMSCW
*
*=== comscw  ===========================================================*
*
      DOUBLE PRECISION FUNCTION COMSCW ( IJ    , XA    , YA    , ZA    ,
     &                                   MREG  , RULL  , LLO   , ICALL )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Special multiplication factors for activity scoring.                 *
*                                                                      *
* Note: Data files 'LE-CH.dat' and 'LE-EC.dat' read on unit 20.        *
*                                                                      *
*                                                                      *
* SDUM = 'abcdefgh'                                              Lesco *
*                                                                      *
*   'ab' = 'LE'  division by exemption limits activated                *
*                                                                      *
*        'c' = 'C' : division by CERN zonage limits                    *
*            'def' = '10 '    : division by 1/10 of the value      10  *
*            'def' = '100'    : division by 1/100 of the value    100  *
*                    otherwise: division by the value itself        1  *
*        'c' = 'S' : division by Swiss LE values                    2  *
*        'c' = 'E' : division by European concentration limits      3  *
*                                                                      *
*                Comscw = 1. for selected isotope                      *
*                       = 0. otherwise                                 *
*                                                                      *
*   'ab' = 'IS'  selection of individual isotops activated       9999  *
*                                                                      *
*        'cd'  = symbol of isotope                                     *
*        'efg' = mass number of isotope                                *
*        'h'   = 'm' : metastable state                                *
*                                                                      *
*                Comscw = 1. for selected isotope                      *
*                       = 0. otherwise                                 *
*                                                                      *
*   otherwise    Comscw = 1.                                       -1  *
*                                                                      *
*   Example(s):                                                        *
*                                                                      *
*    SDUM = 'LEC10xxx' : division by 1/10 of LE(EU) or by Swiss LE     *
*         = 'LESxxxxx' : division by Swiss LE                          *
*         = 'ISV  48x' : scoring of activity for V48                   *
*         = 'ISMn 52m' : scoring of activity for mMn52                 *
*                                                                      *
*    Note: 'x' refers to characters without specific meaning           *
*                                                                      *
*     Input variables:                                                 *
*                                                                      *
*           Ij = (generalized) particle code                           *
*     Xa,Ya,Za = position                                              *
*         Mreg = region number                                         *
*         Rull = amount to be deposited                                *
*          Llo = particle generation                                   *
*        Icall = call id                                               *
*                                                                      *
*     Output variables:                                                *
*                                                                      *
*       Comscw = factor the scored amount will be multiplied by        *
*       Lsczer = logical flag, if true no amount will be scored        *
*                regardless of Comscw                                  *
*                                                                      *
*     Useful variables (common SCOHLP):                                *
*                                                                      *
*     Energy/Star binnings/scorings (Comscw):                          *
*          ISCRNG = 1 --> Energy density  binning                      *
*          ISCRNG = 2 --> Star   density  binning                      *
*          ISCRNG = 3 --> Residual nuclei scoring                      *
*          ISCRNG = 4 --> Momentum transfer density binning            *
*          ISCRNG = 5 --> Activity density binning                     *
*          JSCRNG = # of the binning                                   *
*                                                                      *
*     Useful variables (common SOUEVT):                                *
*                                                                      *
*          X,Y,Zsoevt(i) = position    of the i_th source particle     *
*          TX,Y,Zsoev(i) = direction   of the i_th source particle     *
*              Wtsoev(i) = weight      of the i_th source particle     *
*              Pmsoev(i) = momentum    of the i_th source particle     *
*              Tksoev(i) = kin. energy of the i_th source particle     *
*              Agsoev(i) = age         of the i_th source particle     *
*              Aksoev(i) = Kaon ampl.  of the i_th source particle     *
*              Ussoev(i) = user var.   of the i_th source particle     *
*              Ijsoev(i) = identity    of the i_th source particle     *
*              Nrsoev(i) = region      of the i_th source particle     *
*              Nlsoev(i) = lattice     of the i_th source particle     *
*                Npsoev  = number of the source particles              *
*----------------------------------------------------------------------*
*
      INCLUDE '(FLKMAT)'
      INCLUDE '(SCOHLP)'
      INCLUDE '(SOUEVT)'
*
      INCLUDE '(RSNCCM)'
      INCLUDE '(USRBIN)'
      DIMENSION IZSCO(MXUSBN),IASCO(MXUSBN),ISSCO(MXUSBN),
     &          LESCO(MXUSBN)
*
      CHARACTER CISOIN*2,CISO*2,CSET*10,CDUM*4,CA*3
      PARAMETER (IZMAX  = 109,
     &           IAZMIN = 2,
     &           IAZMAX = 160)
      DIMENSION CISO(IZMAX)
      DIMENSION XLESWS(IZMAX,IAZMIN:IAZMAX,2),
     &          XLEECO(IZMAX,IAZMIN:IAZMAX,2),
     &          XLEO10(IZMAX,IAZMIN:IAZMAX,2)
      DATA CISO /
     &    'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na',
     &    'Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti',
     &    'V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As',
     &    'Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru',
     &    'Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs',
     &    'Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
     &    'Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir',
     &    'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
     &    'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es',
     &    'Fm','Md','No','Lr','Rf','Ha','Sg','Ns','Hs','Mt'/
*
      LOGICAL LFIRST,LFOUND
      DATA LFIRST /.TRUE./
*
      LSCZER = .FALSE.
      COMSCW = ONEONE

      IF (LFIRST) THEN

         WRITE(LUNOUT,'(A)') ' COMSCW: activity weighting activated'
         LFIRST = .FALSE.
*
*----------------------------------------------------------------------*
* Initialization                               
*
* LE data taken from Appendix 3 column 9,
* Ordonnance sur la radioprotection (ORaP) du 22 juin 1994 
* (etat au 4 avril 2000)
*
* (data file in Bq/kg!)
*
         DO 104 IZ=1,IZMAX
            DO 105 IAZ=IAZMIN,IAZMAX
               XLESWS(IZ,IAZ,1) = ZERZER
               XLESWS(IZ,IAZ,2) = ZERZER
  105       CONTINUE
  104    CONTINUE
*
         OPEN(20,FILE='LE-CH.dat',STATUS='UNKNOWN')
         NISO   = 0
         IZ0    = 1
  100    CONTINUE
         READ(20,*,END=101) CISOIN,IA,XLIMIT
* convert from Bq/kg into Bq/g
         XLIMIT = 1.0D-3*XLIMIT
         LFOUND = .FALSE.
         DO 102 IZ=IZ0,IZMAX
            IF (CISOIN.EQ.CISO(IZ)) THEN
               IF (IA.LT.0) THEN
                  IS = 2
                  IA = ABS(IA)
               ELSE
                  IS = 1
               ENDIF
               IAZ = IA-IZ
               IF ((IAZ.LT.IAZMIN).OR.(IAZ.GT.IAZMAX)) THEN
                  WRITE(LUNOUT,*)
     &               ' COMSCW: warning! Iaz out of allowed range: ',
     &               IAZ, IAZMIN,IAZMAX
                  STOP
               ENDIF
               IF (XLESWS(IZ,IAZ,IS).GT.ZERZER) THEN
                  WRITE(LUNOUT,*) 
     &               ' COMSCW: warning! two entries for this isotope: ',
     &               CISOIN,IA,XLIMIT
                  STOP
               ELSE
                  XLESWS(IZ,IAZ,IS) = XLIMIT
               ENDIF
               NISO = NISO+1
               IZ0  = IZ
               LFOUND = .TRUE.
               GOTO 103
            ENDIF
  102    CONTINUE
         WRITE(LUNOUT,*) 
     &      ' COMSCW: isotope not recognized: ',CISOIN,IA,XLIMIT
         STOP
  103    CONTINUE
         GOTO 100
  101    CONTINUE
         CLOSE(20)
*
* LE data taken from 
* Official journal of the European Communities L159, 29 June 1996
* Council Directive 96/29/Euratom
*
* (data file in Bq/g!)
*
         DO 204 IZ=1,IZMAX
            DO 205 IAZ=IAZMIN,IAZMAX
               XLEECO(IZ,IAZ,1) = ZERZER
               XLEECO(IZ,IAZ,2) = ZERZER
  205       CONTINUE
  204    CONTINUE
*
         OPEN(20,FILE='LE-EC.dat',STATUS='UNKNOWN')
         NISO   = 0
         IZ0    = 1
  200    CONTINUE
         READ(20,*,END=201) CISOIN,IA,XLIMIT,IFLAG
         LFOUND = .FALSE.
         DO 202 IZ=IZ0,IZMAX
            IF (CISOIN.EQ.CISO(IZ)) THEN
               IF (IA.LT.0) THEN
                  IS = 2
                  IA = ABS(IA)
               ELSE
                  IS = 1
               ENDIF
               IAZ = IA-IZ
               IF ((IAZ.LT.IAZMIN).OR.(IAZ.GT.IAZMAX)) THEN
                  WRITE(LUNOUT,*)
     &               ' COMSCW: warning! Iaz out of allowed range: ',
     &               IAZ, IAZMIN,IAZMAX
                  STOP
               ENDIF
               IF (XLEECO(IZ,IAZ,IS).GT.ZERZER) THEN
                  WRITE(LUNOUT,*) 
     &               ' COMSCW: warning! two entries for this isotope: ',
     &               CISOIN,IA,XLIMIT
                  STOP
               ELSE
                  XLEECO(IZ,IAZ,IS) = XLIMIT
* zero entries with Swiss values
                  IF (IFLAG.EQ.1) XLEECO(IZ,IAZ,IS) = ZERZER
               ENDIF
               NISO = NISO+1
               IZ0  = IZ
               LFOUND = .TRUE.
               GOTO 203
            ENDIF
  202    CONTINUE
         WRITE(LUNOUT,*) 
     &      ' COMSCW: isotope not recognized: ',CISOIN,IA,XLIMIT
         STOP
  203    CONTINUE
         GOTO 200
  201    CONTINUE
         CLOSE(20)
         DO 404 IZ=1,IZMAX
            DO 405 IAZ=IAZMIN,IAZMAX
               DO 406 IS=1,2
                  XLEO10(IZ,IAZ,IS) = ZERZER
                  IF (XLEECO(IZ,IAZ,IS).GT.ZERZER) THEN
                     XLEO10(IZ,IAZ,IS) = XLEECO(IZ,IAZ,IS)/10.0D0
                  ELSE
                     XLEO10(IZ,IAZ,IS) = XLESWS(IZ,IAZ,IS)
                  ENDIF
  406          CONTINUE
  405       CONTINUE
  404    CONTINUE
*
         DO 500 I=1,MXUSBN
            IZSCO(I) = 0
            IASCO(I) = 0
            ISSCO(I) = 0
            LESCO(I) = 0
  500    CONTINUE

      ENDIF
*
*----------------------------------------------------------------------*
* Online weighting                             
*
      IF ( ISCRNG .EQ. 5 ) THEN
*
* determine type of weighting
         IF (LESCO(JSCRNG).EQ.0) THEN
            CSET = TITUSB(JSCRNG)
            IF (CSET(1:2).EQ.'IS') THEN
               LESCO(JSCRNG) = 9999
               DO 700 IZ=1,IZMAX
                  IF (CSET(3:4).EQ.CISO(IZ)) THEN
                     IZSCO(JSCRNG) = IZ
                     READ(CSET,'(A4,A3)') CDUM,CA
                     READ(CA,'(I3)') IASCO(JSCRNG)
                     IF (CSET(8:8).EQ.'m') THEN
                        ISSCO(JSCRNG) = 2
                     ELSE
                        ISSCO(JSCRNG) = 1
                     ENDIF
                  ENDIF
  700          CONTINUE
               IF ((IZSCO(JSCRNG).LE.0).OR.(IASCO(JSCRNG).LE.0)
     &                                 .OR.(ISSCO(JSCRNG).LE.0)) THEN
                  WRITE(LUNOUT,*) ' COMSCW: unknown isotope, Z,A,S = ',
     &               IZSCO(JSCRNG),IASCO(JSCRNG),ISSCO(JSCRNG)
                  STOP
               ENDIF
            ELSEIF(CSET(1:2).EQ.'LE') THEN
               IF (CSET(3:3).EQ.'C') THEN
                  IF (CSET(4:6).EQ.'100') THEN
                     LESCO(JSCRNG) = 100
                  ELSEIF (CSET(4:6).EQ.'10') THEN
                     LESCO(JSCRNG) = 10
                  ELSE
                     LESCO(JSCRNG) = 1
                  ENDIF
               ELSEIF (CSET(3:3).EQ.'S') THEN
                  LESCO(JSCRNG) = 2
               ELSEIF (CSET(3:3).EQ.'E') THEN
                  LESCO(JSCRNG) = 3
               ELSE
                  WRITE(LUNOUT,*) ' COMSCW: unknown LE set ',CSET(3:3)
                  STOP
               ENDIF
            ELSE
               LESCO(JSCRNG) = -1
            ENDIF
            WRITE(LUNOUT,1000) ' COMSCW: scoring ',JSCRNG,CSET,
     &         ' weighted with properties ',
     &         LESCO(JSCRNG),IZSCO(JSCRNG),IASCO(JSCRNG),ISSCO(JSCRNG)
 1000       FORMAT(A,I3,2A,I5,3I4)
         ENDIF
*
* obtain present isotope from common block
         JA  = IARSDL(1)
         JZ  = IZRSDL(1)
         JS  = ISRSDL(1)+1
         IF (JS.GT.2) JS = 2
         JAZ = JA-JZ
*
         COMSCW = ZERZER
         IF (LESCO(JSCRNG).EQ.-1) THEN
            COMSCW = ONEONE
         ELSEIF (LESCO(JSCRNG).EQ.9999) THEN
            IF ((JA.EQ.IASCO(JSCRNG)).AND.(JZ.EQ.IZSCO(JSCRNG)).AND.
     &                                    (JS.EQ.ISSCO(JSCRNG))) THEN
               COMSCW = ONEONE
            ELSE
               COMSCW = ZERZER
            ENDIF
         ELSEIF ((LESCO(JSCRNG).EQ.1).OR.(LESCO(JSCRNG).EQ.10).OR.
     &           (LESCO(JSCRNG).EQ.100)) THEN
            FACT = 10.0D0/DBLE(LESCO(JSCRNG))
            IF (XLEO10(JZ,JAZ,JS).GT.ZERZER) 
     &         COMSCW = ONEONE/(FACT*XLEO10(JZ,JAZ,JS))
         ELSEIF (LESCO(JSCRNG).EQ.2) THEN
            IF (XLESWS(JZ,JAZ,JS).GT.ZERZER)
     &         COMSCW = ONEONE/XLESWS(JZ,JAZ,JS)
         ELSEIF (LESCO(JSCRNG).EQ.3) THEN
            IF (XLEECO(JZ,JAZ,JS).GT.ZERZER)
     &         COMSCW = ONEONE/XLEECO(JZ,JAZ,JS)
         ELSE
            WRITE(LUNOUT,*) ' COMSCW: invalid option ',LESCO(JSCRNG)
            STOP
         ENDIF

      ENDIF

      RETURN
*=== End of function Comscw ===========================================*
      END
