*
** 21-feb-2001 get from polar.ihep.su@DATA$STAR:[PAVLINOV.RSC]SGPDGE.CAR
*
*-- Author :    David Adams   03/11/90
      SUBROUTINE SGPDGE( I, PDGGEA )
CA)========================================(begin of short description)
C Converts a PDG particle number to a GEANT particle number.
CB)-----------------------------------------(begin of long description)
C
C    Purpose and methods :
C Translates the PDG code I to the GEANT particle code PDGGEA.
C Zero is returned if the PDG code cannot be translated.
C
C    Called by: SGGENE
C
C    Input arguments : I - PDG particle code
C
C    Input/output arguments : none
C
C    Output arguments: PDGGEA - GEANT particle code
C
C    Side effects (other variables/banks modified, etc.)
C
C    Created:   dd-mmm-1990   David Adams
C    Modified:  19-FEB-1991   David Adams - fixed bugs in particle codes.
C
CC)-------------------------------------------(end of long description)
CD)----------------------------------(declare common blocks, variables)
C  Common blocks
C  Variables in argument list
C  Save variables
C  Local variables
      INTEGER I, IABS, PDGGEA
C  External references
C  Statement functions
C  Local data statements
CE)===============================================(end of declarations)
C
      PDGGEA = 0
      IABS = ABS(I)

C Leptons
      IF ( IABS .LT. 20 ) THEN

C . e-, e+, nue
       IF     ( I .EQ.  11 ) THEN
        PDGGEA = 3
       ELSEIF ( I .EQ. -11 ) THEN
        PDGGEA = 2
       ELSEIF ( IABS .EQ.  12 ) THEN
        PDGGEA = 4

C . mu-, mu+, numu
       ELSEIF ( I .EQ.  13 ) THEN
        PDGGEA = 6
       ELSEIF ( I .EQ. -13 ) THEN
        PDGGEA = 5
       ELSEIF ( IABS .EQ.  14 ) THEN
        PDGGEA = 4

C . tau-, tau+, nutau
       ELSEIF ( I .EQ.  15 ) THEN
        PDGGEA = 34
       ELSEIF ( I .EQ. -15 ) THEN
        PDGGEA = 33
       ELSEIF ( IABS .EQ.  16 ) THEN
        PDGGEA = 4

       ENDIF

       RETURN

C Bosons.
      ELSEIF ( IABS .LT. 100 ) THEN

C . gamma, W+, W-, Z
       IF     ( I .EQ.  22 ) THEN
        PDGGEA = 1
       ELSEIF ( I .EQ.  24 ) THEN
        PDGGEA = 42
       ELSEIF ( I .EQ. -24 ) THEN
        PDGGEA = 43
       ELSEIF ( I .EQ.  23 ) THEN
        PDGGEA = 44
       ENDIF

       RETURN


C Mesons
      ELSEIF ( IABS .LT. 1000 ) THEN

C . pi0, pi+, pi-
       IF     ( I .EQ.  111 ) THEN
        PDGGEA = 7
       ELSEIF ( I .EQ.  211 ) THEN
        PDGGEA = 8
       ELSEIF ( I .EQ. -211 ) THEN
        PDGGEA = 9

C . K0L, K0S, K+, K-
       ELSEIF ( I .EQ.  130 ) THEN
        PDGGEA = 10
       ELSEIF ( I .EQ.  310 ) THEN
        PDGGEA = 16
       ELSEIF ( I .EQ.  321 ) THEN
        PDGGEA = 11
       ELSEIF ( I .EQ. -321 ) THEN
        PDGGEA = 12

C . K0L, K0S, K+, K-
       ELSEIF ( I .EQ.  130 ) THEN
        PDGGEA = 10
       ELSEIF ( I .EQ.  310 ) THEN
        PDGGEA = 16
       ELSEIF ( I .EQ.  321 ) THEN
        PDGGEA = 11
       ELSEIF ( I .EQ. -321 ) THEN
        PDGGEA = 12

C . eta
       ELSEIF ( I .EQ.  221 ) THEN
        PDGGEA = 17

C . D+, D-, D0, D0bar
       ELSEIF ( I .EQ.  411 ) THEN
        PDGGEA = 35
       ELSEIF ( I .EQ. -411 ) THEN
        PDGGEA = 36
       ELSEIF ( I .EQ.  421 ) THEN
        PDGGEA = 37
       ELSEIF ( I .EQ. -421 ) THEN
        PDGGEA = 38

C . DS+, DS- (previously F+, F-)
       ELSEIF ( I .EQ.  431 ) THEN
        PDGGEA = 39
       ELSEIF ( I .EQ. -431 ) THEN
        PDGGEA = 40

       ENDIF

      RETURN


C Hadrons
      ELSEIF ( IABS .LT. 10000 ) THEN

C . n, nbar, p, pbar
       IF     ( I .EQ.  2112 ) THEN
        PDGGEA = 13
       ELSEIF ( I .EQ. -2112 ) THEN
        PDGGEA = 25
       ELSEIF ( I .EQ.  2212 ) THEN
        PDGGEA = 14
       ELSEIF ( I .EQ. -2212 ) THEN
        PDGGEA = 15

C . lambda, antilambda
       ELSEIF ( I .EQ.  3122 ) THEN
        PDGGEA = 18
       ELSEIF ( I .EQ. -3122 ) THEN
        PDGGEA = 26

C . sigma +, 0, -, antisigma -, 0, +
       ELSEIF ( I .EQ.  3222 ) THEN
        PDGGEA = 19
       ELSEIF ( I .EQ.  3212 ) THEN
        PDGGEA = 20
       ELSEIF ( I .EQ.  3112 ) THEN
        PDGGEA = 21
       ELSEIF ( I .EQ. -3222 ) THEN
        PDGGEA = 27
       ELSEIF ( I .EQ. -3212 ) THEN
        PDGGEA = 28
       ELSEIF ( I .EQ. -3112 ) THEN
        PDGGEA = 29

C . xi 0, xi -, antixi 0, antixi +
       ELSEIF ( I .EQ.  3322 ) THEN
        PDGGEA = 22
       ELSEIF ( I .EQ.  3312 ) THEN
        PDGGEA = 23
       ELSEIF ( I .EQ. -3322 ) THEN
        PDGGEA = 30
       ELSEIF ( I .EQ. -3312 ) THEN
        PDGGEA = 31

C . omega, antiomega
       ELSEIF ( I .EQ.  3334 ) THEN
        PDGGEA = 24
       ELSEIF ( I .EQ. -3334 ) THEN
        PDGGEA = 32

C . lambda C +
       ELSEIF ( I .EQ. 4122 ) THEN
        PDGGEA = 41

       ENDIF

       RETURN


      ENDIF

      RETURN
      END

