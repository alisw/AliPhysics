*CMZ :          17/07/98  15.44.31  by  Federico Carminati
*-- Author :
C*********************************************************************
C*********************************************************************
C*                                                                  **
C*                                                      May 1990    **
C*                                                                  **
C*   The Lund Monte Carlo for Jet Fragmentation and e+e- Physics    **
C*                                                                  **
C*                        JETSET version 7.3                        **
C*                                                                  **
C*                        Torbjorn Sjostrand                        **
C*                                                                  **
C*                    CERN/TH, CH-1211 Geneva 23                    **
C*                BITNET/EARN address TORSJO@CERNVM                 **
C*                       Tel. +22 - 767 28 20                       **
C*                                                                  **
C*          LUSHOW is written together with Mats Bengtsson          **
C*                                                                  **
C*           A complete manual exists on a separate file            **
C*         Please report any program errors to the author!          **
C*                                                                  **
C*                   Copyright Torbjorn Sjostrand                   **
C*                                                                  **
C* 	 Modified by F. Antinori 5.12.91 for high multiplicity      **
C*********************************************************************
C*********************************************************************
C                                                                    *
C  List of subprograms in order of appearance, with main purpose     *
C  (S = subroutine, F = function, B = block data)                    *
C                                                                    *
C  S   LU1ENT   to fill one entry (= parton or particle)             *
C  S   LU2ENT   to fill two entries                                  *
C  S   LU3ENT   to fill three entries                                *
C  S   LU4ENT   to fill four entries                                 *
C  S   LUJOIN   to connect entries with colour flow information      *
C  S   LUGIVE   to fill (or query) commonblock variables             *
C  S   LUEXEC   to administrate fragmentation and decay chain        *
C  S   LUPREP   to rearrange showered partons along strings          *
C  S   LUSTRF   to do string fragmentation of jet system             *
C  S   LUINDF   to do independent fragmentation of one or many jets  *
C  S   LUDECY   to do the decay of a particle                        *
C  S   LUKFDI   to select parton and hadron flavours in fragm        *
C  S   LUPTDI   to select transverse momenta in fragm                *
C  S   LUZDIS   to select longitudinal scaling variable in fragm     *
C  S   LUSHOW   to do timelike parton shower evolution               *
C  S   LUBOEI   to include Bose-Einstein effects (crudely)           *
C  F   ULMASS   to give the mass of a particle or parton             *
C  S   LUNAME   to give the name of a particle or parton             *
C  F   LUCHGE   to give three times the electric charge              *
C  F   LUCOMP   to compress standard KF flavour code to internal KC  *
C  S   LUERRM   to write error messages and abort faulty run         *
C  F   ULALPS   to give the alpha_strong value                       *
C  F   ULANGL   to give the angle from known x and y components      *
C  F   RLU      to provide a random number generator                 *
C  S   RLUGET   to save the state of the random number generator     *
C  S   RLUSET   to set the state of the random number generator      *
C  S   LUROBO   to rotate and/or boost an event                      *
C  S   LUEDIT   to remove unwanted entries from record               *
C  S   LULIST   to list event record or particle data                *
C  S   LUUPDA   to update particle data                              *
C  F   KLU      to provide integer-valued event information          *
C  F   PLU      to provide real-valued event information             *
C  S   LUSPHE   to perform sphericity analysis                       *
C  S   LUTHRU   to perform thrust analysis                           *
C  S   LUCLUS   to perform three-dimensional cluster analysis        *
C  S   LUCELL   to perform cluster analysis in (eta, phi, E_T)       *
C  S   LUJMAS   to give high and low jet mass of event               *
C  S   LUFOWO   to give Fox-Wolfram moments                          *
C  S   LUTABU   to analyze events, with tabular output               *
C                                                                    *
C  S   LUEEVT   to administrate the generation of an e+e- event      *
C  S   LUXTOT   to give the total cross-section at given CM energy   *
C  S   LURADK   to generate initial state photon radiation           *
C  S   LUXKFL   to select flavour of primary qqbar pair              *
C  S   LUXJET   to select (matrix element) jet multiplicity          *
C  S   LUX3JT   to select kinematics of three-jet event              *
C  S   LUX4JT   to select kinematics of four-jet event               *
C  S   LUXDIF   to select angular orientation of event               *
C  S   LUONIA   to perform generation of onium decay to gluons       *
C                                                                    *
C  S   LUHEPC   to convert between /LUJETS/ and /HEPEVT/ records     *
C  S   LUTEST   to test the proper functioning of the package        *
C  B   LUDATA   to contain default values and particle data          *
C                                                                    *
C*********************************************************************

      SUBROUTINE LU1ENT(IP,KF,PE,THE,PHI)

C...Purpose: to store one parton/particle in commonblock LUJETS.
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEEP,LUDAT2.
      COMMON /LUDAT2/ KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT2/
*KEND.

C...Standard checks.
      MSTU(28)=0
      IF(MSTU(12).GE.1) CALL LULIST(0)
      IPA=MAX(1,IABS(IP))
      IF(IPA.GT.MSTU(4)) CALL LUERRM(21,
     &'(LU1ENT:) writing outside LUJETS memory')
      KC=LUCOMP(KF)
      IF(KC.EQ.0) CALL LUERRM(12,'(LU1ENT:) unknown flavour code')

C...Find mass. Reset K, P and V vectors.
      PM=0.
      IF(MSTU(10).EQ.1) PM=P(IPA,5)
      IF(MSTU(10).GE.2) PM=ULMASS(KF)
      DO 100 J=1,5
      K(IPA,J)=0
      P(IPA,J)=0.
  100 V(IPA,J)=0.

C...Store parton/particle in K and P vectors.
      K(IPA,1)=1
      IF(IP.LT.0) K(IPA,1)=2
      K(IPA,2)=KF
      P(IPA,5)=PM
      P(IPA,4)=MAX(PE,PM)
      PA=SQRT(P(IPA,4)**2-P(IPA,5)**2)
      P(IPA,1)=PA*SIN(THE)*COS(PHI)
      P(IPA,2)=PA*SIN(THE)*SIN(PHI)
      P(IPA,3)=PA*COS(THE)

C...Set N. Optionally fragment/decay.
      N=IPA
      IF(IP.EQ.0) CALL LUEXEC

      RETURN
      END
