*/////////////////////////////////////////////////////////////////////////////////////
*//                                                                                 //
*//  !!!!!!! WARNING!!!!!   This source may be  agressive !!!!                      //
*//                                                                                 //
*//  Due to short common block names it may owerwrite variables in other parts      //
*//  of the code.                                                                   //
*//                                                                                 //
*//  One should add suffix c_Photos_ to names of all commons as soon as possible!!  //
*//                                                                                 //
*/////////////////////////////////////////////////////////////////////////////////////

C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOtos CDE's
C.
C.    Purpose:  Keep definitions  for PHOTOS QED correction Monte Carlo.
C.
C.    Input Parameters:   None
C.
C.    Output Parameters:  None
C.
C.    Author(s):  Z. Was, B. van Eijk             Created at:  29/11/89
C.                                                Last Update: 10/08/93
C.
C. =========================================================
C.    General Structure Information:                       =
C. =========================================================
C:   ROUTINES:
C.             1) INITIALIZATION:
C.                                      PHOCDE
C.                                      PHOINI
C.                                      PHOCIN
C.                                      PHOINF
C.             2) GENERAL INTERFACE:
C.                                      PHOTOS
C.                                      PHOTOS_GET
C.                                      PHOTOS_SET
C.                                      PHOTOS_MAKE
C.                                      PHOBOS
C.                                      PHOIN
C.                                      PHOTWO (specific interface
C.                                      PHOOUT
C.                                      PHOCHK
C.                                      PHTYPE (specific interface
C.                                      PHOMAK (specific interface
C.             3) QED PHOTON GENERATION:
C.                                      PHINT
C.                                      PHOPRE
C.                                      PHOOMA
C.                                      PHOENE
C.                                      PHOCOR
C.                                      PHOFAC
C.                                      PHODO
C.             4) UTILITIES:
C.                                      PHOTRI
C.                                      PHOAN1
C.                                      PHOAN2
C.                                      PHOBO3
C.                                      PHORO2
C.                                      PHORO3
C.                                      PHORIN
C.                                      PHORAN
C.                                      PHOCHA
C.                                      PHOSPI
C.                                      PHOERR
C.                                      PHOREP
C.                                      PHLUPA
C.                                      PHCORK
C.   COMMONS:
C.   NAME     USED IN SECT. # OF OCC.     Comment
C.   PHOQED   1) 2)            3      Flags whether emisson to be gen. 
C.   PHOLUN   1) 4)            6      Output device number
C.   PHOCOP   1) 3)            4      photon coupling & min energy
C.   PHPICO   1) 3) 4)         5      PI & 2*PI
C.   PHSEED   1) 4)            3      RN seed 
C.   PHOSTA   1) 4)            3      Status information
C.   PHOKEY   1) 2) 3)         7      Keys for nonstandard application
C.   PHOVER   1)               1      Version info for outside
C.   HEPEVT   2)               2      PDG common
C.   PH_HEPEVT2)               8      PDG common internal
C.   PHOEVT   2) 3)           10      PDG branch
C.   PHOIF    2) 3)            2      emission flags for PDG branch 
C.   PHOMOM   3)               5      param of char-neutr system
C.   PHOPHS   3)               5      photon momentum parameters
C.   PHOPRO   3)               4      var. for photon rep. (in branch)
C.   PHOCMS   2)               3      parameters of boost to branch CMS
C.   PHNUM    4)               1      event number from outside         
C.----------------------------------------------------------------------
      SUBROUTINE PHOINI
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays INItialisation
C.
C.    Purpose:  Initialisation  routine  for  the  PHOTOS  QED radiation
C.              package.  Should be called  at least once  before a call
C.              to the steering program 'PHOTOS' is made.
C.
C.    Input Parameters:   None
C.
C.    Output Parameters:  None
C.
C.    Author(s):  Z. Was, B. van Eijk             Created at:  26/11/89
C.                                                Last Update: 12/04/90
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER INIT
      SAVE INIT
      DATA INIT/ 0/
C--
C--   Return if already initialized...
      IF (INIT.NE.0) RETURN
      INIT=1
C--
C--   Preset parameters in PHOTOS commons
      CALL PHOCIN
C--
C--   Print info
      CALL PHOINF
C--
C--   Initialization of kinematic correction against rounding errors.
C--   Default parameter is 1 (no correction) optionally 2, 3, 4
      CALL PHCORK(1)
C--
C--   Initialize Marsaglia and Zaman random number generator
      CALL PHORIN
      RETURN
      END
