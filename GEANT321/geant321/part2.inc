*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_PART2_INC
#define CERNLIB_GEANT321_PART2_INC
*
*
* part2.inc
*
*
*=== part2 ============================================================*
*
*----------------------------------------------------------------------*
*     Include file: part2 copy        Revised on 20-7-90 by A. Ferrari *
*     Note: see also part copy and part3 copy                          *
*     Changes: none                                                    *
*     Included in the following subroutines or functions: not updated  *
*                                                                      *
*     Description of the common block(s) and variable(s)               *
*                                                                      *
*         Kptoip = conversion from part to paprop numbering            *
*         Iptokp = conversion from paprop to part numbering            *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER (IDMAX8=183)
      CHARACTER*8  ANAME
      COMMON / FKPART / AAM  (IDMAX8), GA     (IDMAX8), TAU    (IDMAX8),
     &                AAMDSC (IDMAX8), IICH   (IDMAX8), IIBAR  (IDMAX8),
     &                K1     (IDMAX8), K2     (IDMAX8), KPTOIP (IDMAX8),
     &                IPTOKP (39)
      COMMON / FKCHPA / ANAME (IDMAX8)
 

#endif
