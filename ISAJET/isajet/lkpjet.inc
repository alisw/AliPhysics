#ifndef CERNLIB_ISAJET_LKPJET_INC
#define CERNLIB_ISAJET_LKPJET_INC
*
*
* lkpjet.inc
*
C----------------------------------------------------------------------
C-
C-   Name     LKPJET.INC
C-   Purpose  Temporary link area for PJET banks
C-   Created  5-DEC-1989  CHIP STEWART (HBP)
C-   Updated  13-JAN-1990   Harrison B. Prosper
C-   Updated   6-NOV-1990   Chip Stewart  - ADDED ISP1,ISV1
C-
C----------------------------------------------------------------------
C
C ****  JPJET(1)        User flag
C ****  JPJET(2)        System word
C ****  JPJET(3)        First link in area (= KPJET(1))
C
      INTEGER PJLON
      PARAMETER( PJLON = 1 )        ! Activate link area
C
      INTEGER PJLOFF
      PARAMETER( PJLOFF= 0 )        ! De-activate link area
C
      INTEGER MXPJET
      PARAMETER( MXPJET = 8 )
      INTEGER LPJHD,LPJET,LPJPT,LISAQ,LISAJ,LISP1,LISV1
      INTEGER JPJET,KPJET(MXPJET)
      EQUIVALENCE ( LPJHD, KPJET(1)  )
      EQUIVALENCE ( LPJET, KPJET(2)  )
      EQUIVALENCE ( LPJPT, KPJET(3)  )
      EQUIVALENCE ( LISAQ, KPJET(5)  )
      EQUIVALENCE ( LISAJ, KPJET(6)  )
      EQUIVALENCE ( LISP1, KPJET(7)  )
      EQUIVALENCE ( LISV1, KPJET(8)  )
      COMMON /LKPJET/ JPJET(2),KPJET
      SAVE /LKPJET/

#endif
