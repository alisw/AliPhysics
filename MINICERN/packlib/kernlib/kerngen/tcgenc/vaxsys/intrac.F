*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:40  mclareni
* Kernlib
*
*
      LOGICAL FUNCTION INTRAC(DUMMY)
C
C CERN PROGLIB# Z044    INTRAC          .VERSION KERNVAX  2.42  930609
C ORIG. 08/06/93  FCA + JZ
C
      COMMON /SLATE/ NDSLAT,NESLAT,SLDUMMY(38)

      INTEGER      DVI$_DEVCHAR, DEV$V_TRM
      PARAMETER   (DVI$_DEVCHAR  = '002'X)
      PARAMETER   (DEV$V_TRM     = '002'X)

      INTEGER*4    MVALUE
      INTEGER*4    LISTF(4)
      INTEGER*2    LISTH(8)
      EQUIVALENCE (LISTH(1),LISTF(1))
      EQUIVALENCE (MVALUE, NDSLAT)
      INTEGER      SYS$GETDVIW

C--    System service item list

      DATA LISTH(1),LISTH(2) / 4 , DVI$_DEVCHAR /
      DATA LISTF(3),LISTF(4) / 0 , 0 /


      LISTF(2) = %LOC(MVALUE)
      ISTAT    = SYS$GETDVIW(,,'SYS$INPUT',LISTF,,,,)

      MASK = MVALUE .AND. DEV$V_TRM
      IF (MASK.NE.0)               GO TO 91
      INTRAC = .FALSE.
      RETURN

   91 INTRAC = .TRUE.
      RETURN
      END
