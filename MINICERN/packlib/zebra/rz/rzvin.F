*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:27:20  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:27  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZVIN(V,NDIM,N,KEYU,ICYCLE,CHOPT)
*
************************************************************************
*
*           Read Fortran array V from a RZ file
* Input:
*   NDIM*   Declared dimension of array V
*  *N*      With option O: N contains offset of first word to be read
*   KEYU    Keyword vector of the information to be read
*   ICYCLE  Cycle number of the key to be read
*           ICYCLE > highest cycle number means read the highest cycle
*           ICYCLE = 0 means read the lowest cycle
*   CHOPT   Character variable specifying the options selected.
*           data structure
*             default
*                   Same as 'D' below
*             'A'   Read continuation of the previously read data structure
*                   with identifier KEYU,ICYCLE
*                   Given that option implies that the record was written with
*                   the same option by a call to RZOUT.
*             'C'   Provide   information   about   the   cycle   numbers
*                   associated with KEY.
*                   The  total number  of  cycles  and the  cycle  number
*                   identifiers of the 19 highest  cycles are returned in
*                   IQUEST(50) and IQUEST(51..89) respectively
*             'D'   Read the  Data structure  with the  (key,cycle)  pair
*                   specified.
*             'N'   Read the neighbouring. keys (i.e. those preceding and
*                   following KEY).
*                   The  key-vectors of  the previous  and  next key  are
*                   available   respectively   as   IQUEST(31..35)    and
*                   IQUEST(41..45), see below.
*             'O'   Argument N contains offset of first word of vector
*                   to be read
* Output:
*   V       Array containing the contents of the record
*   N       Actual length of array V
*           If the current length is greater than the declared dimension,then
*           only NDIM are returned in V.
*           In that case,IQUEST(11) contains the number of words in the record
*
*  IQUEST(13)=IFORM
*
*  Integer value CHOPT for RZVOUT Meaning
*  ------------- ---------------- -------
*  3                              Array V contains only floating points
*  2          'I'                 Array V contains integers
*  1          'B'                 Array V contains bit-patterns
*  5          'H'                 Array V contains Holleriths
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 03.04.86
*  Last mod: 22.06.92 JDS - return IFORM in IQUEST(13)
*          : 04.11.94 JDS - option O - N contains offset of first
*                           word of vector to be read
*
************************************************************************
#include "zebra/rzcl.inc"
      CHARACTER*(*) CHOPT
      DIMENSION KEYU(*),V(*)
      EQUIVALENCE (IOPTA,IQUEST(91)), (IOPTC,IQUEST(92))
     +,           (IOPTD,IQUEST(93)), (IOPTN,IQUEST(94))
     +,           (IOPTO,IQUEST(95))
*
*-----------------------------------------------------------------------
*

#include "zebra/q_jbyt.inc"

*           Find KEY,CYCLE
*
      CALL RZINK(KEYU,ICYCLE,CHOPT)
      IF(IQUEST(1).NE.0)GO TO 99
*
*!! options set by call to UOPTC in RZINK !!
*
      IF(IOPTC.NE.0.AND.IOPTD.EQ.0)GO TO 99
      IF(IOPTN.NE.0.AND.IOPTD.EQ.0)GO TO 99
      ISTART = 1
      CALL UOPTC(CHOPT,'O',IOPTO)
      IF(IOPTO.NE.0) ISTART = N
*
      NW = IQUEST(12)
      N  = NW
      IF(NW.GT.NDIM)THEN
         IQUEST(11)=NW
         NW=NDIM
      ENDIF
*
      IFORM=JBYT(IQUEST(14),1,3)
      IF(IFORM.LE.0)THEN
         IQUEST(1)=2
         GO TO 99
      ENDIF
*
*           Read record into array V
*
      CALL RZREAD(V,NW,ISTART,IFORM)
*
*     Store I/O descriptor in IQUEST(13)
*
      IQUEST(13) = IFORM
*
  99  RETURN
      END
