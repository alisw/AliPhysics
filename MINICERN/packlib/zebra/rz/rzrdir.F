*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZRDIR(MAXDIR,CHDIR,NDIR)
*
************************************************************************
*
*       Returns the list of subdirectories of the CWD
* Input:
*   MAXDIR  Length of the character array CHDIR
* Output:
*   CHDIR*  Character  array  which  will  contain  the  directory  names
*           attached to the CWD.  If the  length of the directory name is
*           greater then the length of one  element of CHDIR (as obtained
*           by the LEN function),  only as many characters as will fit in
*           the array element are returned, and an error code will be set
*           in IQUEST(1).
*   NDIR*   Actual number of subdirectories attached to the CWD
*           If this number is greater than MAXDIR,  only the first MAXDIR
*           directory names will be returned in CHDIR (see IQUEST(12))
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 11.04.86
*  Last mod: 19.08.86
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzk.inc"
      DIMENSION IHDIR(4)
      CHARACTER*(*) CHDIR(*)
*
*-----------------------------------------------------------------------
*
      IQUEST(1)=0
      NDIR=0
      IF(LQRS.EQ.0) GO TO 99
      IF(LCDIR.EQ.0)GO TO 99
      NCHMAX= LEN(CHDIR(1))
      NCH   = MIN(NCHMAX,16)
      LS    = IQ(KQSP+LCDIR+KLS)
      NSDIR = IQ(KQSP+LCDIR+KNSD)
      IF(NSDIR.LE.MAXDIR)THEN
         NDIR=NSDIR
      ELSE
         NDIR=MAXDIR
         IQUEST(1)=1
      ENDIF
      IQUEST(11)=NSDIR
*
      DO 10 I=1,NDIR
         CALL ZITOH(IQ(KQSP+LCDIR+LS+7*(I-1)),IHDIR,4)
         CALL UHTOC(IHDIR,4,CHDIR(I),NCH)
  10  CONTINUE
*
  99  RETURN
      END
