/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:35:13  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:12  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:36  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
/*>    ROUTINE CFSIZE
  CERN PROGLIB# Z310    CFSIZE          .VERSION KERNFOR  4.29  910718
  ORIG. 12/01/91, JZ
      CALL CFSIZE (LUNDES, MEDIUM, NWREC, JRECL, ISTAT)
      get the position of the end-of-file and position to it :
       LUNDES  file descriptor
       MEDIUM  = 0,1,2,3 : primary disk/tape, secondary disk/tape
       NWREC   number of words per record
      *JRECL   number of records before end-of-file
      *ISTAT   status, =zero if success
*/
#include "kerngen/cf_seek.h"
#include "kerngen/cf_xaft.h"
#include "kerngen/wordsizc.h"
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QX_SC)
void type_of_call cfsize_(lundes, medium, nwrec, jrecl, stat)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call cfsize(lundes, medium, nwrec, jrecl, stat)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CFSIZE(lundes, medium, nwrec, jrecl, stat)
#endif
      int  *lundes, *medium, *nwrec, *jrecl, *stat;
{
      int  fildes;
      int  nboff;

/*        position the file to the end     */

      fildes  = *lundes;
      nboff = lseek (fildes, 0, 2);
      if (nboff < 0)               goto trouble;

/*        get position of the file        */

      nboff = nboff / NBYTPW;
      nboff = nboff / *nwrec;
      *jrecl = nboff;
      *stat = 0;
      return;

trouble:  *stat = -1;
          perror (" error in CFSIZE");
          return;
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_CFSIZE
#undef CERNLIB_TCGEN_CFSIZE
#endif
