/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:35:21  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:20  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:39  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CISIZE
  CERN PROGLIB# Z311    CISIZE          .VERSION KERNFOR  4.31  911111
  ORIG. 12/10/91, JZ
      CALL CISIZE (LUNDES, JBYTL, ISTAT)
      get the position of the end-of-file and position to it :
       LUNDES  file descriptor
      *JBYTL   number of bytes before end-of-file
      *ISTAT   status, =zero if success
*/
#include "kerngen/cf_seek.h"
#include "kerngen/cf_xaft.h"
#if defined(CERNLIB_QX_SC)
void type_of_call cisize_(lundes, jbytl, stat)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call cisize(lundes, jbytl, stat)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CISIZE(lundes, jbytl, stat)
#endif
      int  *lundes, *jbytl, *stat;
{
      int  fildes;
      int  nboff;

/*        position the file to the end     */

      fildes  = *lundes;
      nboff = lseek (fildes, 0, 2);
      if (nboff < 0)               goto trouble;

/*        get position of the file        */

      *jbytl = nboff;
      *stat = 0;
      return;

trouble:  *stat = -1;
          perror (" error in CISIZE");
          return;
}
/*> END <----------------------------------------------------------*/
