/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:35:22  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:21  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:39  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CITELL
  CERN PROGLIB# Z311    CITELL          .VERSION KERNFOR  4.31  911111
  ORIG. 12/10/91, JZ
      CALL CITELL (LUNDES, JCBYT, ISTAT)
      get the current position of the file :
       LUNDES  file descriptor
      *JCBYT   number of bytes before current
      *ISTAT   status, =zero if success
*/
#include "kerngen/cf_seek.h"
#include "kerngen/cf_xaft.h"
#if defined(CERNLIB_QX_SC)
void type_of_call citell_(lundes, jcbyt, stat)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call citell(lundes, jcbyt, stat)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CITELL(lundes, jcbyt, stat)
#endif
      int  *lundes, *jcbyt, *stat;
{
      int  fildes;
      int  nboff;

/*        get position of the file        */

      fildes = *lundes;
      nboff  = lseek (fildes, 0, 1);
      if (nboff < 0)               goto trouble;
      *jcbyt = nboff;
      *stat = 0;
      return;

trouble:  *stat = -1;
          perror (" error in CITELL");
          return;
}
/*> END <----------------------------------------------------------*/
