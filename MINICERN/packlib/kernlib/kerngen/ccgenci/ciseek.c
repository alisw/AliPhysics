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

/*>    ROUTINE CISEEK
  CERN PROGLIB# Z311    CISEEK          .VERSION KERNFOR  4.31  911111
  ORIG. 12/10/91, JZ
      CALL CISEEK (LUNDES, JCBYT, ISTAT)
      reposition the file :
       LUNDES  file descriptor
       JCBYT   number of bytes before current
      *ISTAT   status, =zero if success
*/
#include "kerngen/cf_seek.h"
#include "kerngen/cf_xaft.h"
#if defined(CERNLIB_QX_SC)
void type_of_call ciseek_(lundes, jcbyt, stat)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call ciseek(lundes, jcbyt, stat)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CISEEK(lundes, jcbyt, stat)
#endif
      int  *lundes, *jcbyt, *stat;
{
      int  fildes;
      int  nbdo;
      int  isw;

/*        position the file        */

      fildes = *lundes;
      nbdo   = *jcbyt;
      isw    = lseek (fildes, nbdo, 0);
      if (isw <  0)                goto trouble;
      *stat = 0;
      return;

trouble:  *stat = -1;
          perror (" error in CISEEK");
          return;
}
/*> END <----------------------------------------------------------*/
