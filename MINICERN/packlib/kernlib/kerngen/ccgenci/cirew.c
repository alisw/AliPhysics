/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/09/16 09:45:00  mclareni
 * Typing error affecting VMS
 *
 * Revision 1.2  1997/02/04 17:35:20  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:19  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:39  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CIREW
  CERN PROGLIB# Z311    CIREW           .VERSION KERNFOR  4.31  911111
  ORIG. 12/10/91, JZ
      CALL CIREW (LUNDES)
      rewind the file :
       LUNDES  file descriptor
*/
#include "kerngen/cf_seek.h"
#include "kerngen/cf_xaft.h"

#if defined(CERNLIB_QX_SC)
void type_of_call cirew_(lundes)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call  cirew(lundes)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CIREW(lundes)
#endif
      int  *lundes;
{
      int  fildes;
      int  newpos;

      fildes = *lundes;
      newpos = lseek (fildes, 0, 0);
      if (newpos < 0)              goto trouble;
      return;

trouble:  perror (" error in CIREW");
          return;
}
/*> END <----------------------------------------------------------*/
