/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:35:16  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:15  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:37  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CICLOS
  CERN PROGLIB# Z311    CICLOS          .VERSION KERNFOR  4.31  911111
  ORIG. 12/10/91, JZ
      CALL CICLOS (LUNDES)
      close the file :
       LUNDES  file descriptor
*/
#include "kerngen/cf_clos.h"
#include "kerngen/cf_xaft.h"

#if defined(CERNLIB_QX_SC)
void type_of_call ciclos_(lundes)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call ciclos(lundes)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CICLOS(lundes)
#endif
      int  *lundes;
{
      int  fildes;

      fildes = *lundes;
      close (fildes);
      return;
}
/*> END <----------------------------------------------------------*/
