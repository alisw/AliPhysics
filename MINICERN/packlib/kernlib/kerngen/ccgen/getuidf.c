/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:20  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:28  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:22  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QMOS9)
#include "os9gs/getuidf.c"
#else
/*>    ROUTINE GETUIDF
  CERN PROGLIB# Z265    GETUIDF         .VERSION KERNFOR  4.38  931108
  ORIG. 01/04/93, JS
  Fortran interface routine to getuid
*/
#include <sys/types.h>
#if defined(CERNLIB_QX_SC)
void type_of_call getuidf_(uid)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call getuidf(uid)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call GETUIDF(uid)
#endif
      int  *uid;
{
#ifndef WIN32
      uid_t  getuid();

      *uid = getuid();
#endif
      return;
}
/*> END <----------------------------------------------------------*/
#endif
