/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:18  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:27  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:22  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QMOS9)
#include "os9gs/getgidf.c"
#else
/*>    ROUTINE GETGIDF
  CERN PROGLIB# Z265    GETGIDF         .VERSION KERNFOR  4.38  931108
  ORIG. 01/04/93, JS
  Fortran interface routine to getgid
*/
#include <sys/types.h>
#if defined(CERNLIB_QX_SC)
void type_of_call getgidf_(gid)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call getgidf(gid)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call GETGIDF(gid)
#endif
      int *gid;
{
#ifndef WIN32
      gid_t  getgid();

      *gid = getgid();
#endif
      return;
}
/*> END <----------------------------------------------------------*/
#endif
