/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:19  mclareni
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
#include "os9gs/getpidf.c"
#else
/*>    ROUTINE GETPIDF (IPID)
  CERN PROGLIB# Z265    GETPIDF         .VERSION KERNFOR  4.38  931108
  ORIG. 22/02/91, JZ
  Fortran interface routine to getpid
*/
#ifdef WIN32
#include <process.h>
#endif
#if defined(CERNLIB_QX_SC)
void type_of_call getpidf_(pid)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call getpidf(pid)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call GETPIDF(pid)
#endif
      int  *pid;
{
      int getpid();
      *pid = getpid();
      return;
}
/*> END <----------------------------------------------------------*/
#endif
