/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:43  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:43  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:27  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
/*>    ROUTINE SIGUNBL
  CERN PROGLIB#         SIGUNBL         .VERSION KERNFOR  4.42  951011
  ORIG. 10/10/95, JZ
  unblock all signals    */
#include <stdio.h>
#include <signal.h>
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QX_SC)
void type_of_call sigunbl_()
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call sigunbl()
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call SIGUNBL()
#endif
{
#ifndef CERNLIB_WINNT
      sigset_t   newmask;

      sigemptyset(&newmask);
      sigprocmask (SIG_SETMASK, &newmask, NULL);
#else
      DoAttention(" Attention !!! SIGUNBL is not implemented for Windows NT\n");
#endif
      return;
}
/*> END <----------------------------------------------------------*/
