/* 
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:42  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:42  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:27  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
/*>    ROUTINE SIGPRNT
  CERN PROGLIB#         SIGPRNT         .VERSION KERNFOR  4.42  951011
  ORIG. 10/10/95, JZ
  print the mask of blocked signals    */
#include <stdio.h>
#include <signal.h>
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QX_SC)
void type_of_call sigprnt_()
#endif
#if defined(CERNLIB_QXNO_SC)
void  type_of_call sigprnt()
#endif
#if defined(CERNLIB_QXCAPT)
void  type_of_call SIGPRNT()
#endif
{
#ifndef CERNLIB_WINNT
      sigset_t   oldmask;

      sigprocmask (NULL, NULL, &oldmask);

      printf (" blocked signals are: %x\n", oldmask);
#else
      printf ("Printing the mask of blocked signal isn't implemented for Windows\n");
#endif
      return;
}
/*> END <----------------------------------------------------------*/
