/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  1996/12/06 08:49:18  gunter
 * Correct problem with ucopy2 on OSF (f77 - version 4) by replcing ucopy2
 * with a c routine ucopy2c.c calling memmove. This is also used on some
 * other unix - where it was faster.
 *
 * Kernlib
 *
 *  C version of ucopy2 using memmove; memmove is POSIX
 */
#include "kerngen/pilot.h"

#include <string.h>

#if defined(CERNLIB_QXNO_SC)
#define ucopy2_ ucopy
#endif

void ucopy2_(int *from, int *to, int *nwords)
{
	if ( *nwords > 1 ) {
	   memmove( (void *)to, (void *)from,(size_t) 4*(*nwords));
	} else {
	   if ( *nwords == 1 ) {
	      *to = *from;
	   }
        /* else  nothing to do */
	}
}
