/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:38  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:39  mclareni
 * All mods for Winnt 96a on winnt branch
 * 
 * Revision 1.1.1.1  1996/02/15 17:49:29  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#if defined(CERNLIB_QENVBSD)
#include "setenvbsd.c"
#else
#include "setenvsy5.c"
#endif
