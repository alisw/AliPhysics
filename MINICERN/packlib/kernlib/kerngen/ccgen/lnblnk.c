/*
 * $Id$
 *
 * $Log$
 * Revision 1.4  1997/09/02 14:26:37  mclareni
 * WINNT correction
 *
 * Revision 1.3  1997/02/04 17:34:31  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.2  1996/09/20 14:51:17  cernlib
 * Linux added
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:35  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:28  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMHPX)||defined(CERNLIB_QMOS9)
#include "hpxgs/lnblnk.c"
#elif defined(CERNLIB_QMIRT)||defined(CERNLIB_QMIRTD)
#include "irtgs/lnblnk.c"
#elif defined(CERNLIB_QMVAOS)||defined(CERNLIB_QMVMI)||defined(CERNLIB_LINUX)||defined(CERNLIB_MSSTDCALL)
#include "allgs/lnblnk.c"
#endif 
