/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:12  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:22  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:20  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE ABEND
  CERN PROGLIB# Z035    ABEND           .VERSION KERNFOR  4.31  911111
*/
#if defined(CERNLIB_QX_SC)
void type_of_call abend_()
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call abend()
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call ABEND()
#endif
{
    exit(7);
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_ABEND
#undef CERNLIB_TCGEN_ABEND
#endif
