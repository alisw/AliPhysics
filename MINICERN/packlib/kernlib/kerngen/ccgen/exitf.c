/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:16  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:27  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:21  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE EXITF
  CERN PROGLIB# Z035    EXITF           .VERSION KERNFOR  4.39  940228
*/
#if defined(CERNLIB_QX_SC)
void type_of_call exitf_(st)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call exitf(st)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call EXITF(st)
#endif
      int  *st;
{
      exit(*st);
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_EXITF
#undef CERNLIB_TCGEN_EXITF
#endif
