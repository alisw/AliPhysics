/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:35:11  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:10  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:36  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CFPERM
  CERN PROGLIB# Z311    CFPERM          .VERSION KERNFOR  4.34  930114
  ORIG. 03/06/92, JZ
      CALL CFPERM (NPERM)
      set permission mask NPERM to be used in next call to CFOPEN
*/
#if defined(CERNLIB_QX_SC)
void type_of_call cfperm_(nperm)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call cfperm(nperm)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CFPERM(nperm)
#endif
      int  *nperm;
{
      extern int cfopen_perm;

      cfopen_perm = *nperm & 0777;
      return;
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_CFPERM
#undef CERNLIB_TCGEN_CFPERM
#endif
