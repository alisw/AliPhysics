/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:35:18  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:18  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:39  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CIPERM
  CERN PROGLIB# Z311    CIPERM          .VERSION KERNFOR  4.34  930114
  ORIG. 03/06/92, JZ
      CALL CIPERM (NPERM)
      set permission mask NPERM to be used in next call to CIOPEN
*/
#if defined(CERNLIB_QX_SC)
void type_of_call ciperm_(nperm)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call ciperm(nperm)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CIPERM(nperm)
#endif
      int  *nperm;
{
      extern int ciopen_perm;

      ciopen_perm = *nperm & 0777;
      return;
}
/*> END <----------------------------------------------------------*/
