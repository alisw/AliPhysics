/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/10/23 16:25:09  mclareni
 * NT mods, mostly C Fortran interface
 *
 * Revision 1.2  1997/02/04 17:34:22  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:30  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:23  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"

/*>    ROUTINE JMPLONG
  CERN PROGLIB#         JMPLONG         .VERSION KERNFOR  4.36  930602
  Fortran interface routine to longjmp for JMPSET
*/
#if defined(CERNLIB_QX_SC)
void type_of_call jmplong_(area,fnum)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call jmplong(area,fnum)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call JMPLONG(area,fnum)
#endif
      char *area;
      int  *fnum;
{
      int  num;

      num = *fnum;
#if defined(CERNLIB_QSIGJMP)
      siglongjmp(area,num);
#endif
#if !defined(CERNLIB_QSIGJMP)
      longjmp(area,num);
#endif
}
/*> END <----------------------------------------------------------*/
