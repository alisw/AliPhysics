/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:40  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:40  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:26  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
/*>    ROUTINE SETJMPF (AREA, ROUT)

  CERN PROGLIB#         SETJMPF         .VERSION KERNFOR  4.26  910313

       The function  setjmp  cannot be implemented
            by a Fortran interface routine

       Instead we provide now the function JMPSET;
       the present SETJMPF is obsolete and will be
       removed in the next update.
*/
#include <setjmp.h>
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QX_SC)
void type_of_call setjmpf_(area,ufun)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call setjmpf(area,ufun)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call SETJMPF(area,ufun)
#endif
   jmp_buf area;
#if defined(CERNLIB_QCCINDAD)
   void (** ufun)();
#endif
#if !defined(CERNLIB_QCCINDAD)
   void (* ufun)();
#endif
{
   void (* unext)();

#if defined(CERNLIB_QCCINDAD)
   unext = *ufun;
#endif
#if !defined(CERNLIB_QCCINDAD)
   unext = ufun;
#endif

   if ( setjmp(area) )       return;

   (* unext)();
}
/*> END <----------------------------------------------------------*/
