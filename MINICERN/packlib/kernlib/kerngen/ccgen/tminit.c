/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:48  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:46  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:28  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE TMINIT
  CERN PROGLIB#         TMINIT          .VERSION KERNFOR  4.36  930602
  ORIG. 20/07/90, RH + JZ
  Fortran interface routine to initialize TMPRO / TMREAD
      CALL TMINIT (INIT)
*/
#include <stdio.h>
#if defined(CERNLIB_QX_SC)
void type_of_call tminit_(ptinit)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call tminit(ptinit)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call TMINIT(ptinit)
#endif
      int  *ptinit;
{
      *ptinit = 7;
/*    setbuf (stdout,NULL);        */
      return;
}
/*> END <----------------------------------------------------------*/
