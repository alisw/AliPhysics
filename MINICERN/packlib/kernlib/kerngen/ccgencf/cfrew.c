/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:35:12  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:11  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:36  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CFREW
  CERN PROGLIB# Z310    CFREW           .VERSION KERNFOR  4.29  910718
  ORIG. 12/01/91, JZ
      CALL CFREW (LUNDES,MEDIUM)
      rewind the file :
       LUNDES  file descriptor
       MEDIUM  = 0,1,2,3 : primary disk/tape, secondary disk/tape
*/
#include "kerngen/cf_seek.h"
#include "kerngen/cf_xaft.h"

#if defined(CERNLIB_QX_SC)
void type_of_call cfrew_(lundes, medium)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call cfrew(lundes, medium)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call CFREW(lundes, medium)
#endif
      int  *lundes, *medium;
{
      int  fildes;
      int  newpos;

      fildes = *lundes;
      newpos = lseek (fildes, 0, 0);
      if (newpos < 0)              goto trouble;
      return;

trouble:  perror (" error in CFREW");
          return;
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_CFREW
#undef CERNLIB_TCGEN_CFREW
#endif
