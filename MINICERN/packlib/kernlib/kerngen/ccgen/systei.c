/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/09/02 14:26:40  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/04 17:34:45  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:44  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:27  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE SYSTEI
  CERN PROGLIB# Z265    SYSTEI          .VERSION KERNFOR  4.31  911111
  ORIG. 22/02/91, JZ
  Fortran interface routine to system

      ISTAT =  SYSTEMF (TEXT)

          TEXT  the command to be executed                  .
         ISTAT  returns zero if successful
*/
#include <stdio.h>
#include "kerngen/fortchar.h"

#if defined(CERNLIB_QX_SC)
int type_of_call systei_(ftext,nsize)
#endif

#if defined(CERNLIB_QXNO_SC)
int type_of_call systei(ftext,nsize)
#endif

#if defined(CERNLIB_QXCAPT)
# ifndef CERNLIB_MSSTDCALL
   int type_of_call SYSTEI(ftext,nsize)
# else
   int type_of_call SYSTEI(ftext,len_ftext,nsize)
# endif
#endif

#if defined(CERNLIB_QMCRY)
      _fcd ftext;
#endif
#if !defined(CERNLIB_QMCRY)
      char *ftext;
#endif
      int  *nsize;
#ifdef CERNLIB_MSSTDCALL
      int len_ftext; 
#endif
{
      char *ptname, *fchtak();
      int  system();
      int  istat;

/*        get memory and copy TEXT terminated  */

      ptname = fchtak(ftext,*nsize);
      if (ptname == NULL)           goto bad;

      istat = system (ptname);
      free (ptname);
      return istat;

bad:  return -1;
}
/*> END <----------------------------------------------------------*/
