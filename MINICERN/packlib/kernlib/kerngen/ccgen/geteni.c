/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/09/02 14:26:36  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/04 17:34:17  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:27  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:22  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE GETENI
  CERN PROGLIB# Z265    GETENI          .VERSION KERNFOR  4.31  911111
  ORIG. 22/02/91, JZ
  Fortran interface routine to getenv

      CALL GETENVF (NAME, TEXT*)

          NAME  the name of the environment variable,
          TEXT  returns its value
                ISLATE(1) returns its length
*/
#include <stdio.h>
#include "kerngen/fortchar.h"
#if defined(CERNLIB_QX_SC)
void type_of_call geteni_(fname, ftext, lgtext, lgname)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call geteni(fname, ftext, lgtext, lgname)
#endif
#if defined(CERNLIB_QXCAPT)
#  ifdef CERNLIB_MSSTDCALL
    void type_of_call GETENI(fname, len_fname, ftext, len_ftext, lgtext, lgname)
     int len_fname, len_ftext;
#  else
    void type_of_call GETENI(fname, ftext, lgtext, lgname)
# endif
#endif
#if defined(CERNLIB_QMCRY)
      _fcd  fname,  ftext;
#endif
#if !defined(CERNLIB_QMCRY)
      char *fname, *ftext;
#endif
      int  *lgtext, *lgname;
{
      char *ptname, *fchtak();
      char *pttext, *getenv();
      int  fchput();

      pttext = NULL;
      ptname = fchtak(fname,*lgname);
      if (ptname == NULL)          goto out;
      pttext = getenv (ptname);
      free(ptname);

out:  *lgtext = fchput (pttext,ftext,*lgtext);
      return;
}
/*> END <----------------------------------------------------------*/
