/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:39  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:40  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:26  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE SETENI
  CERN PROGLIB# Z265    SETENI          .VERSION KERNFOR  4.36  930602
  ORIG. 22/02/91, JZ
  Fortran interface routine to setenv

      CALL SETENVF (NAME, TEXT*)

          NAME  the name of the environment variable,
          TEXT  the value to be assigned
*/
#include <stdio.h>
#include "kerngen/fortchar.h"
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QX_SC)
int type_of_call seteni_(fname, ftext, lgname, lgtext)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call seteni(fname, ftext, lgname, lgtext)
#endif
#if defined(CERNLIB_QXCAPT)
int type_of_call SETENI(fname, ftext, lgname, lgtext)
#endif
#if defined(CERNLIB_QMCRY)
      _fcd  fname,  ftext;
#endif
#if !defined(CERNLIB_QMCRY)
      char *fname, *ftext;
#endif
      int  *lgtext, *lgname;
{
      char *ptname, *pttext, *fchtak();
      int  nname, ntext, istat, setenv();

      istat = -1;
      nname = *lgname;
      ntext = *lgtext;

      ptname = fchtak(fname,nname);
      if (ptname == NULL)          goto out1;
      pttext = fchtak(ftext,ntext);
      if (pttext == NULL)          goto out2;
      istat = setenv (ptname, pttext);
      free(pttext);
out2: free(ptname);
out1: return istat;
}
/*> END <----------------------------------------------------------*/
