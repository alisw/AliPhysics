/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/09/02 14:26:36  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/04 17:34:20  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:29  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:23  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE GETWDI
  CERN PROGLIB# Z265    GETWDI          .VERSION KERNFOR  4.38  931108
  ORIG. 22/02/91, JZ
  Fortran interface routine to getwd

      CALL GETWDF (TEXT*)

      returns the name of the c.w.d. in TEXT
      ISLATE(1) returns its lenth NTEXT
*/
#include <stdio.h>
#ifdef WIN32
#include <direct.h>
# ifndef getcwd
#   define getcwd _getcwd
# endif
#endif
#include "kerngen/fortchar.h"
#if defined(CERNLIB_QX_SC)
void type_of_call getwdi_(fname, lgname)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call getwdi(fname, lgname)
#endif
#if defined(CERNLIB_QXCAPT)
# ifdef CERNLIB_MSSTDCALL
    void type_of_call GETWDI(fname, len_fname, lgname)
    int len_fname;
#  else
    void type_of_call GETWDI(fname, lgname)
#  endif
#endif
#if defined(CERNLIB_QMCRY)
      _fcd fname;
#endif
#if !defined(CERNLIB_QMCRY)
      char *fname;
#endif
      int  *lgname;
{
      char *malloc();
      char *ptalc, *pttext;
      int  fchput();
      int  nalc;
#if !defined(CERNLIB_QGETCWD)
      char *getwd();
#endif
#if defined(CERNLIB_QGETCWD)
      char *getcwd();
      int  nsize;
#endif

      pttext = NULL;
      nalc   = 2048;
      ptalc  = malloc(nalc);
      if (ptalc == NULL)           goto out;

#if !defined(CERNLIB_QGETCWD)
      pttext = getwd (ptalc);
#endif
#if defined(CERNLIB_QGETCWD)
      nsize  = nalc;
      pttext = getcwd (ptalc, nsize);
#endif

out:  *lgname = fchput (pttext,fname,*lgname);
      if (ptalc != NULL)   free(ptalc);
      return;
}
/*> END <----------------------------------------------------------*/
