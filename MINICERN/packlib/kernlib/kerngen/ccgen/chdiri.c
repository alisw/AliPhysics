/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/09/02 14:26:34  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/04 17:34:14  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:23  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:21  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE CHDIRI
  CERN PROGLIB# Z265    CHDIRI          .VERSION KERNFOR  4.38  931108
  ORIG. 22/02/91, JZ
  Fortran interface routine to chdir

      ISTAT =  CHDIRF (NAME)

          NAME  the name of the new current working directory
         ISTAT  returns zero if successful
*/
#include <stdio.h>
#ifdef WIN32
#include <direct.h>
# ifndef chdir
#   define chdir _chdir
# endif
#endif
#include "kerngen/fortchar.h"
#if defined(CERNLIB_QX_SC)
int type_of_call chdiri_(fname,lgname)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call chdiri(fname,lgname)
#endif
#if defined(CERNLIB_QXCAPT)
#  ifdef CERNLIB_MSSTDCALL
    int type_of_call CHDIRI(fname, len_fname, lgname)
    int len_fname;
#  else
    int type_of_call CHDIRI(fname,lgname)
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
      char *ptname, *fchtak();
      int  istat, chdir();

/*        get memory and copy NAME terminated  */

      ptname = fchtak(fname,*lgname);
      if (ptname == NULL)           goto bad;

      istat = chdir (ptname);
      free (ptname);
      return istat;

bad:  return -1;
}
/*> END <----------------------------------------------------------*/
