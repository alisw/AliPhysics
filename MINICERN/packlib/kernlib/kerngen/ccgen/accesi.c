/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/10/23 16:25:08  mclareni
 * NT mods, mostly C Fortran interface
 *
 * Revision 1.2  1997/02/04 17:34:13  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:22  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:21  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#ifdef CERNLIB_WINNT
#  include <io.h>
#endif
/*>    ROUTINE ACCESI
  CERN PROGLIB# Z265    ACCESI          .VERSION KERNFOR  4.34  930114
  ORIG. 06/10/92, RDM + JZ
  Fortran interface to access

     R_OK    4   test for read permission
     W_OK    2   test for write permission
     X_OK    1   test for execute (search) permission
     F_OK    0   test for presence of file

     accessible = access(path, mode)
     int accessible;
     char *path;
     int mode;

  access checks the given file path for accessibility according to mode,
  which is an inclusive or of the bits R_OK, W_OK, and X_OK. Specifying
  mode as F_OK (that is, 0) tests whether the directories leading to the
  file can be searched and the file exists.
*/
#include <stdio.h>
#include "kerngen/fortchar.h"
#if defined(CERNLIB_QX_SC)
int type_of_call accesi_(fname, mode, lgname)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call accesi(fname, mode, lgname)
#endif
#if defined(CERNLIB_QXCAPT)
int type_of_call ACCESI(fname,
#ifdef CERNLIB_MSSTDCALL
                       lfname,
#endif 
                              mode, lgname)
#endif
#ifdef CERNLIB_MSSTDCALL
     int  lfname;
#endif 
#if defined(CERNLIB_QMCRY)
      _fcd  fname;
#endif
#if !defined(CERNLIB_QMCRY)
      char *fname;
#endif
      int  *lgname, *mode;
{
      char   *ptf, *fchtak();
      int     istat, umode;

      istat = -1;
      ptf = fchtak(fname, *lgname);
      if (ptf == NULL)                goto exit;

      umode = *mode & 7;
      istat = access(ptf, umode);
      free(ptf);

exit:
      return istat;
}
