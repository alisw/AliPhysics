/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/09/02 14:26:40  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/04 17:34:50  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:47  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:28  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE UNLINI
  CERN PROGLIB# Z265    UNLINI          .VERSION KERNFOR  4.31  911111
  ORIG. 15/10/91, JZ
  Fortran interface routine to unlink

      ISTAT =  UNLINKF (NAME)

          NAME  the name of the file to be deleted
         ISTAT  returns zero if successful
*/
#include <stdio.h>
#include "kerngen/fortchar.h"
#if defined(CERNLIB_QX_SC)
int type_of_call unlini_(fname,lgname)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call unlini(fname,lgname)
#endif
#if defined(CERNLIB_QXCAPT)
#  ifndef CERNLIB_MSSTDCALL
     int type_of_call UNLINI(fname,lgname)
#  else
     int type_of_call UNLINI(fname,lfname,lgname)
     int lfname;
#  endif
#endif
#if defined(CERNLIB_QMCRY)
      _fcd  fname;
#endif
#if !defined(CERNLIB_QMCRY)
      char *fname;
#endif
      int  *lgname;
{
      char *ptname, *fchtak();
      int  istat, unlink();

/*        get memory and copy NAME terminated  */

      ptname = fchtak(fname,*lgname);
      if (ptname == NULL)           goto bad;

      istat = unlink (ptname);
      free (ptname);
      return istat;

bad:  return -1;
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_UNLINKF
#undef CERNLIB_TCGEN_UNLINKF
#endif
