/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:38  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:39  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:25  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE RENAMI
  CERN PROGLIB# Z265    RENAMI          .VERSION KERNFOR  4.31  911111
  ORIG. 22/02/91, JZ
  Fortran interface routine to rename

      ISTAT = RENAMEF (FROM, TO)

          FROM  old file name
            TO  new file name
         ISTAT  zero if successful
*/
#include <stdio.h>
#include "kerngen/fortchar.h"
#if defined(CERNLIB_QX_SC)
int type_of_call renami_(frpath, topath, lgfr, lgto)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call renami(frpath, topath, lgfr, lgto)
#endif
#if defined(CERNLIB_QXCAPT)
int type_of_call RENAMI(frpath, topath, lgfr, lgto)
#endif
#if defined(CERNLIB_QMCRY)
      _fcd  frpath,  topath;
#endif
#if !defined(CERNLIB_QMCRY)
      char *frpath, *topath;
#endif
      int  *lgfr, *lgto;
{
      char *ptfr, *ptto, *fchtak();
      int  istat, rename();

      istat = -1;
      ptfr  = fchtak(frpath,*lgfr);
      if (ptfr == NULL)            goto bad;

      ptto  = fchtak(topath,*lgto);
      if (ptto == NULL)            goto pre;

      istat = rename (ptfr, ptto);

      free (ptto);
pre:  free (ptfr);
bad:  return istat;
}
/*> END <----------------------------------------------------------*/
