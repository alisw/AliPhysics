/*
* $Id$
*
* $Log$
* Revision 1.2  1997/02/04 17:34:37  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:29:38  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:29  mclareni
* Kernlib
*
*/
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE READLNI
  CERN PROGLIB# Z265    READLNI         .VERSION KERNFOR  4.36  930602
  ORIG. 24/03/93, JZ
  Fortran interface routine to readlink

      NCH = READLNF (PATH, TEXT*)

          PATH  the path name of the link
          TEXT  returns its value
          NCH   returns the length of the value
*/
#include <stdio.h>
#include "kerngen/fortchar.h"
#if defined(CERNLIB_QX_SC)
int type_of_call readlni_(fname, ftext, lgtext, lgname)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call readlni(fname, ftext, lgtext, lgname)
#endif
#if defined(CERNLIB_QXCAPT)
int type_of_call READLNI(fname, ftext, lgtext, lgname)
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
      int  readlink();
      int  fchput();
      int  nch;

      nch = -1;
      ptname = fchtak(fname,*lgname);
      if (ptname == NULL)          goto out;
      nch = readlink (ptname,ftext,*lgtext);
      free(ptname);

out:  return nch;
}
/*> END <----------------------------------------------------------*/
