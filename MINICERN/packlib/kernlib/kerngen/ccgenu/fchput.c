/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:40  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
/*>    ROUTINE FCHPUT
  CERN PROGLIB#         FCHPUT          .VERSION KERNFOR  4.31  911111
  ORIG. 22/02/91, JZ

      Copy a zero-terminated C character string
      to a Fortran character string of length NTEXT,
      return length and blank-fill
*/
#include <stdio.h>
#include "kerngen/fortchar.h"
int fchput(pttext,ftext,lgtext)
      char *pttext;
#if defined(CERNLIB_QMCRY)
      _fcd ftext;
#endif
#if !defined(CERNLIB_QMCRY)
      char *ftext;
#endif
      int  lgtext;
{
      char *utext;
      int  limit, jcol;
      int  nhave;

      limit = lgtext;
      jcol  = 0;
#if defined(CERNLIB_QMCRY)
      utext = _fcdtocp(ftext);
#endif
#if !defined(CERNLIB_QMCRY)
      utext = ftext;
#endif
      if (pttext == NULL)          goto out;

/*--      copy the text to the caller   */
      for (jcol = 0; jcol < limit; jcol++)
      {   if (*pttext == '\0')  break;
          *utext++ = *pttext++;
        }

out:  nhave = jcol;
      for (; jcol < limit; jcol++)   *utext++ = ' ';
      return nhave;
}
/*> END <----------------------------------------------------------*/
