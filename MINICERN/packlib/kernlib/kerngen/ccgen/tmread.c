/*
* $Id$
*
* $Log$
* Revision 1.2  1997/02/04 17:34:49  mclareni
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

/*>    ROUTINE TMREAD
  CERN PROGLIB#         TMREAD          .VERSION KERNFOR  4.37  930715
  ORIG. 20/07/90, JZ
      read the next line from stdin :
      CALL TMREAD (MAXCH, LINE, NCH, ISTAT)
          MAXCH   maxim. # of characters into LINE
          NCH     actual # of characters read into LINE
          ISTAT   status return, zero : OK  -ve : EoF
*/
#include <stdio.h>
#if defined(CERNLIB_QMVAX)
#include descrip
#endif
#if defined(CERNLIB_QX_SC)
void type_of_call tmread_(alim, cols, anch, astat)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call tmread(alim, cols, anch, astat)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call TMREAD(alim, cols, anch, astat)
#endif
#if defined(CERNLIB_QMVAX)
      struct dsc$descriptor_s  *cols;
#endif
#if !defined(CERNLIB_QMVAX)
      char *cols;
#endif
      int  *alim, *anch, *astat;
{
      char *ubuf;
      int ch, jcol, lim;

#if defined(CERNLIB_QMVAX)
      ubuf = cols->dsc$a_pointer;
#endif
#if !defined(CERNLIB_QMVAX)
      ubuf = cols;
#endif

/*--      read the text   */
      lim  = *alim;
      jcol = 0;
      while (lim-- > 0)
      {   ch = getchar();
          if (ch == EOF)           goto endf;
          if (ch == '\n')          goto endl;
          *ubuf++ = ch;
          jcol = jcol + 1;
       }
/*        discard excess characters   */
loop: ch = getchar();
      if (ch == '\n')          goto endl;
      if (ch != EOF)           goto loop;

endf: *anch  = jcol;
      *astat = -1;
      clearerr(stdin);
      return;

endl: *anch  = jcol;
      *astat = 0;
      return;
}
/*> END <----------------------------------------------------------*/
