/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/10/23 16:25:11  mclareni
 * NT mods, mostly C Fortran interface
 *
 * Revision 1.2  1996/09/20 14:51:15  cernlib
 * Linux added
 *
 * Revision 1.1.1.1  1996/02/15 17:49:20  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMHPX)||defined(CERNLIB_QMOS9)
#include "hpxgs/lenocc.c"
#elif defined(CERNLIB_QMIRT)||defined(CERNLIB_QMIRTD)
#include "irtgs/lenocc.c"
#elif 1
/*>    ROUTINE LENOCC
  CERN PROGLIB# M507    LENOCC          .VERSION KERNVMI  1.06  920511
  ORIG. 30/04/92, RDM + JZ

  N = LENOCC (CHLINE)   find last non-blank character in CHLINE
*/
#if defined(CERNLIB_QXNO_SC)
      int type_of_call lenocc(chline, len)
#elif defined(CERNLIB_QXCAPT)
      int type_of_call LENOCC(chline, len)
#else
      int type_of_call lenocc_(chline, len)
#endif
      char  *chline;
      int   len;
{
      static unsigned int blnk = 0x20202020;
      unsigned int *wdcur;
      char  *chcur;
      int   ntail, i;

      chcur = chline + len;
      if (len <= 24)               goto small;

/* ----        handle long string             */

/*        look at the last ntail characters   */

      ntail = ((int)chcur & 3);
      for (i = ntail; i > 0; i--)
        { if (*--chcur != ' ')       goto exit; }

/*        look for trailing blank words   */

      wdcur = (unsigned int*) (chcur-4);
      while (wdcur >= (unsigned int*)chline )
        {  if (*wdcur != blnk)   break;  wdcur--; }

/*        find last non-blank character   */

      chcur = (char*) (wdcur+1);
      while (chcur > chline)
        {  if (*--chcur != ' ')      goto exit; }
      return 0;

exit: return chcur+1 - chline;

/* ----        handle short string            */

small:
      while (chcur > chline)
        {  if (*--chcur != ' ')      goto exit; }
      return 0;
}
/*> END <----------------------------------------------------------*/
#endif
