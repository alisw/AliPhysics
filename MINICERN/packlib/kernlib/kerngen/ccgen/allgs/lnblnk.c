/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/10/23 16:25:15  mclareni
 * NT mods, mostly C Fortran interface
 *
 * Revision 1.2  1997/02/04 17:35:00  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:30:06  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:34  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE LNBLNK
  CERN PROGLIB# M432    LNBLNK          .VERSION KERNVMI  1.06  920511
  ORIG. 30/04/92, RDM + JZ

  N = LNBLNK (CHLINE)   find last non-blank character in CHLINE
*/
#ifndef CERNLIB_MSSTDCALL
      int lnblnk_(chline, len)
#else
      int __stdcall LNBLNK(chline, len)
#endif
      char  *chline;
      int   len;
{
      char  *chcur;

      chcur = chline + len;
      while (chcur > chline)
        {  if (*--chcur != ' ')      goto exit; }
      return 0;

exit: return chcur+1 - chline;
}
/*> END <----------------------------------------------------------*/
