/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:31  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE LNBLNK
  CERN PROGLIB# M507    LNBLNK          .VERSION KERNOS9  1.01  940718
  ORIG. 30/04/92, RDM + JZ
 
  N = LNBLNK (CHLINE)   find last non-blank character in CHLINE
*/
#if defined(CERNLIB_QX_SC)
      int lnblnk_(chline, len)
#endif
#if defined(CERNLIB_QXNO_SC)
      int lnblnk(chline, len)
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
