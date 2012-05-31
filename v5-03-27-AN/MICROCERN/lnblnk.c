/*
 * $Id$
 *
 */
/*>    ROUTINE LNBLNK
  CERN PROGLIB# M432    LNBLNK          .VERSION KERNVMI  1.06  920511
  ORIG. 30/04/92, RDM + JZ

  N = LNBLNK (CHLINE)   find last non-blank character in CHLINE
*/
      int lnblnk_(chline, len)
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
