/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:31  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE LOCF
  CERN PROGLIB# N100    LOCF            .VERSION KERNIRT  1.06  930811
*/
*    Number of ADdress Units Per Word
#define NADUPW 8   /* Number of ADdress Units Per Word */
#define LADUPW 3   /* Logarithm base 2 of ADdress Units Per Word */
#if defined(CERNLIB_QX_SC)
unsigned int locf_(iadr)
#endif
#if defined(CERNLIB_QXNO_SC)
unsigned int locf(iadr)
#endif
#if defined(CERNLIB_QXCAPT)
unsigned int LOCF(iadr)
#endif
   char *iadr;
{
   return( ((unsigned) iadr) >> LADUPW );
}
/*> END <----------------------------------------------------------*/
