/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:49:32  mclareni
 * Kernlib
 *
 */
/*>    ROUTINE LOCF
  CERN PROGLIB# N100    LOCF            .VERSION KERNVMI  1.08  930527
  ORIG. 11/05/93, JZ
*/
unsigned int locf_(iadr)
   char *iadr;
{
   return( ((unsigned) iadr) >> 2 );
}
/*> END <----------------------------------------------------------*/
