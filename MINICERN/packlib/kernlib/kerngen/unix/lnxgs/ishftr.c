/*
 * $Id$
 *
 * $Log$
 * Revision 1.1.1.1  1996/02/15 17:50:07  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
/*>    ROUTINE ISHFT
  CERN PROGLIB#         ISHFTR          .VERSION KERNLNX  1.02  940511

  Logical right shift by *len (+ve) places
*/
unsigned int ishftr_(arg,len)
unsigned int *arg;
int *len;
{
   return(*arg >> *len);
}
/*> END <----------------------------------------------------------*/
