/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  1997/02/04 17:34:50  mclareni
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

/*>    ROUTINE VXINVB
  CERN PROGLIB#M434     VXINVB          .VERSION KERNFOR  4.42  951011
  ORIG. 10/07/95, JZ

      CALL VXINVB (IXV,N)

      byte inversion N words of IXV, in situ
*/

#if defined(CERNLIB_QX_SC)
void type_of_call vxinvb_(ixv, n)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call vxinvb(ixv, n)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call VXINVB(ixv, n)
#endif
      int  *ixv, *n;
{
      int  limit, jloop;
      int  in;

      limit = *n;

/*--          swop:   1 | 2 | 3 | 4
                to:   4 | 3 | 2 | 1     */

      for (jloop = 0; jloop < limit; jloop++)
    { in = ixv[jloop];
      ixv[jloop] =
          ((in >> 24) & 0x000000ff) |
          ((in >>  8) & 0x0000ff00) |
          ((in <<  8) & 0x00ff0000) |
          ((in << 24) & 0xff000000);
        }
      return;
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_VXINVB
#undef CERNLIB_TCGEN_VXINVB
#endif
