/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/09/02 14:26:37  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/04 17:34:31  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:35  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:24  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

#if defined(CERNLIB_MSSTDCALL) && defined(CERNLIB_LOCF_CHARACTER)
# define Dummy2LocPar  ,_dummy
# define DummyDef     int _dummy;
#else
# define Dummy2LocPar  
# define DummyDef
#endif


/*>    ROUTINE LOCB
  CERN PROGLIB# N101    LOCB            .VERSION KERNFOR  4.36  930602
*/
#if defined(CERNLIB_QX_SC)
int type_of_call locb_(iadr Dummy2LocPar)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call locb(iadr Dummy2LocPar)
#endif
#if defined(CERNLIB_QXCAPT)
int type_of_call LOCB(iadr Dummy2LocPar)
#endif
   char *iadr;
#ifdef DummyDef
   DummyDef
#endif
{
   return( (int) iadr );
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_LOCB
#undef CERNLIB_TCGEN_LOCB
#endif
