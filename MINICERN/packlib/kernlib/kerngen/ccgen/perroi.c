/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/09/02 14:26:38  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/04 17:34:36  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:37  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:25  mclareni
 * Kernlib
 *
 */ 
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"

/*>    ROUTINE PERROI
  CERN PROGLIB# Z265    PERROI          .VERSION KERNFOR  4.31  911111
  ORIG. 22/02/91, JZ
  Fortran interface routine to perror

      CALL PERRORF (TEXT)

          TEXT  the text to be printed before the error message
*/
#include <stdio.h>
#include "kerngen/fortchar.h"
#if defined(CERNLIB_QX_SC)
void type_of_call perroi_(ftext, lgtext)
#endif
#if defined(CERNLIB_QXNO_SC)
void type_of_call perroi(ftext, lgtext)
#endif
#if defined(CERNLIB_QXCAPT)
void type_of_call PERROI(ftext, 
#  ifdef CERNLIB_MSSTDCALL
                        len_ftext,
#  endif
                        lgtext)
#endif
#if defined(CERNLIB_QMCRY)
      _fcd ftext;
#endif
#if !defined(CERNLIB_QMCRY)
      char *ftext;
#endif

#ifdef CERNLIB_MSSTDCALL
      int  len_ftext;
#endif
      int  *lgtext;
{
      char *pttext, *fchtak();

      pttext = fchtak(ftext,*lgtext);
      perror (pttext);
      if (pttext != NULL)   free (pttext);
      return;
}
/*> END <----------------------------------------------------------*/
#ifdef CERNLIB_TCGEN_PERRORF
#undef CERNLIB_TCGEN_PERRORF
#endif
