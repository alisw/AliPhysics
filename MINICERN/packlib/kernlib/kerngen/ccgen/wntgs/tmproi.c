/*
* $Id$
*
* $Log$
* Revision 1.1  1997/02/04 17:35:07  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1  1996/02/15 17:49:32  mclareni
* Kernlib
*
*/
/*>    ROUTINE TMPROI
  CERN PROGLIB#         TMPROI          .VERSION KERNFOR  4.36  930602
  ORIG. 30/05/91, JZ
  Fortran interface routine to print a prompt string
      CALL TMPRO (TEXT)
*/
#ifdef WIN32
#include <io.h>
#endif
#include <stdio.h>
#include "kerngen/fortranc.h"

#if defined(CERNLIB_QX_SC)
void type_of_call tmproi_(ftext, lgtext)
#endif
#if defined(CERNLIB_QXNO_SC)
void  type_of_call tmproi(ftext, lgtext)
#endif
#if defined(CERNLIB_QXCAPT)
void  type_of_call TMPROI(ftext, lgtext)
#endif
      char *ftext;
      int  *lgtext;
{
      write (1, ftext, *lgtext);
      return;
}
/*> END <----------------------------------------------------------*/

