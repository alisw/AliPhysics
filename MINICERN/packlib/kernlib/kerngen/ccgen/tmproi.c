/*
* $Id$
*
* $Log$
* Revision 1.3  1997/10/23 16:25:12  mclareni
* NT mods, mostly C Fortran interface
*
* Revision 1.2  1997/02/04 17:34:48  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:29:46  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:28  mclareni
* Kernlib
*
*/
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMDOS)
#include "wntgs/tmproi.c"
#else
/*>    ROUTINE TMPROI
  CERN PROGLIB#         TMPROI          .VERSION KERNFOR  4.39  940228
  ORIG. 30/05/91, JZ
  Fortran interface routine to print a prompt string
      CALL TMPRO (TEXT)
*/
#ifdef WIN32
#include <io.h>
#endif
#include <stdio.h>
#include "kerngen/fortchar.h"
#if defined(CERNLIB_QX_SC)
void tmproi_(ftext, lgtext)
#endif
#if defined(CERNLIB_QXNO_SC)
void tmproi(ftext, lgtext)
#endif
#if defined(CERNLIB_QXCAPT)
void TMPROI(ftext, lgtext)
#endif
#if defined(CERNLIB_QMCRY)
      _fcd ftext;
#endif
#if !defined(CERNLIB_QMCRY)
      char *ftext;
#endif
      int  *lgtext;
{
      write (1, ftext, *lgtext);
      return;
}
/*> END <----------------------------------------------------------*/
#endif
