/*
* $Id$
*
* $Log$
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
#include "kerngen/fortchar.inc"
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

