/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  1997/10/23 16:25:10  mclareni
 * NT mods, mostly C Fortran interface
 *
 * Revision 1.2  1997/02/04 17:34:22  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:30  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:23  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"

/*>    ROUTINE JMPSET (AREA, ROUT)

  CERN PROGLIB#         JMPSET          .VERSION KERNFOR  4.29  910718

       The function  setjmp  cannot be implemented
            by a Fortran interface routine

      Instead, we provide  IRETN = JMPSET (AREA,ROUT) which
      will dump the environment into AREA and call ROUT.
      IRTN = 0 on return signals normal return from ROUT,
      IRTN = n signals return from JMPLONG (AREA,n)

Usage :
          PROGRAM TOP
              COMMON /JMP/ AREA(32)
              EXTERNAL  XQT

           12 IRTN = JMPSET (AREA,XQT)
              GO TO 12
              END

          SUBROUTINE XQT
              CALL DOWN
              END

          SUBROUTINE DOWN
              COMMON /JMP/ AREA(32)

              IF (HOME)  CALL JMPLONG (AREA,1)
              END
*/
#include <setjmp.h>
#if defined(CERNLIB_QX_SC)
int type_of_call jmpset_(area,ufun)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call jmpset(area,ufun)
#endif
#if defined(CERNLIB_QXCAPT)
int type_of_call JMPSET(area,ufun)
#endif
#if defined(CERNLIB_QSIGJMP)
      sigjmp_buf area;
#endif
#if !defined(CERNLIB_QSIGJMP)
      jmp_buf area;
#endif
#if defined(CERNLIB_QCCINDAD)
      void (type_of_call ** ufun)();
#endif
#if !defined(CERNLIB_QCCINDAD)
      void (type_of_call * ufun)();
#endif
{
      void (type_of_call * unext)();
      int  irtn;

#if defined(CERNLIB_QCCINDAD)
      unext = *ufun;
#endif
#if !defined(CERNLIB_QCCINDAD)
      unext = ufun;
#endif

#if defined(CERNLIB_QSIGJMP)
      irtn = sigsetjmp(area,7);
#endif
#if !defined(CERNLIB_QSIGJMP)
      irtn = setjmp(area);
#endif

      if (irtn != 0)               return irtn;
      (* unext)();
      return 0;
}
/*> END <----------------------------------------------------------*/
