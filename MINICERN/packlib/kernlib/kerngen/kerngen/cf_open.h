/*
* $Id$
*
* $Log$
* Revision 1.2  1997/02/04 17:35:35  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:30:25  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:17  mclareni
* Kernlib
*
*
*
* cf#open.inc
*/
#if defined(CERNLIB_QMAPO)
#include <sys/file.h>        /*  Apollo                     */
#elif defined(CERNLIB_QMAMX)
#include <sys/types.h>       /*  AMX                        */
#include <sys/fcntl.h>
#include <sys/file.h>
#elif defined(CERNLIB_QMCRY)||defined(CERNLIB_QMHPX)||defined(CERNLIB_QMIBX)||defined(CERNLIB_QMSUN)
#include <fcntl.h>           /*  CRAY HPX IBX SUN           */
#elif defined(CERNLIB_QMIRT)||defined(CERNLIB_QMIRTD)
#include <fcntl.h>           /*  IRT                        */
#include <sys/mode.h>
#elif defined(CERNLIB_QMVAX)
#include <file.h>            /*  VAX/VMS                    */
#elif defined(CERNLIB_QMVMI)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>           /*  VMI Decstation             */
#elif defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
#ifdef __GNUC__
#include <sys/file.h>
#else
# ifdef WIN32
#  ifdef __STDC__
#   undef __STDC__
#  endif
#  include <io.h>
#  include <stdio.h>
#  include <fcntl.h>
# endif
#endif
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/fcntl.h>       /*  default Posix              */
#endif
