/*
 * $Id$
 *
 * $Log$
 * Revision 1.4  1997/09/02 14:26:39  mclareni
 * WINNT correction
 *
 * Revision 1.3  1997/02/04 17:34:45  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.2  1996/10/11 07:49:51  cernlib
 * Sgi now has  st_blksize and st_blocks, so store and return these.
 * Use #else  instead of complicated #if
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:44  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:27  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
/*>    ROUTINE STATI
  CERN PROGLIB# Z265    STATI           .VERSION KERNFOR  4.40  940929
C ORIG. 14/03/91, RDM
  Fortran interface routine to stat
*/
#include <stdio.h>
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMOS9)
#include <types.h>
#include <stat.h>
#endif
#if (!defined(CERNLIB_QMVAX))&&(!defined(CERNLIB_QMOS9))
#include <sys/types.h>
#include <sys/stat.h>
#include "kerngen/fortchar.h"
#endif
#include "kerngen/fortranc.h"
#if defined(CERNLIB_QX_SC)
int type_of_call stati_(fname, info, lgname)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call stati(fname, info, lgname)
#endif
#if defined(CERNLIB_QXCAPT)
# ifndef CERNLIB_MSSTDCALL
    int type_of_call STATI(fname, info, lgname)
# else
    int type_of_call STATI(fname, lfname, info, lgname)
    int lfname;
# endif
#endif
#if defined(CERNLIB_QMCRY)
      _fcd  fname;
#endif
#if !defined(CERNLIB_QMCRY)
      char *fname;
#endif
      int  *lgname;
      int  *info;
{
#ifndef WIN32
      struct stat *buf;
#else
      struct _stat *buf;
#endif
      char *ptname, *fchtak();
      int  istat, stat();

      istat  = -1;
      ptname = fchtak(fname,*lgname);
      if (ptname == NULL)          goto out1;

#ifndef WIN32
      buf = (struct stat *) malloc(sizeof (struct stat));
#else
      buf = (struct _stat *) malloc(sizeof (struct _stat));
#endif

      if (buf == NULL)             goto out2;

      istat = stat(ptname, buf);

      if (!istat) {
         info[0] = (int) buf->st_dev;
         info[1] = (int) buf->st_ino;
         info[2] = (int) buf->st_mode;
         info[3] = (int) buf->st_nlink;
         info[4] = (int) buf->st_uid;
         info[5] = (int) buf->st_gid;
         info[6] = (int) buf->st_size;
         info[7] = (int) buf->st_atime;
         info[8] = (int) buf->st_mtime;
         info[9] = (int) buf->st_ctime;
#if defined(CERNLIB_QMDOS)||defined(CERNLIB_QMVAX)||defined(CERNLIB_QMOS9)||defined(CERNLIB_WINNT)
         info[10] = 0;
         info[11] = 0;
#else
         info[10] = (int) buf->st_blksize;
         info[11] = (int) buf->st_blocks;
#endif
       };

      free(buf);
out2: free(ptname);
out1: return istat;

}
/*> END <----------------------------------------------------------*/
