/*
* $Id$
*
* $Log$
* Revision 1.2  1997/02/04 17:34:35  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:29:36  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:29  mclareni
* Kernlib
*
*/
#include "kerngen/pilot.h" 
#if defined(CERNLIB_WINNT)
#include "wntgs/lstati.c"
#elif defined(CERNLIB_QMDOS)
#include "dosgs/lstati.c"
#elif defined(CERNLIB_QMIRTD)
#include "irtdgs/lstati.c"
#else
/*>    ROUTINE LSTATI
  CERN PROGLIB# Z265    LSTATI          .VERSION KERNFOR  4.39  940228
C ORIG. 24/03/91, RDM + JZ
  Fortran interface routine to lstat
*/
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "kerngen/fortchar.h"
#if defined(CERNLIB_QX_SC)
int lstati_(fname, info, lgname, slate)
#endif
#if defined(CERNLIB_QXNO_SC)
int lstati(fname, info, lgname, slate)
#endif
#if defined(CERNLIB_QXCAPT)
int LSTATI(fname, info, lgname, slate)
#endif
#if defined(CERNLIB_QMCRY)
      _fcd  fname;
#endif
#if !defined(CERNLIB_QMCRY)
      char *fname;
#endif
      int  *lgname;
      int  *info;
      int  *slate;
{
      struct stat *buf;
      char *ptname, *fchtak();
      int  istat, lstat();

      istat  = -1;
      ptname = fchtak(fname,*lgname);
      if (ptname == NULL)          goto out1;

      buf = (struct stat *) malloc(sizeof (struct stat));
      if (buf == NULL)             goto out2;

      istat = lstat(ptname, buf);

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
#if (!defined(CERNLIB_QMSGI))&&(!defined(CERNLIB_QMDOS))
         info[10] = (int) buf->st_blksize;
         info[11] = (int) buf->st_blocks;
#endif
#if defined(CERNLIB_QMSGI)||defined(CERNLIB_QMDOS)
         info[10] = 0;
         info[11] = 0;
#endif
         *slate++ = (buf->st_mode & S_IFMT) ^ S_IFREG;
         *slate++ = (buf->st_mode & S_IFMT) ^ S_IFLNK;
         *slate++ = (buf->st_mode & S_IFMT) ^ S_IFDIR;
       };

      free(buf);
out2: free(ptname);
out1: return istat;

}
/*> END <----------------------------------------------------------*/
#endif
