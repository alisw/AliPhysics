*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:31  mclareni
* Kernlib
*
*
/*>    ROUTINE LSTATI
  CERN PROGLIB# Z265    LSTATI          .VERSION KERNIRT  1.06  930811
C ORIG. 24/03/91, RDM + JZ
  Fortran interface routine to lstat
*/
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "kerngen/fortchar.inc"
#if defined(CERNLIB_QX_SC)
int lstati_(fname, info, lgname, slate)
#endif
#if defined(CERNLIB_QXNO_SC)
int lstati(fname, info, lgname, slate)
#endif
#if defined(CERNLIB_QXCAPT)
int LSTATI(fname, info, lgname, slate)
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

/     if (!istat) {
         info[0] = (int) buf->st_dev;
         info[2] = (int) buf->st_ino;
         info[4] = (int) buf->st_mode;
         info[6] = (int) buf->st_nlink;
         info[8] = (int) buf->st_uid;
         info[10] = (int) buf->st_gid;
         info[12] = (int) buf->st_size;
         info[14] = (int) buf->st_atime;
         info[16] = (int) buf->st_mtime;
         info[18] = (int) buf->st_ctime;
         info[20] = (int) buf->st_blksize;
         info[22] = (int) buf->st_blocks;
         *slate++ = (buf->st_mode & S_IFMT) ^ S_IFREG;
         *slate++ = (buf->st_mode & S_IFMT) ^ S_IFLNK;
         *slate++ = (buf->st_mode & S_IFMT) ^ S_IFDIR;
       };

      free(buf);
out2: free(ptname);
out1: return istat;

}
/*> END <----------------------------------------------------------*/
