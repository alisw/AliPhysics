/*>    ROUTINE LSTATI
  CERN PROGLIB# Z265    LSTATI          .VERSION KERNFOR  4.38  931108
C ORIG. 24/03/91, RDM + JZ
  Fortran interface routine to lstat
  Version for Windows NT/Windows 95 by Valery Fine 30/05/96 (fine@vxcern.cern.ch)
*/
#ifdef __STDC__
# undef __STDC__
#endif

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "kerngen/fortranc.h"

#if defined(CERNLIB_QX_SC)
int  type_of_call lstati_(fname, info, lgname, slate)
#endif
#if defined(CERNLIB_QXNO_SC)
int  type_of_call lstati(fname, info, lgname, slate)
#endif
#if defined(CERNLIB_QXCAPT)
int  type_of_call LSTATI(fname, info, lgname, slate)
#endif
      char *fname;
      int  *lgname;
      int  *info;
      int  *slate;
{
      struct _stat *buf;
      char *ptname, *fchtak();
      int  istat, lstat();

      istat  = -1;
      ptname = fchtak(fname,*lgname);
      if (ptname == NULL)          goto out1;

      buf = (struct _stat *) malloc(sizeof (struct _stat));
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
         info[10] = 0;
         info[11] = 0;
         *slate++ = (buf->st_mode & S_IFMT) ^ S_IFREG; 
         *slate++ = 0;
         *slate++ = (buf->st_mode & S_IFMT) ^ S_IFDIR;
       };

      free(buf);
out2: free(ptname);
out1: return istat;

}
/*> END <----------------------------------------------------------*/
