/*
 * $Id$
 *
 * $Log$
 * Revision 1.5  1997/11/05 11:05:35  mclareni
 * Remove malloc and free, which caused optimisation problems on NT; should be faster on Unix too.
 *
 * Revision 1.4  1997/10/23 16:33:20  mclareni
 * NT mods
 *
 * Revision 1.3  1997/09/02 14:26:47  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/20 16:41:48  gunter
 * Mods for WNT; ie. transcribe the mods done to stati.c for WNT.
 *
 * Revision 1.1  1996/10/16 12:57:35  cernlib
 * Add cfstat. cfstati is used by cfstat. This uses rfio_stat if CERNLIB_SHIFT
 * is set.
 *
 * Kernlib
 *
 */

#ifdef CERNLIB_WINNT
/*#pragma optimize( "", off ) */
#endif

#include "kerngen/pilot.h"
/*>    ROUTINE STATI
 *  CERN PROGLIB# Z310    CFSTATI
 * ORIG. stolen with mods from stati.c, 11-Oct-96; GF.
 *  Routine used by cfstat; interface to stat or rfio_stat ( if shift software 
 *    is in use)
*/
#include <stdio.h>
#if defined(CERNLIB_QMVAX)||defined(CERNLIB_QMOS9)
#include <types.h>
#include <stat.h>
#endif
#if (!defined(CERNLIB_QMVAX))&&(!defined(CERNLIB_QMOS9))
#ifndef WIN32
#include <sys/types.h>
#include <sys/stat.h>
#include "kerngen/cf_xaft.h"
#else
# include <sys\types.h>
# include <sys\stat.h>
#endif
#include "kerngen/fortchar.h"
#endif
#include "kerngen/fortranc.h"
#if defined(CERNLIB_QX_SC)
int type_of_call cfstati_(fname, info, lgname)
#endif
#if defined(CERNLIB_QXNO_SC)
int type_of_call cfstati(fname, info, lgname)
#endif
#if defined(CERNLIB_QXCAPT)
#ifndef CERNLIB_MSSTDCALL
  int type_of_call CFSTATI(fname, info, lgname)
#else
  int type_of_call CFSTATI(fname, lfname, info, lgname )
  int lfname;
#endif
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

      struct stat buf;
      char *ptname, *fchtak();
      int  istat=-1, stat();

      ptname = fchtak(fname,*lgname);
      if (ptname == NULL) return -1;


#if defined(CERNLIB_PROJSHIFT)

      istat = rfio_stat(ptname, &buf);

#else

      istat = stat(ptname, &buf);

#endif

      if (!istat) {
         info[0] = (int) buf.st_dev;
         info[1] = (int) buf.st_ino;
         info[2] = (int) buf.st_mode;
         info[3] = (int) buf.st_nlink;
         info[4] = (int) buf.st_uid;
         info[5] = (int) buf.st_gid;
         info[6] = (int) buf.st_size;
         info[7] = (int) buf.st_atime;
         info[8] = (int) buf.st_mtime;
         info[9] = (int) buf.st_ctime;
#if defined(CERNLIB_QMDOS)||defined(CERNLIB_QMVAX)||defined(CERNLIB_QMOS9) \
  ||defined(CERNLIB_WINNT)
         info[10] = 0;
         info[11] = 0;
#else
         info[10] = (int) buf.st_blksize;
         info[11] = (int) buf.st_blocks;
#endif
       };

      free(ptname);
      return istat;

}
/*> END <----------------------------------------------------------*/
