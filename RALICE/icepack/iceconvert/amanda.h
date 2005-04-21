#ifndef _RDMC_AMANDA_H
#define _RDMC_AMANDA_H

#define _USE_MATH_DEFINES
#include <math.h>

#include <malloc.h>

#define PID180 (M_PI/180.)

#define SEC_PER_DAY 86400



int rdmc_amanda_V(const char *s, mcfile *fp);

int rdmc_rhd_amanda(mcfile *fp);   /* reads the format info from the file */
int rdmc_rarr_amanda(mcfile *fp, array *ar);   /* (amanda/f2000 ascii file) */

int rdmc_warr_amanda(const mcfile *fp, const array *ar);      /* (amanda) */
int rdmc_wrend_amanda(const mcfile *fp);    /* write end of file (F2000 format) */

#if 0
/* writes the arry of coment lines essentially as they are */
int rdmc_wrcomment_amanda(const mcfile *fp, const char *s);           /* (amanda) */
/* writes an array of lines but puts prefix (e.g. "HI ") before each line */
int rdmc_wrhist_amanda(const mcfile *fp, const char *s, const char *prefix);
#endif

int rdmc_skipevt_amanda(mcfile *fp);                                 /* (amanda) */
int rdmc_revt_amanda(mcfile *fp, mevt *ev, array *ar);/* (amanda/f2000 ascii) */
int rdmc_wevt_amanda(const mcfile *fp, const mevt *ev, const array *ar);/* amanda */


#endif


