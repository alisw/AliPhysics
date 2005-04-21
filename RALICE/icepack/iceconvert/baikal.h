
#ifndef _RDMC_BAIKAL_H
#define _RDMC_BAIKAL_H

int rdmc_rarr_baikal_mc(mcfile *fp,array *ar);         /* (baikal mc/data file) */

int rdmc_revt_baikal(mcfile *fp, mevt *ev, const array *ar);        /* (baikal) */
int rdmc_revt_baikal_mc(mcfile *fp, mevt *ev, const array *ar);  /* (baikal mc) */
int rdmc_revt_baikal_data(mcfile *fp, mevt *ev, const array *ar);/* (baikal data)*/
int rdmc_revt_baikal_zeu(mcfile *fp, mevt *ev, const array *ar);/* (zeuthen format) */
int rdmc_revt_baikal_mos(mcfile *fp, mevt *ev, const array *ar); /* (moscow format) */

int rdmc_skipevt_baikal_mc(mcfile *fp);                             /* (baikal mc/data) */

int rdmc_warr_baikal_mc(mcfile *fp, const array *ar);               /* (baikal mc/data) */

int rdmc_wevt_baikal_mc(mcfile *fp, const mevt *ev, const array *ar);   /* (baikal mc) */
int rdmc_wevt_baikal_data(mcfile *fp, const mevt *ev, const array *ar);  /* (baikal data) */
int rdmc_wevt_baikal_mos(mcfile *fp, const mevt *ev, const array *ar);  /* (moscow format) */

int rdmc_wrcomment_baikal_mc(const mcfile *fp, const char *s); /* (baikal, just a dummy) */

int rdmc_rhd_baikal_mc(mcfile *fp);  /* reads the format info from the file */


long rdmc_mcgeti(mcfile *fp);                            /* reads a int (4 byte) */
long rdmc_mcputi(long i, mcfile *fp);                   /* writes a int (4 byte) */
double rdmc_mcgetf(mcfile *fp);                /* reads a float (4 byte) from fp */
double rdmc_mcputf(double r, mcfile *fp);       /* writes a float (4 byte) to fp */
int rdmc_mcseek(mcfile *fp, long n);                /* seeks over n 4-byte-words */
long rdmc_swap(long i);                                       /* swaps the bytes */


#endif
