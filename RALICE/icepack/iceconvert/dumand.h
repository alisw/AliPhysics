#ifndef _RDMC_DUMAND_H
#define _RDMC_DUMAND_H

#include <malloc.h>

int rdmc_a_V(const char *s, mcfile *fp);                     /* ascii lines */
int rdmc_rhd_ascii(mcfile *fp);           /* reads the format info from the file */
int rdmc_rarr_ascii(mcfile *fp, array *ar);            /* (dumand/cw ascii file) */
int rdmc_warr_ascii(const mcfile *fp, const array *ar);        /* (dumand) */
int rdmc_wrcomment_ascii(const mcfile *fp, const char *s);   /* (dumand) */

int rdmc_revt_ascii(mcfile *fp, mevt *ev, const array *ar); /* (dumand/cw ascii) */

int rdmc_skipevt_ascii(mcfile *fp);                             /* (dumand) */

int rdmc_wevt_ascii(const mcfile *fp, const mevt *ev, const array *ar);

#endif


