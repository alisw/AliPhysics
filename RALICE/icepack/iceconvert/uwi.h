#ifndef _RDMC_UWI_H
#define _RDMC_UWI_H

int rdmc_uwi_FH(const char *s, mcfile *fp); /* file header line */
int rdmc_rhd_uwi(mcfile *fp);             /* reads the format info from the file */
int rdmc_warr_uwi(char *geofile,const array *geo);                      /* (UWI) */
int rdmc_whd_uwi(const mcfile *fp);                  /* write header to UWI file */

/**** in UWI format we need to read the special geometry file separately ***/
/** this geo file is read by the following routine */
int rdmc_rarr_uwi(char *geofile, array *ar);                 /* (UWI ascii file) */

int rdmc_revt_uwi(mcfile *fp, mevt *ev, const array *ar);  /* (uwi/ascii format) */
int rdmc_skipevt_uwi(mcfile *fp);                                       /* (UWI) */
int rdmc_wevt_uwi(const mcfile *fp,const mevt *event, const array *ar); /* (UWI) */


#endif
