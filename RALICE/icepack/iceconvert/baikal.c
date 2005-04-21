/*
 * Read/write files in the BAIKAL formats
 * Zeuthen: Krabi/Wischnewski
 * Moscow: an older Nora/Belolaptikov format (only reading)
 */

char *baikal_rdmc_cvsid = 
"$Header: /net/local/cvsroot/siegmund/rdmc/baikal.c,v 1.20 2000/02/07 13:02:17 ole Exp $";

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

////////////#include <unistd.h>


#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>

#include "rdmc.h"

#include "rdmc_local.h"
#include "baikal.h"                /* function prototyping */


#ifndef DEBUG
#define DEBUG 0                       /* DEBUG 1 writes some info to stderr */
#endif

#define ERRLINE (-__LINE__)

#define PID180 (M_PI/180.)

#ifdef BAIKAL_BIN_F

#define BAIKAL_PMTPERCHANNEL 2




static void  rdmc_add_fitdef_baikal(array *ar);


/****************************************************************************/
/* The function rarr_mc() reads the header of a baikal like file            */
/****************************************************************************/

int rdmc_rarr_baikal_mc(mcfile *fp, array *ar)
{
  
  int nw;                                       /* number of words to read */
  int ptr;                    /* pointer to the current word in the header */
  int i,j;                                               /* loop variables */
  int *iclust;              /* pointer to current channel number in string */
  
  rdmc_init_array(ar);                                  /* reset the array */

/**************************** the array info *******************************/
  ar->nrun=fp->info.bai.nrun;

  ptr = 0;
  if (fp->info.bai.fortran_recs) 
    rdmc_mcgeti(fp);                              /* skip fortran record header */
  nw = rdmc_mcgeti(fp);                                /* 1: number of words    */

  if (nw <= 2) {             /* too less words: nonreconstructed data file */
    ar->id = NT_36;                               /* default: baikal array */
    ar->nch = 0;                              /* no telescope info aviable */
    ar->nstr = 0;
    rdmc_mcseek(fp, nw - ptr);                        /* read the unused fields */
    if (fp->info.bai.fortran_recs)
      rdmc_mcgeti(fp);                                          /* unused dummy */
    return 0;                                                /* and return */
  } /* if nw <= 2 */

  rdmc_mcseek(fp,2);                                      /* here unused values */
  ar->nch = rdmc_mcgetf(fp);                           /* 3: number of channels */
  ar->nstr = rdmc_mcgetf(fp);                          /* 4: number of strings  */
  if (ar->nch <= 18)
    ar->id = NT_36;
  else if (ar->nch <= 36)
    ar->id = NT_72;
  else if (ar->nch <= 48)
    ar->id = NT_96;
  else if (ar->nch <= 72)
    ar->id = NT_144;
  else if (ar->nch <= 96)
    ar->id = NT_192;
  else
    ar->id = NT_200;
  ptr += 4;                                                /* 3 words read */
  
  if (nw < ptr + ar->nch * 5) return ERRLINE;            /* invalid format */
  if (ar->nch > RDMC_MAXCHANNELS) return ERRLINE;
  if (ar->nch <= 0) return ERRLINE;
  if (ar->nstr <= 0) return ERRLINE;

  iclust = (int *)calloc(ar->nstr,sizeof(int)); 
                                         /* alloc mem for channel nr array */

/*************************** Channel position and parameters ***************/

  for (i = 0; i < ar->nch; i++) {
    ar->x[i] = rdmc_mcgetf(fp)/100.0;                            /* coordinates */
    ar->y[i] = rdmc_mcgetf(fp)/100.0;                       /* in meters of the */
    ar->z[i] = rdmc_mcgetf(fp)/100.0;                           /* i-th channel */
    ar->costh[i] = rdmc_mcgetf(fp);                                /* cos theta */
    ar->str[i] = rdmc_mcgetf(fp);                              /* string number */
    ar->type[i] = QUASAR;                              /* default pmt type */
    ar->clust[i] = iclust[ar->str[i] - 1]++;             /* channel number */
  } /* for i */
  ptr += 5 * ar->nch;
  ar->is_calib.geo=1;

  free (iclust);                       /* free the mem for the ch nr array */

/****************************** Pmt parameters *****************************/

  if (nw < ptr + ar->nch * 2 * BAIKAL_PMTPERCHANNEL) return ERRLINE;
  for (i = 0; i < ar->nch; i++) {
    ar->thresh[i] = rdmc_mcgetf(fp);                  /* threshold for each pmt */
    ar->sensit[i] = rdmc_mcgetf(fp);                 /* sensitivity of each pmt */
    for (j = 1; j < BAIKAL_PMTPERCHANNEL; j++) {
      rdmc_mcgetf(fp);                                /* threshold for each pmt */
      rdmc_mcgetf(fp);                               /* sensitivity of each pmt */
    } /* for j */
    ar->serial[i] = 0;                                 /* no serial number */
  } /* for i */
  ptr += ar->nch * 2 * BAIKAL_PMTPERCHANNEL;

  rdmc_mcseek(fp, nw - ptr);                      /* read the last unused fields */

  if (fp->info.bai.fortran_recs)
    rdmc_mcgeti(fp);                                            /* unused dummy */

  rdmc_add_fitdef_baikal(ar);

  return 0;

} /* function rarr_mc() */

/****************************************************************************/
/* revt_baikal() tries to read the different BAIKAL dialects                */
/****************************************************************************/
int rdmc_revt_baikal(mcfile *fp, mevt *ev, const array *ar)
{
  long filepos;                         /* shows the current file position */
  int r;

  filepos = ftell(fp->fp);                    /* store the file position */

  if (fp->info.bai.mc) {                   /* if it is (probably) a mc file */
    r = rdmc_revt_baikal_mc(fp,ev,ar);                    /* try to read as an mc event */
    if (r < EOF) {
      fseek(fp->fp,filepos,SEEK_SET);        /* set to old file position */
      r = rdmc_revt_baikal_data(fp,ev,ar);             /* try to read as data file */
      if (r == 0)                               /* if it was successfull */
	fp->info.bai.mc = 0;                      /* set file type to data */
      else                                                     /* if not */
	fseek(fp->fp,filepos,SEEK_SET);         /* restore file position */
    } /* if r == ILF */
  } /* if mc file */
  
  else {                              /* if it is (probably) a data file */
    r = rdmc_revt_baikal_data(fp,ev,ar);             /* try to read as data event */
    if (r < EOF) {                                        /* if it fails */
      fseek(fp->fp,filepos,SEEK_SET);        /* set to old file position */
      r = rdmc_revt_baikal_mc(fp,ev,ar);           /* try to read as an mc event */
      if (r == 0)                               /* if it was successfull */
	fp->info.bai.mc = 1;                     /* set file type to mc */
      else                                                     /* if not */
	fseek(fp->fp,filepos,SEEK_SET);         /* restore file position */
    } /* if r == ILF */
  } /* if data file */
  
  fp->fpos = ftell(fp->fp);         /* store file position for debugging */

  /* now update stuff which was not included in the data */
  /* now uopdate the channel and srtring reference of the event */
  ev->nch = rdmc_count_nch(ev);         /* calc the number of hit channels */
  rdmc_fill_mhit_str(ev,ar);            /* fill the mhit->str fields */
  ev->nstr = rdmc_count_nstr(ev);      /* calc the number of hit channels */



  return r;

} /* revt_baikal() */

/****************************************************************************/
/* revt_mc() reads the next *mc*-event                                      */
/* it is just a switch to the correct routines (moscow/zeuthen format)      */
/****************************************************************************/

int rdmc_revt_baikal_mc(mcfile *fp, mevt *ev, const array *ar)
{
  switch(fp->fmajor) {
  case BAIKAL_OLD:
  case BAIKAL_NEW:
    return rdmc_revt_baikal_zeu(fp,ev,ar);
  case BAIKAL_MOS:
    return rdmc_revt_baikal_mos(fp,ev,ar);
  default:
    return ERRLINE;
  }
  
} /* function revt_mc() */

/****************************************************************************/
/* revt_zeu() reads the zeuthen formats (BAIKAL_OLD,BAIKAL_NEW)             */
/****************************************************************************/

int rdmc_revt_baikal_zeu(mcfile *fp, mevt *ev, const array *ar)
{

  int nw;                                       /* number of words to read */
  int ptr = 0;                /* pointer to the current word in the header */
  int i;                                                  /* loop variable */

  rdmc_clear_mevt(ev);                                  /* clear old event info */

  if (fp->info.bai.fortran_recs)
    rdmc_mcseek(fp, 1);                           /* skip fortran record header */

  nw = rdmc_mcgeti(fp) - 1;                    /* number of words in the record */
  if (feof(fp->fp)) return EOF;                       /* file end reached? */

  if (nw < 8)  return ERRLINE;

/***************** read the first (common) informations ********************/

  ev->t_offset = 0.0;                     /* time offset in baikal is zero */

  if (fp->fmajor == BAIKAL_NEW) {                    /* for the new format */
    i = rdmc_mcgeti(fp);
    ev->nhits = i % 1000;           /* number of hits (higher part is nrun)*/
    ev->nrun = i/1000;
    ev->enr = rdmc_mcgeti(fp);                                  /* event number */
    fp->info.bai.enr = ev->enr;
    ptr += 2;
  } /* if new format */
  else {                                             /* for the old format */
    ev->nrun = ar->nrun;
    ev->enr = ++(fp->info.bai.enr);                         /* set dummies */
  } /* if old format */
  
/************ allocate memory for generated track and fill it **************/

  { 
    mtrack gen;
    rdmc_init_mtrack(&gen);
    gen.costh = rdmc_mcgetf(fp);                   /* cos theta + 10 * nhits */

    if (fp->fmajor == BAIKAL_OLD) {                    /* for the old format */
      ev->nhits = (fabs(gen.costh)+2.0)/10;      /* fill the nr of hits */
      if (ev->nhits == 0) return ERRLINE;         /* there must be >= 1 hit! */
    } /* if old format */
    else                                               /* for the new format */
      if ((int)((fabs(gen.costh)+2.0)/10.0) != ev->nhits) return ERRLINE; 
                                                         /* test for nhits */ 
    gen.costh = fmod(gen.costh,10.);
                                           /* this is now really cos theta */
    gen.phi = rdmc_mcgetf(fp)*PID180;                          /* phi in rad */
    gen.x = rdmc_mcgetf(fp)/100.0;                               /* read the */
    gen.y = rdmc_mcgetf(fp)/100.0;                            /* coordinates */
    gen.z = rdmc_mcgetf(fp)/100.0;                           /* of the track */
    gen.e = rdmc_mcgetf(fp)*1000.0;                         /* Energy in MeV */
    gen.id = MUON_PLUS;                           /* we generated muons */
    gen.tag = 1;
    gen.length = RDMC_BIG;                          /* infinity track length */
    rdmc_tau_tr(&gen);                      /* calc the direction cosinuus */
    rdmc_add_gen(ev,&gen,0);
  }
  ptr += 5;                                     /* increase record pointer */

  if (nw < ptr + ev->nhits * 3 + 2)       /* if there are not enough words */
    return ERRLINE;

/************* allocate memory for the hits and fill them ******************/

  if (ev->nhits > 0)
    ev->h = malloc(ev->nhits * sizeof(mhit)); 
                                              /* allocate mem for the hits */
  for (i = 0; i < ev->nhits; i++)
    rdmc_init_mhit( &(ev->h[i]));
  for (i = 0; i < ev->nhits; i++) {
    ev->h[i].ch = rdmc_mcgetf(fp);                         /* contains ch and m */
    ev->h[i].ma = ev->h[i].ch/10000;                 /* the first 2 digits */
    ev->h[i].mt = fmod(ev->h[i].ch/100,100.);        /* the digits 3 and 4 */
    ev->h[i].ch = fmod(ev->h[i].ch,100.);             /* the last 2 digits */
    ev->h[i].ch--;                          /* for "C" style (index 0,...) */
    ev->h[i].amp = rdmc_mcgetf(fp);                                /* amplitude */
    ev->h[i].t = rdmc_mcgetf(fp);                                       /* time */
    ev->h[i].tot = -1.0;                             /* no TOT information */
  } /* for i */


  rdmc_mcseek(fp, 1);                                         /* one dummy word */
  ev->gen->t = rdmc_mcgetf(fp);                  /* time corresponding gen->xyz */
  ev->gen->nmuon = ev->ntrack;                           /* number of muons */

  if (ev->gen->costh < -1.5) {            /* if costh = -2 (no gen. track) */
    free(ev->gen);                                  /* free the gen. track */
    ev->gen = NULL;                                   /* reset the pointer */
    ev->ntrack = 0;                    /* set the number of tracks to zero */
  } /* if costh < -1.5 */

  ptr += ev->nhits * 3 + 2;

  if (nw <= ptr + 2) {                             /* if no reconstruction */
    ev->nfit = 0;                       /* number of reconstr. tracks is 0 */
    rdmc_mcseek(fp, nw-ptr);                   /* read the last (unused) fields */
    if (fp->info.bai.fortran_recs)
      rdmc_mcseek(fp, 1);             /* read the field used for fortran format */
    return 0;
  } /* if no reco */

  if (nw < ptr + 10) return ERRLINE;

/*********** allocate memory for reconstructed track and fill it ***********/

  {
    mtrack rec;
    mevt_special_t fresult;
    rdmc_init_mtrack(&rec);
    rdmc_init_fit_jk(&fresult,1);

    rec.costh = rdmc_mcgetf(fp);             /* cos theta, or -99 if no reco */
    rec.phi = rdmc_mcgetf(fp)*PID180;                   /* reconstructed phi */
    rec.x = rdmc_mcgetf(fp)/100.0;                   /* reconstructed vertex */
    rec.y = rdmc_mcgetf(fp)/100.0;
    rec.z = rdmc_mcgetf(fp)/100.0;
    rec.t = rdmc_mcgetf(fp);                           /* corresponding time */
    rec.id = MUON_PLUS;                       /* we reconstructed muons */
    rec.e = 0.0;                             /* no energy reconstructed */
    rec.length = RDMC_BIG;                          /* infinity track length */
    rec.nmuon = 1;                                 /* default: one muon */
    rec.tag = 1;
    rdmc_tau_tr(&rec);                       /* calc the direction cosinuus */

    /************** read jaanus' reconstruction results ************************/
    
    fresult.val[JK_SIGTH] = rdmc_mcgetf(fp);
    fresult.val[JK_COVMAX] = rdmc_mcgetf(fp);
    fresult.val[JK_COVMIN] = fmod(fresult.val[JK_COVMAX],1000.)/1000.;
    fresult.val[JK_COVMAX] = (fresult.val[JK_COVMAX] /1000. 
			      - fresult.val[JK_COVMIN])/1000.;
    fresult.val[JK_CUTFLAG] = rdmc_mcgetf(fp);
    fresult.val[JK_FITID] = 
      fp->info.bai.rec_vers + 10000 * fp->info.bai.min_vers;
    ptr += 9;
    
    rdmc_mcseek(fp, nw-ptr);                     /* read the last (unused) fields */

    fresult.val[JK_RCHI2] = -1.0; /* in the baikal format is no rchi2  */
    fresult.val[JK_CHI2] = -1.0;    /* in the baikal format is no chi2 */
    fresult.val[JK_PROB] = -1.0;    /* in the baikal format is no prob */
    
    if (rec.costh >= -2.0)  {                     /* if cos theta is -99 */
      rdmc_add_fit(ev,&rec,&fresult,0);
    }
  } /* read reco */
  
  if (fp->info.bai.fortran_recs)
    rdmc_mcseek(fp, 1);           /* read the field used for fortran format */
  
  return 0;
  
} /* function revt_zeu() */

/****************************************************************************/
/* revt_mos() reads the next moscow format event (BAIKAL_MOS)               */
/****************************************************************************/

int rdmc_revt_baikal_mos(mcfile *fp, mevt *ev, const array *ar)
{

  int nw;                                       /* number of words to read */
  int ptr = 0;                /* pointer to the current word in the header */
  int i,j;                                                /* loop variable */
  int nint;                            /* number of interactions of a muon */
  mtrack gentrack;                       /* struct to hold the first track */

  rdmc_clear_mevt(ev);                                  /* clear old event info */

  if (fp->info.bai.fortran_recs)
    rdmc_mcseek(fp, 1);                           /* skip fortran record header */

  nw = rdmc_mcgeti(fp) - 1;                    /* number of words in the record */
  if (feof(fp->fp)) return EOF;                       /* file end reached? */

  if (nw < 9)  return ERRLINE;

/***************** read the first (common) informations ********************/

  ev->t_offset = 0.0;                     /* time offset in baikal is zero */


  ev->nhits = rdmc_mcgetf(fp);                                /* number of hits */
  ev->nrun = 0;                                 /* no run number is stored */
  ev->enr = ++(fp->info.bai.enr);         /* event number from file pointer */
  ptr += 1;

/************ read the first generated track (central muon) ****************/

  rdmc_init_mtrack(&gentrack);
  gentrack.costh = -cos(rdmc_mcgetf(fp)*PID180);                   /* cos theta */
  gentrack.phi = (180.0+rdmc_mcgetf(fp))*PID180;                  /* phi in rad */
  if (gentrack.phi > 360.0) gentrack.phi -= 360.0;       /* phi in (0,360) */
  gentrack.x = rdmc_mcgetf(fp);                                  /* read the */
  gentrack.y = rdmc_mcgetf(fp);                               /* coordinates */
  gentrack.z = rdmc_mcgetf(fp);                              /* of the track */
  gentrack.t = rdmc_mcgetf(fp);                        /* time of the vertex */
  gentrack.id = MUON_PLUS;                           /* we generated muons */
  gentrack.length = RDMC_BIG;                       /* infinity track length */
  gentrack.nmuon = 1;                      /* one muon produces this track */
  gentrack.tag=1;
  rdmc_tau_tr(&gentrack);                     /* calc the direction cosinuus */
  gentrack.e = rdmc_mcgetf(fp) * 1e3;                       /* Energy in MeV */
  ev->ntrack = (int)rdmc_mcgetf(fp) % 100;               /* number of tracks */
                                                /* higher part is nmuonall */
  rdmc_mcgetf(fp);                          /* one dummy ("EMMAX + 1000* ESUM") */
  ptr += 8;                                     /* increase record pointer */

/******* allocate memory for generated track(s) and fill them **************/

#if 0 /* i dont know if I really need the first track (axis) */
  ev->ntrack++;                 /* one additional track for the group axis */
#endif

  if (ev->ntrack > 0)
    ev->gen = malloc(ev->ntrack*sizeof(mtrack));
                                       /* allocate mem for generated track */

#if 0 /* i dont know if I really need the first track (axis) */
  memcpy(ev->gen,&gentrack,sizeof(mtrack)); /* fill the first track (axis) */

  for (i = 1; i < ev->ntrack; i++){                      /* for all tracks */
#else
  for (i = 0; i < ev->ntrack; i++){                      /* for all tracks */
#endif
    rdmc_init_mtrack(ev->gen+i);
    nint = rdmc_mcgetf(fp);               /* read the number of interactions */
    ev->gen[i].costh = gentrack.costh;                            /* theta */
    ev->gen[i].phi = gentrack.phi;                 /* and phi are the same */
    ev->gen[i].t = gentrack.t;                        /* and the time, too */
    ev->gen[i].x = rdmc_mcgetf(fp);                           /* koordinates */
    ev->gen[i].y = rdmc_mcgetf(fp);                                /* of the */
    ev->gen[i].z = rdmc_mcgetf(fp);                                /* vertex */
    rdmc_tau_tr(ev->gen+i);                   /* calc the direction cosinuus */
    ev->gen[i].nmuon = 1;                  /* one muon produces this track */
    ev->gen[i].length = RDMC_BIG;                        /*infinity length */
    ev->gen[i].id = MUON_PLUS;
    ev->gen[i].tag = i+1;
    ev->gen[i].e = 0.0;                           /* reset the muon energy */
    ptr += 4;                                              /* 4 words read */
    for (j = 0; j < nint; j++) {                   /* for all interactions */
      ev->gen[i].e += rdmc_mcgetf(fp)*1e3;        /* add the interaction energy */
      rdmc_mcseek(fp,3);                          /* skip the interaction point */
      ptr += 4;                                            /* 4 words read */
    } /* for j */

  } /* for i */
  
  ev->gen[0].e = gentrack.e;       /* the first track gets all muon energy */

  if (nw < ptr + ev->nhits * 5)           /* if there are not enough words */
    return ERRLINE;

/************* allocate memory for the hits and fill them ******************/

  if (ev->nhits > 0)
    ev->h = malloc(ev->nhits * sizeof(mhit)); 
  for (i = 0; i < ev->nhits; i++) {
    rdmc_init_mhit(&(ev->h[i]));
  } 
            /* allocate mem for the hits */
  for (i = 0; i < ev->nhits; i++) {
    ev->h[i].ch = rdmc_mcgetf(fp);                        /* channel number */
    ev->h[i].ma = ((int)rdmc_mcgetf(fp)+1000) % 1000;   /* amp-defining muon */
    ev->h[i].mt = ((int)rdmc_mcgetf(fp)+1000) % 1000;  /* time-defining muon */
    ev->h[i].ch--;                          /* for "C" style (index 0,...) */
    ev->h[i].amp = rdmc_mcgetf(fp);                            /* amplitude */
    ev->h[i].t = rdmc_mcgetf(fp);                                  /* time */
    ev->h[i].tot = -1.0;                             /* no TOT information */
  } /* for i */

  ptr += ev->nhits * 5;

  rdmc_mcseek(fp, nw-ptr);                 /* read the last (unused) fields */

  if (fp->info.bai.fortran_recs)
    rdmc_mcseek(fp, 1);         /* read the field used for fortran format */

  return 0;

} /* function revt_mos() */

/****************************************************************************/
/* revt_data() reads the next *data* event                                  */
/****************************************************************************/

int rdmc_revt_baikal_data(mcfile *fp, mevt *ev, const array *ar)
{

  int nw;                                       /* number of words to read */
  int ptr = 0;                /* pointer to the current word in the header */
  int i;                                                  /* loop variable */

  rdmc_clear_mevt(ev);                                  /* clear old event info */

  if (fp->info.bai.fortran_recs)
    rdmc_mcseek(fp, 1);                           /* skip fortran record header */

  nw = rdmc_mcgeti(fp);                        /* number of words in the record */
  if (feof(fp->fp)) return EOF;                       /* file end reached? */

  if (nw < 2)  return ERRLINE;

/***************** read the first (common) informations ********************/

  i = rdmc_mcgeti(fp);                     /* number of channels and run number */
  ev->nhits = fmod(i,1000.);                         /* number of channels */
  ev->nrun = i/1000;                                         /* run number */

  ev->enr = rdmc_mcgeti(fp);                                    /* event number */
  fp->info.bai.enr = ev->enr;
  ev->ntrack = 0;                                         /* no generation */

  ptr += 2;
  if (ptr + 3*ev->nhits > nw)/* if there are not enough words for hit info */
    return ERRLINE;

/************* allocate memory for the hits and fill them ******************/

  if (ev->nhits > 0)
    ev->h = malloc(ev->nhits * sizeof(mhit)); 
  for (i = 0; i < ev->nhits; i++) { /* init */
    rdmc_init_mhit( &(ev->h[i]));
  }

  for (i = 0; i < ev->nhits; i++) {                  /* read the hit infos */
    ev->h[i].ch = rdmc_mcgetf(fp);                        /* channel nr (1,...) */
    ev->h[i].ch--;                          /* for "C" style (index 0,...) */
    ev->h[i].amp = rdmc_mcgetf(fp);                    /* amplitude of this hit */
    ev->h[i].t = rdmc_mcgetf(fp);                           /* time of this hit */
    ev->h[i].mt = ev->h[i].ma = RDMC_NA;       /* there is no particle info */
    ev->h[i].tot = -1.0;                             /* no TOT information */
  } /* for i */


  ptr += 3*ev->nhits;
  if (nw <= ptr + 2) {                             /* if no reconstruction */
    ev->nfit = 0;                                        /* reset fit flag */
    rdmc_mcseek(fp, nw-ptr);                   /* read the last (unused) fields */
    if (fp->info.bai.fortran_recs)
      rdmc_mcseek(fp, 1);             /* read the field used for fortran format */
    return 0;
  } /* if no reco */

  if (nw < ptr + 9) return ERRLINE;

/*********** allocate memory for reconstructed track and fill it ***********/
  {
    mtrack rec;
    mevt_special_t fresult;

    rdmc_init_mtrack(&rec);
    rdmc_init_fit_jk(&fresult,1);

    rec.costh = rdmc_mcgetf(fp);             /* cos theta, or -99 if no reco */
    rec.phi = rdmc_mcgetf(fp)*PID180;                   /* reconstructed phi */
    rec.x = rdmc_mcgetf(fp)/100.0;              /* reconstructed coordinates */
    rec.y = rdmc_mcgetf(fp)/100.0;                        /* of track origin */
    rec.z = rdmc_mcgetf(fp)/100.0;
    rec.t = rdmc_mcgetf(fp);           /* reconstructed time of track origin */
    rec.id = MUON_PLUS;                       /* we reconstructed muons */
    rec.e = 0.0;                             /* no energy reconstructed */
    rec.length = RDMC_BIG;                       /* infinity track length */
    rec.nmuon = 1;                                 /* default: one muon */
    rec.tag = 1;
    rdmc_tau_tr(&rec);                    /* calc the direction cosinuus */


/************** read jaanus' reconstruction results ************************/

    fresult.val[JK_SIGTH] = rdmc_mcgetf(fp);
    fresult.val[JK_COVMAX] = rdmc_mcgetf(fp);
    fresult.val[JK_COVMIN] = fmod(fresult.val[JK_COVMAX],1000.)/1000.;
    fresult.val[JK_COVMAX] = (fresult.val[JK_COVMAX] /1000. 
				 - fresult.val[JK_COVMIN])/1000.;
    fresult.val[JK_CUTFLAG] = rdmc_mcgetf(fp);
    ptr += 9;

    rdmc_mcseek(fp, nw-ptr);            /* read the last (unused) fields */
    
    fresult.val[JK_RCHI2] = -1.0; /* in the baikal format is no rchi2  */
    fresult.val[JK_CHI2] = -1.0;    /* in the baikal format is no chi2 */
    fresult.val[JK_PROB] = -1.0;    /* in the baikal format is no prob */
    
    if (ev->rec->costh >= -2.0)  {                 /* if cos theta is -99 */
      rdmc_add_fit(ev,&rec,&fresult,0);
    } /* if no reco */
  }
  if (fp->info.bai.fortran_recs)
    rdmc_mcseek(fp, 1);       /* read the field used for fortran format */
  return 0;
  
} /* function revt_data() */

/****************************************************************************/
/* warr_mc() writes the array info to the baikal-like file                  */
/****************************************************************************/

int rdmc_warr_baikal_mc(mcfile *fp, const array *ar)
{
  int nw;
  int i,j;
  double mc_ind_muon, mc_id_spec;                       /* Jaanus' mc flags */
  double mc_vers ;

  nw = 14 + ar->nch * (5 + 2 * BAIKAL_PMTPERCHANNEL);    /* number of words */
  if (fp->info.bai.fortran_recs)           /* if it is a formatted fortran file */
    rdmc_mcputi(4 * nw + 4, fp);                      /* write the record length */

  rdmc_mcputi(nw, fp);                 /* write the number of words as first int */
  
  rdmc_mcputf(0.0, fp);                                   /* write something (?) */

  mc_vers=atof(fp->info.bai.mc_vers);
  mc_vers = ((int) mc_vers) % 10000 ;
  rdmc_mcputf(mc_vers + 10000.0 * fp->fmajor, fp);         /* format version */
  rdmc_mcputf((double)ar->nch, fp);
  rdmc_mcputf((double)ar->nstr, fp);

  for (i = 0; i < ar->nch; i++) {                       /* for all channels */
    rdmc_mcputf(100.0*ar->x[i], fp);
    rdmc_mcputf(100.0*ar->y[i], fp);
    rdmc_mcputf(100.0*ar->z[i], fp);
    rdmc_mcputf(ar->costh[i], fp);
    rdmc_mcputf((double)ar->str[i], fp);
  } /* for i */

  for (i = 0; i < ar->nch; i++) {
    for (j = 0; j < BAIKAL_PMTPERCHANNEL; j++) {
      rdmc_mcputf(ar->thresh[i], fp);               /* threshold for each pmt */
      rdmc_mcputf(ar->sensit[i], fp);              /* sensitivity of each pmt */
   } /* for j */
  } /* for i */

  rdmc_mcputf((double)rdmc_o_rdateconv(fp->info.bai.time), fp); /* write time in YYMMDD format */

#if 0
  /* igen does contain the random seed since ASCII version -6 */
  mc_ind_muon = fp->info.bai.igen /10000;
  mc_id_spec = fp->info.bai.igen % 10000;
#else
  mc_ind_muon=mc_id_spec=0.;
#endif

  rdmc_mcputf(mc_ind_muon, fp);
  rdmc_mcputf(mc_id_spec, fp);
  rdmc_mcputf((double)fp->info.bai.igtrack, fp);
  rdmc_mcputf(0.0, fp);
  rdmc_mcputf(0.0, fp);

  rdmc_mcputf((double)fp->info.bai.nw_rec_arr, fp);
  rdmc_mcputf((double)fp->info.bai.nw_rec_evt, fp);
  rdmc_mcputf((double)fp->info.bai.rec_vers, fp);
  rdmc_mcputf((double)fp->info.bai.min_vers, fp);

  if (fp->info.bai.fortran_recs)                     /* if it is a formatted fortran file */
    rdmc_mcputi(4 * nw + 4, fp);                      /* write the record length */

#if (DEBUG == 1)
    fprintf(stderr, "<Write><Header: nch=%i nstr=%i>\n",ar->nch,ar->nstr);
#endif

  return errno;

} /* function warr_mc() */

/****************************************************************************/
/* function wevt_mc() writes an event to  file                              */
/****************************************************************************/

int rdmc_wevt_baikal_mc(mcfile *fp, const mevt *ev, const array *ar)
{
  int nw;                                                /* number of words */
  int i;                                                   /* loop variable */
  
  nw = 20 + 3*ev->nhits;

  if (fp->info.bai.fortran_recs)            /* record length for formatted files */
    rdmc_mcputi(4*nw+4,fp);
  
/***************** write the first (common) informations ********************/

  rdmc_mcputi(nw,fp);                           /* first word is number of words */
  if (fp->fmajor == BAIKAL_NEW) {                     /* for the new format */
    rdmc_mcputi(ev->nhits+1000*ev->nrun,fp);
    rdmc_mcputi(ev->enr,fp);
  } /* if new format */

/******************* write the generating track *****************************/

  if (ev->ntrack > 0) {                   /* if there is a generating track */
    if (ev->gen[0].costh > 0)
      rdmc_mcputf(ev->nhits*10.0+ev->gen[0].costh, fp);
    else
      rdmc_mcputf(-ev->nhits*10.0+ev->gen[0].costh, fp);
    rdmc_mcputf(ev->gen[0].phi/PID180, fp);
    rdmc_mcputf(ev->gen[0].x*100.0, fp);
    rdmc_mcputf(ev->gen[0].y*100.0, fp);
    rdmc_mcputf(ev->gen[0].z*100.0, fp);
    rdmc_mcputf(ev->gen[0].e/1000.0, fp);               /* baikal: energy in GeV */
  } /* if gen */
  else {                                          /* if no generating track */
    rdmc_mcputf(-ev->nhits*10.0-2.0, fp);            /* seems to be not good ole */
/*    mcputf(ev->nhits*10.0, fp); */        /* this could be an alternative */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
  } /* if not gen */

/***************************** write the hits *******************************/

  for (i = 0; i < ev->nhits; i++) {                           /* store hits */
    rdmc_mcputf(ev->h[i].ch+1.0+ev->h[i].mt*100.0+ev->h[i].ma*10000.0, fp);
    rdmc_mcputf((double)ev->h[i].amp, fp);
    rdmc_mcputf((double)ev->h[i].t, fp);
  }
  
  rdmc_mcputf(1.0, fp);                                /* I dont know what it is */

  if (ev->ntrack > 0)                     /* if there is a generating track */
    rdmc_mcputf(ev->gen[0].t, fp);         /* time for origin of generated track */
  else
    rdmc_mcputf(0.0, fp);

/******************* write the reconstructed track **************************/
  
  if (ev->nfit > 0) {                      /* if there is a reconstrunction */
    rdmc_mcputf(ev->rec[0].costh, fp);
    rdmc_mcputf(ev->rec[0].phi/PID180, fp);
    rdmc_mcputf(ev->rec[0].x*100.0, fp);
    rdmc_mcputf(ev->rec[0].y*100.0, fp);
    rdmc_mcputf(ev->rec[0].z*100.0, fp);
    rdmc_mcputf(ev->rec[0].t, fp);
    if ( (ev->fresult != NULL)              /* if there is an jk record */
	 && (0 <= ev->fresult->id )
	 && (ar->n_fit > ev->fresult->id )  /* there is a fit defined  */
	 && (rdmc_is_this_jk(&(ar->def_fit[ev->fresult->id])
				,ev->fresult) ) )    {
      rdmc_mcputf(ev->fresult->val[JK_SIGTH], fp);
      rdmc_mcputf(ev->fresult->val[JK_COVMIN]*1e3 
	     + 1e6*ev->fresult->val[JK_COVMAX], fp);
      rdmc_mcputf(ev->fresult->val[JK_CUTFLAG], fp);
    } /* if res */
    else {                                      /* if there is no jk record */
      rdmc_mcputf(0.0, fp);                                   /* set all to zero */
      rdmc_mcputf(0.0, fp);                                   /* set all to zero */
      rdmc_mcputf(0.0, fp);                                   /* set all to zero */
    } /* if not res */
  } /* if reconstruction */
  else {                                   /* if there is no reconstruction */
    rdmc_mcputf(-99., fp);
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
    rdmc_mcputf(0.0, fp);                                     /* set all to zero */
  } /* if no reconstruction */

  if (fp->fmajor == BAIKAL_NEW)
    rdmc_mcputi(0, fp);                           /* flag for further subrecords */

  if (fp->info.bai.fortran_recs)            /* record length for formatted files */
    rdmc_mcputi(4*nw+4,fp);
  
  return 0;

} /* function wevt() */


/****************************************************************************/
/* function wevt_data() writes an event to file                             */
/****************************************************************************/

int rdmc_wevt_baikal_data(mcfile *fp, const mevt *ev, const array *ar)
{

  return rdmc_wevt_baikal_mc(fp,ev,ar);         /* we use mc format instead of data format */

} /* function wevt_data() */

/****************************************************************************/
/* rhd_mc() reads the format relevant informations for baikal-like formats, */
/*          and then rewinds                                                */
/****************************************************************************/

int rdmc_rhd_baikal_mc(mcfile *fp) 
{
  long i1,i2;
  float f1,f2,f3,f4;

  i1 = rdmc_mcgeti(fp);
  i2 = rdmc_mcgeti(fp);
  
  if ((i1 == i2*4+4) || (rdmc_swap(i1) == rdmc_swap(i2) * 4 + 4)){/* if the first int */
    fp->info.bai.fortran_recs = 1;  /* contains the record length -> fortran file */
  }
  rewind(fp->fp);                               /* rewind file to the begin */
  if (fp->info.bai.fortran_recs) rdmc_mcseek(fp,1);  /* skip the fortran record header */
  i1 = rdmc_mcgeti(fp);                /* contains the number of words in header */
  if (i1 == 2) {                             /* if there are only two words */
    fp->info.bai.mc = 0;              /* it is a non-reconstrcted data file */
    fp->info.bai.nrun = rdmc_mcgetf(fp);                  /* first word is the run number */
    fp->info.bai.time = 0;                                /* no time stored */
    fp->info.bai.rec_vers = 0;                         /* no reconstruction */
    fp->info.bai.min_vers = 0;                           /* no minimization */
    rewind(fp->fp);                                      /* rewind the file */
    return 0;                                                 /* and return */
  } /* if i1 = 2 */
  f1 = rdmc_mcgetf(fp);                                             /* dummy a(1)*/
  f2 = rdmc_mcgetf(fp);                        /* contains mc and format version */
  f3 = rdmc_mcgetf(fp);                                    /* number of channels */
  f4 = rdmc_mcgetf(fp);                                 /* the number of strings */

  if ((f4 != (int)f4)                                /* should be a integer */
    ||(f4 <= 0.0)                                     /* and greater than 0 */
    ||(f4 > 20.0)) {    /* and less than 20 (we havn't such big telescopes) */
      fp->info.bai.swap = 1;              /* perhaps the bytes are swapped */
    }                                         /* try it with swapped bytes! */
  rewind(fp->fp);                               /* rewind file to the begin */
  if (fp->info.bai.fortran_recs) rdmc_mcseek(fp,1); /* skip the fortran record header */
  i1 = rdmc_mcgeti(fp);                /* contains the number of words in header */
  if (i1 <= 2) {                             /* if there are too less words */
    fp->info.bai.mc = 0;              /* it is a non-reconstrcted data file */
    fp->info.bai.nrun = rdmc_mcgetf(fp);       /* first word is the run number */
    fp->info.bai.time = 0;                               /* no time stored */
    fp->info.bai.rec_vers = 0;                         /* no reconstruction */
    fp->info.bai.min_vers = 0;                           /* no minimization */
    rewind(fp->fp);                                      /* rewind the file */
    return 0;                                                 /* and return */
  } /* if i1 = 2 */
  f1 = rdmc_mcgetf(fp);                                            /* dummy a(1) */
  f2 = rdmc_mcgetf(fp);                        /* contains mc and format version */
  f3 = rdmc_mcgetf(fp);                                    /* number of channels */
  f4 = rdmc_mcgetf(fp);                                 /* the number of strings */

  fp->fmajor = (int) f2 / 10000;           /* format version in second word */
  fp->fminor = 0;
  sprintf(fp->info.bai.mc_vers,"%i", ((int) f2) % 10000);

  if ((f1 >= 0.9) && (f1 <= 10.0)&& (fp->fmajor == BAIKAL_OLD))
    fp->fmajor = BAIKAL_MOS;              /* 1.0 as 1st word: moscow format */

  rdmc_mcseek(fp, (int)f3*9);                               /* skip the pmt data */
  fp->info.bai.time = rdmc_o_dateconv(rdmc_mcgetf(fp));            /* time in sec since 1.1.70 */
  fp->info.bai.igen = 10000 * rdmc_mcgetf(fp);               /* k+2 is mc model */
  fp->info.bai.igen += rdmc_mcgetf(fp);                    /* k+3 is mc spectrum */
  if (fp->info.bai.igen >= 10000 )
    fp->info.bai.igen = 103;
  else
    fp->info.bai.igen = 0;
#if 0
  fp->igtrack = rdmc_mcgetf(fp);                            /* k+4 is mc trigger */
#else
  fp->info.bai.nrun = rdmc_mcgetf(fp); /*ignore them */
#endif
  rdmc_mcseek(fp, 2);                                             /* two dummies */
#if 0
  fp->info.bai.nrun = fp->igtrack;                      /* it seems to be so */
#endif

  i1 -= 10 + f3 * 9;              /* decrease the header length s a pointer */

  if (i1-- > 0) fp->info.bai.nw_rec_arr = rdmc_mcgetf(fp);
  if (i1-- > 0) fp->info.bai.nw_rec_evt = rdmc_mcgetf(fp);
  if (i1-- > 0) fp->info.bai.rec_vers = rdmc_mcgetf(fp);
  if (i1-- > 0) fp->info.bai.min_vers = rdmc_mcgetf(fp);

  rewind(fp->fp);                            /* rewind the file so that one */
                                       /* starts reading from the beginning */
  fp->info.bai.mc_id = (fp->info.bai.mc) ? 'U':'D';

  return 0;

} /* function rhd_mc() */


/****************************************************************************/
/* mcgeti() reads an 4 byte int from a mc/data file                         */
/****************************************************************************/

long rdmc_mcgeti(mcfile *fp)
{

  long r;
  fread(&r, sizeof(r), 1, fp->fp);
  return (fp->info.bai.swap)? rdmc_swap(r):r;

} /* function rdmc_mcgeti() */


/****************************************************************************/
/* mcputi() writes an 4 bytes int to a mc/data file                         */
/****************************************************************************/

long rdmc_mcputi(long i, mcfile *fp)
{
  long r;

  r = (fp->info.bai.swap)? rdmc_swap(i) : i;
  fwrite(&r, sizeof(r), 1, fp->fp);
  return r;

} /* function mcputi() */


/****************************************************************************/
/* mcgetf() reads a 4 byte float form a mc/data file                        */
/* attention: altough mcgetf() reads a float it returns double according to */
/* the C standard                                                           */
/****************************************************************************/

double rdmc_mcgetf(mcfile *fp)
{

  union {
    long r;
    float f;
  } u;

  fread(&u.r, sizeof(u.r), 1, fp->fp);

  if (fp->info.bai.swap) u.r = rdmc_swap(u.r);
  return (double)u.f;

} /* function mcgetf() */


/****************************************************************************/
/* mcputf() writes a 4 byte float to a mc/data file                         */
/* attention: altought mcputf writes a 4 byte float, it requires a double   */
/* argument because this is standard in C                                   */
/****************************************************************************/

double rdmc_mcputf(double d, mcfile *fp)
{

  union {
    long r;
    float f;
  } u;


  u.f = d;
  if (fp->info.bai.swap) u.r = rdmc_swap(u.r);
  fwrite(&u.r, sizeof(u.r), 1, fp->fp);

  return d;

} /* function mcgetf() */


/****************************************************************************/
/* mcseek() seeks over a number of 4-byte-words (int, float)                */
/****************************************************************************/

int rdmc_mcseek(mcfile *fp, long n)
{

  if (feof(fp->fp)) return EOF;

  return fseek(fp->fp, n * sizeof(long), SEEK_CUR);

} /* function mcseek() */

/****************************************************************************/
/* skipevt_baikal_mc() skips to the next record (event) in a baikal mc/data file   */
/****************************************************************************/

int rdmc_skipevt_baikal_mc(mcfile *fp)
{

  int nw;                          /* number of words in the current record */
  int r;                                                    /* return value */

  if (fp->info.bai.fortran_recs)
    rdmc_mcgeti(fp);

  nw = rdmc_mcgeti(fp);
  if (nw < 0) return ERRLINE;

  r = rdmc_mcseek(fp,nw);

  fp->info.bai.enr++;                              /* increase event number */

  if (fp->info.bai.fortran_recs)
    rdmc_mcgeti(fp);

  return r;

} /* function skipevt_mc() */

/****************************************************************************/
/* swap() swaps the high and low bytes of a 4 byte int                      */
/****************************************************************************/

long rdmc_swap(long i)
{

  union {
    unsigned char c[sizeof(long)];
    long l;
  } u;                                          /* help unions for swapping */
  unsigned char r;

  int n;

  u.l = i;

  for (n = 0; n < sizeof(long)/2; n++) { /* for all bytes in the first half */
    r = u.c[n];                                      /* swap the first byte */
    u.c[n] = u.c[sizeof(long) - n - 1];                    /* with the last */
    u.c[sizeof(long) - n - 1] = r;
  } /* for n */

  return u.l;

} /* function swap() */

static void  rdmc_add_fitdef_baikal(array *ar){

  array_hdef_t def_fit;

  rdmc_jk_to_fitdef(&def_fit, 1);
  rdmc_add_fit_def(ar,&def_fit, 0);
}


#endif /* BAIKAL_BIN_F */
/****************************************************************************/
/********************************** E O F ***********************************/
/****************************************************************************/
/* 
   This is just for EMACS:
   Local Variables:
   compile-command: "cd .. ; make -k rdmc" 
   End: 
*/
