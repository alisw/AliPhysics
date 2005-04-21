/*
 * Read/write files in the SiEGMuND native format
 */

char *dumand_rdmc_cvsid = 
"$Header: /net/local/cvsroot/siegmund/rdmc/dumand.c,v 1.29 2004/04/15 14:45:51 wiebusch Exp $";

#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

////////////////#include <unistd.h>


#include <fcntl.h>
#include <sys/stat.h>

#if defined(OSF1) || defined(AIX) || defined (SunOS)
#include <errno.h>
#endif

#if defined(OSF1) || defined (SunOS) || defined (IRIX)
#include <alloca.h>
#endif

#include "rdmc.h"

#ifdef DUMAND_ASCII_F

#include "rdmc_local.h"

#include "dumand.h"

#if RDMC_MAXTOKEN_PER_LINE < 100
#define DUMAND_USER_MAX RDMC_MAXTOKEN_PER_LINE
#else
#define DUMAND_USER_MAX 100
#endif

#define DUMAND_MAXUSER_PER_LINE 15


#ifndef DEBUG
#define DEBUG 0                       /* DEBUG 1 writes some info to stderr */
#endif

#define ERRLINE (-__LINE__)

#define PID180 (M_PI/180.)

static int rdmc_a_M(const char *s, mcfile *fp);                /* read the several */
static int rdmc_a_G(const char *s, array *ar);               /* from a dumand-like */
static int rdmc_a_P(const char *s, array *ar);                           /* format */
static int rdmc_a_K(const char *s, array *ar);
static int rdmc_a_Q(const char *s, array *ar);
static int rdmc_a_E(const char *s, mevt *ev);
static int rdmc_a_H(const char *s, mhit *h); 
static int rdmc_a_T(const char *s, mtrack *t);
static int rdmc_a_F(const char *s, mtrack *t, mevt_special_t *fit);   
static int rdmc_a_U(const char *s, mevt_special_t **user, int *nuser);
static int rdmc_a_C(const char *s, char **comment);
#if 0
static int a_skip_C(mcfile *fp);                     /* skip comment lines */
#endif
static int parse_fit_def_ascii(array *ar, const char *s);   
/* parse the fit definitions */
static int  declare_rdmc_dumand_user(array *ar); 


/****************************************************************************/
/* The function rarr_ascii() reads the header of a dumand like file         */
/****************************************************************************/

int rdmc_rarr_ascii(mcfile *fp, array *ar)
{
  char s[RDMC_MAXLINE];                                /* input line */
  int c;                                        /* first char of input line */
  int ich = 0;                                              /* loop variable */
  int *iclust;               /* pointer to current channel number in string */
  int retval;                                                /* returnvalue */

  rdmc_init_array(ar);                                      /* reset the array */
  ar->nrun=fp->info.dum.nrun;
  for (*s=c=getc(fp->fp); strchr("CGQ", *s) && (c != EOF); *s=c=getc(fp->fp)) {
    if (fgets(s+1, RDMC_MAXLINE-1, fp->fp) == NULL)
      return EOF;
    fp->fpos++;

    switch (*s) {
    case 'C':
      if ((retval = rdmc_a_C(s,&(ar->comment))))
	return retval;
      break;
    case 'Q':
      if ((retval = rdmc_a_Q(s, ar)))                 /* read trigger info */
	return retval;
      break;
    case 'G':
      if ((retval = rdmc_a_G(s, ar)))                 /* read array info */
	return retval;
      break;
    } /* switch *s */
  } /* while noe P nor V or E */

  ungetc(c, fp->fp); /* ungetch the last char */
  
  if (ar->nstr == 0) ar->nstr = 1;                    /* minimal one string */

  iclust = calloc(ar->nstr,sizeof(int));  /* alloc mem for channel nr array */

  for (*s = c = getc(fp->fp); (*s != 'E') && (c != EOF); *s=c=getc(fp->fp)) {
                                  /* for all lines *not* beginning with 'E' */
    ungetc(c, fp->fp);
    do {                                     /* read optional comment lines */
      *s = c = getc(fp->fp);
      if (c == EOF)
	return EOF;
      if (toupper((int) *s) != 'C')
	break;

      if (fgets(s+1, RDMC_MAXLINE-1, fp->fp) == NULL)
	return EOF;
      fp->fpos++;
    } while (rdmc_a_C(s,&(ar->comment)) == 0);

    if (*s == 'E') break;
    ungetc(c, fp->fp);

    fp->fpos += 1;
    if (fgets(s,RDMC_MAXLINE,fp->fp) == NULL)              /* try to read a line */
      return EOF;
  
    if ((retval = rdmc_a_P(s,ar)) != 0)                /* read pmt position */
      return retval;

    if ((retval = rdmc_a_K(s,ar)) != 0)                    /* read pmt calib */
      return retval;

    if ((retval = rdmc_a_Q(s,ar)) != 0)                 /* read trigger pars */
      return retval;

#if 0   /* does not work */
    if ((retval = rdmc_a_C(s,&(ar->comment))) != 0)                 /* read trigger pars */
      return retval;
#endif

  } /* while !e */

  if (c != EOF) 
    ungetc(*s, fp->fp);

  for (ich = 0; ich < ar->nch; ich++) {

    if ((ar->str[ich] < 0) && (ar->str[ich] <= ar->nstr)){
                                                       /* str nr is in range */

      if (ar->clust[ich] == -1)      /* if there was no valid channel number */
	ar->clust[ich] = iclust[ar->str[ich] - 1];           /* generate one */

      iclust[ar->str[ich] - 1] = ar->clust[ich] + 1;  /* increase channel nr */
    } /* if str <= nstr */

  } /* for ich */

  free (iclust);                        /* free the mem for the ch nr array */

  /* this is a dirty patch */
  /* at the time of reading a header in dumand format we do not know */
  /* the number of fits ! */
  /* The (only ?) possible solution is to parse the history lines */
  parse_fit_def_ascii(ar,fp->creator);   /* parse the fit definitions */

  /* this durty patch just declares a default USER Tag rdmc when reading
     dumand format */
  declare_rdmc_dumand_user(ar); 

  /* now patch the fp->M creation time into ar->tbegin */
  ar->tbegin=fp->info.dum.time;

  return 0;                                           /* successfull return */

} /* function rarr_ascii() */

/****************************************************************************/
/* revt_ascii() reads the next dumand-like event                            */
/****************************************************************************/

int rdmc_revt_ascii(mcfile *fp, mevt *ev, const array *ar)
{

  char s[RDMC_MAXLINE];                                            /* input line */
  int c;                                                      /* first char */
  int ihit = 0, itrack = 0;                            /* hit/track counter */
  int i;

  rdmc_clear_mevt(ev);                                   /* clear old event info */

  if (feof(fp->fp))                          /* if file end already reached */
    return EOF;

  ev->nrun = ar->nrun;                       /* copy the run number from fp */

  do {                                            /* read optional comments */
    *s = c = getc(fp->fp);
    if (c == EOF)
      return EOF;
    if (toupper((int) *s) != 'C') {
      ungetc(c, fp->fp);
      break;
    }
    if (fgets(s+1, RDMC_MAXLINE-1, fp->fp) == NULL)
      return EOF;
    fp->fpos++;
  } while (rdmc_a_C(s,&(ev->comment)) == 0);

  *s = c = getc(fp->fp);                                     /* read a char */
  if (c == EOF)
    return EOF;
  while (*s != 'E') {                      /* it should be an 'E' ("event") */
    if (fgets(s+1, RDMC_MAXLINE-1, fp->fp) == NULL)
      return EOF;
    fp->fpos++;
    *s = c = getc(fp->fp);
    if (c == EOF)
      return EOF;
  } /* while != 'E' */

  fp->fpos += 1;
  if (fgets(s+1,RDMC_MAXLINE-1,fp->fp) == NULL)        /* try to read first line */
    return EOF;

  if (rdmc_a_E(s, ev) != 0)                      /* scan for the event info */
    return ERRLINE;
    
  if (ev->ntrack > 0){
    ev->gen = malloc(ev->ntrack*sizeof(mtrack));
                                       /* allocate mem for generated tracks */
    for (i=0 ; i < ev->ntrack ; i++)
      rdmc_init_mtrack(ev->gen+i);
  }

  if (ev->nhits > 0){                          /* allocate mem for the hits */
    ev->h = malloc(ev->nhits * sizeof(mhit)); 
    for (i=0 ; i < ev->nhits ; i++)
      rdmc_init_mhit(ev->h + i);
  }

  for (*s = c = getc(fp->fp); (*s != 'E') && (c != EOF); *s=c=getc(fp->fp)) {
                                  /* for all lines *not* beginning with 'E' */
    fp->fpos += 1;
    if (fgets(s+1,RDMC_MAXLINE-1,fp->fp) == NULL) {        /* try to read a line */
      c = EOF;
      break;
    }

    switch (s[0]) {                /* switch for the first char of the line */
    case 'F':                                                   /* fit line */
      {
	mtrack fit;
	mevt_special_t res;

	rdmc_init_mtrack(&fit);
	rdmc_init_fit_jk(&res,RDMC_NA);
	if (rdmc_a_F(s, &fit, &res) != 0)
	  return ERRLINE;
	fit.tag=ev->nfit+1; /* nfit is still 0..nfit -1 but we want 1..nfit*/
	res.id=ev->nfit;   /* reference to fit_def so 0..nfit-1 is OK */
	if (rdmc_add_fit(ev,&fit,&res,ev->nfit))
	  return ERRLINE;
      }
      break;
    case 'H':                                                        /* hit */
      if (++ihit > ev->nhits)         /* dont fill more hits than allocated */
	return ERRLINE;
      if (rdmc_a_H(s, &(ev->h[ihit-1])) != 0)              /* scan this hit */
	return ERRLINE;
      ev->h[ihit-1].id=ihit;
      break;
    case 'T':                                           /* generating track */
      if (++itrack > ev->ntrack)    /* dont fill more tracks than allocated */
	return ERRLINE;
      if (rdmc_a_T(s, &(ev->gen[itrack-1])) != 0)        /* scan this track */
	return ERRLINE;
      break;
    case 'U':                                                 /* user-def'd */
      if (rdmc_a_U(s, &(ev->user), &(ev->nuser)) != 0) 
                                    /* look for rdmc-compatible user fields */
	return ERRLINE;
      break;                                                  /* is ignored */
    case 'C':
      rdmc_a_C(s, &(ev->comment));
      break;
    case 'V':                                               /* header flags */
    case 'M':                                                        /* are */
    case 'G':                                                        /* NOT */
    case 'P':                                                    /* allowed */
      return ERRLINE;
    default:
      break;
    } /* switch s[0] */

  } /* for all lines not beginning with 'E' */

  if (ev->nhits > ihit) ev->nhits = ihit;    /* set the real number of hits */
  if (ev->ntrack > itrack) ev->ntrack = itrack;    /* set real nr of tracks */

  if (c != EOF) 
    ungetc(*s, fp->fp);

  ev->nch=rdmc_count_nch(ev);         /* calc the number of hitted channels */
  ev->nstr=rdmc_count_nstr(ev);      /* calc the number of hitted channels */

  return 0;
  
} /* function revt_ascii() */

/****************************************************************************/
/* function warr_ascii() writes the array info to a dumand-like file        */
/* This function writes out the head of a DUMAND ascii file                 */
/* opposite to reading the input file it writes not only                    */
/* the Geometry banks ('G', 'P') , but also the ('V' and 'M' flags)         */
/* so the function whd_ascii does not exist                                 */
/****************************************************************************/

int rdmc_warr_ascii(const mcfile *fp,const array *geo)
{
  int iom;                                              /* om index in loop */
  int itrigger;                                    /* trigger index in loop */

  fprintf(fp->fp, "V %i\n",                             /* write the V flag */
	   DUMAND_ASCII_VERSION);                                /* version number */

  rdmc_wrcomment_ascii(fp,fp->creator);
  rdmc_wrcomment_ascii(fp,fp->comment);
  /* write M-line */
  fprintf(fp->fp,"M %c %s %i %i %i %i %i %i\n"
	  ,fp->info.dum.mc_id
	  ,fp->info.dum.mc_vers
	  ,geo->nrun
	  ,geo->id
	  ,fp->info.dum.igen 
	  ,fp->info.dum.igtrack
	  ,fp->info.dum.daswrun
	  ,rdmc_o_rdateconv(fp->info.dum.time)       /* write the time in YYMMDD format */
	  );


  for (itrigger = 0; itrigger < geo->n_trigger; itrigger++) {
    fprintf(fp->fp, "Q U %s -1 -1 -1. -1.\n",
	    geo->def_trig[itrigger].tag);
  }

  if (geo->comment != NULL)                 /* if there is a common comment */
    rdmc_wrcomment_ascii(fp, geo->comment);

  if ((geo->nch>0)||(geo->id > 0)) {     /* if there was a detector geometry */
    fprintf(fp->fp,"G %i %i %i %.2f %.2f %.1f\n",
	    geo->id,geo->nch,geo->nstr,
	    geo->lattitude, geo->longitude, geo->depth);         /* 'G' flag */
  }
/************** write now the positions ('P') of all channels ***************/

  for (iom = 0 ; iom < geo->nch ; iom++) {

    int ori;    /* change orientation from cos theta to 1 for up and 2 down */

    ori = 0;
    if (geo->costh[iom] > 0.5) ori = 1;  /* up */
    if (geo->costh[iom] < -0.5) ori = 2; /* down */

    fprintf(fp->fp,"P %i %i %i %i %li %li %li %i %i %f %f\n",
	    iom+1,geo->type[iom],
	    geo->serial[iom],
	    geo->str[iom],
	    rdmc_nint(geo->x[iom] * 1000.),
	    rdmc_nint(geo->y[iom] * 1000.),
	    rdmc_nint(geo->z[iom] * 1000.),
	    ori,
	    geo->clust[iom] + 1,
	    geo->thresh[iom],
	    geo->sensit[iom]);
  } /* for iom */

/**** write now the calibration ('K') of all channels, if there is one *******/

  for (iom = 0 ; iom < geo->nch ; iom++) {
    if ( (geo->is_calib.adc)
	 || (geo->is_calib.tdc) ){
      fprintf(fp->fp,"K %i %.5f %.3f %.3f %.2f %.5f %.2f\n",
	      iom+1,
	      geo->cal[iom].beta_t,
	      geo->cal[iom].t_0,
	      geo->cal[iom].alpha_t,
	      geo->cal[iom].ped,
	      geo->cal[iom].beta_a,
	      geo->cal[iom].kappa);
    }
  } /* for iom */

  return 0;

} /* warr_ascii() */

/****************************************************************************/
/* function wevt_ascii() writes an event to a dumand-like file              */
/****************************************************************************/

int rdmc_wevt_ascii(const mcfile *fp,const mevt *event, const array *ar)
{
  int itra,ifit,ihit;                           /* MC-track looop varialble */
  int iuser;                                 /* loop over user-def'd fields */
                                               
  fprintf(fp->fp,"E %i %i %i %i %i %i %i %li %li %i %s\n", /* the 'E' flag */
	  event->enr,
	  event->nhits,
	  (int)rdmc_mjd_to_unixtime(event->mjd, event->secs),/* the time in secs */
	  event->nsecs,                        /* the nsec part of the time */
	  0,0,0,                             /* dummies for the Dumand time */
	  event->trigger,
	  rdmc_nint(event->t_offset),event->ntrack,"Q1");

  for (itra=0 ; itra<event->ntrack ; itra++)  {  /* write out the 'T' flags */
    fprintf(fp->fp,"T %i %li %li %li %g %g %g %li %.2f %i %li\n"
	    ,event->gen[itra].tag
	    ,rdmc_nint(event->gen[itra].x * 1000 )
	    ,rdmc_nint(event->gen[itra].y * 1000 )
	    ,rdmc_nint(event->gen[itra].z * 1000 )
	    ,event->gen[itra].px
	    ,event->gen[itra].py
	    ,event->gen[itra].pz
	    ,rdmc_nint( event->gen[itra].e )
	    ,event->gen[itra].t
	    ,event->gen[itra].id
	    ,rdmc_nint(event->gen[itra].length * 1000));
  }

  for (ihit=0 ; ihit<event->nhits ; ihit++)   {  /* write now the 'H' flags */
    fprintf(fp->fp,"H %i %i %g %.2f %.2f %s %i"
	    ,event->h[ihit].str
	    ,event->h[ihit].ch+1
	    ,event->h[ihit].tot
	    ,event->h[ihit].amp
	    ,event->h[ihit].t
	    ,"U0"
	    ,event->h[ihit].mt);
    if (event->h[ihit].ma > 0)
      fprintf(fp->fp, " %i\n",event->h[ihit].ma);
    else
      fprintf(fp->fp, "\n");
  } /* for ihit */

  for (ifit=0 ; ifit<event->nfit ; ifit++)   {    /*write now the 'F' flags */

    if ( (event->fresult != NULL)              /* if there is an jk record */
	&& (0 <= event->fresult[ifit].id )
	&& (ar->n_fit > event->fresult[ifit].id )  /* there is a fit defined  */
	&& (rdmc_is_this_jk(&(ar->def_fit[event->fresult[ifit].id])
			       ,&(event->fresult[ifit])) ) ) {

      fprintf(fp->fp
	      ,"F %i %li %li %li %g %g %g %li %.2f %g %g %g %g %g %i %i %li\n"
	      ,(int) rdmc_nint(event->fresult[ifit].val[JK_FITID])
	      ,rdmc_nint(event->rec[ifit].x * 1000. )
	      ,rdmc_nint(event->rec[ifit].y * 1000. )
	      ,rdmc_nint(event->rec[ifit].z * 1000. )
	      ,event->rec[ifit].px
	      ,event->rec[ifit].py
	      ,event->rec[ifit].pz
	      ,rdmc_nint( event->rec[ifit].e )
	      ,event->rec[ifit].t
	      ,event->fresult[ifit].val[JK_RCHI2]
	      ,event->fresult[ifit].val[JK_PROB]
	      ,event->fresult[ifit].val[JK_SIGTH]
	      ,event->fresult[ifit].val[JK_COVMIN]
	      ,event->fresult[ifit].val[JK_COVMAX]
	      ,(int) rdmc_nint(event->fresult[ifit].val[JK_CUTFLAG])
	      ,event->rec[ifit].id
	      ,rdmc_nint(event->rec[ifit].length*1000.));
    }else{
      fprintf(fp->fp
	      ,"F %i %li %li %li %g %g %g %li %.2f %g %g %g %g %g %i %i %li\n"
	      ,RDMC_NA
	      ,rdmc_nint(event->rec[ifit].x * 1000. )
	      ,rdmc_nint(event->rec[ifit].y * 1000. )
	      ,rdmc_nint(event->rec[ifit].z * 1000. )
	      ,event->rec[ifit].px
	      ,event->rec[ifit].py
	      ,event->rec[ifit].pz
	      ,rdmc_nint( event->rec[ifit].e )
	      ,event->rec[ifit].t
	      ,(float) RDMC_NA
	      ,(float) RDMC_NA
	      ,(float) RDMC_NA
	      ,(float) RDMC_NA
	      ,(float) RDMC_NA
	      , RDMC_NA
	      ,event->rec[ifit].id
	      ,rdmc_nint(event->rec[ifit].length*1000.));
    }
  }


  if (event->nuser > 0) {                        /* if there is a user field */
    int ival;
    for (iuser = 0; iuser < event->nuser; iuser++) {  /* for all user fields */
      fprintf(fp->fp, "U");                          /* write the field id */
      for(ival=0 ; ival < event->user[iuser].nval ; ival++){
	fprintf(fp->fp," %g",event->user[iuser].val[ival]); /* write the field */
      }
      fprintf(fp->fp, "\n");                                   /* end the line */
    } /* for iuser */
  } /* if nuser > 0 */

  if (event->comment != NULL)                /* if there is an event comment */
    rdmc_wrcomment_ascii(fp, event->comment);

  return 0;

} /* function wevt_ascii() */

/****************************************************************************/
/* wrcomment_ascii() writes a comment line to an ASCII file                 */
/****************************************************************************/

int rdmc_wrcomment_ascii(const mcfile *fp, const char *s)
{
  static char *line=NULL;                  /* temporary copy of the string */
  char *eol;                                         /* end of current line */
  char *bol;                                       /* beginning of the line */


  if(s){
#ifndef CRAY
    line = alloca((strlen(s)+1)*sizeof(char));
#else
    line = malloc((strlen(s)+1)*sizeof(char));
#endif
    strcpy(line,s);                                       /* copy the string */
    for (bol = eol = line; eol; bol = eol+1) {              /* for each line */
      eol = strchr(bol,'\n');                             /* search line end */
      if (eol) *eol = '\0';                          /* replace it with '\0' */
      if (*bol != '\0')                          /* if it contains something */
	fprintf(fp->fp, "C %s\n", bol+1);                  /* print the line */
    }
  }
#ifdef CRAY
  free(line);   /* Jacobsen added semicolon here */
#endif

  return 0;

} /* wrcomment_ascii() */

/****************************************************************************/
/* rhd_ascii() reads the format relevant informations for dumand-like       */
/*          formats                                                         */
/****************************************************************************/

int rdmc_rhd_ascii(mcfile *fp) 
{
  char s[RDMC_MAXLINE];                                            /* input line */
  int c;                                                      /* first char */

  do {                                              /* read comment lines */
    *s = c = getc(fp->fp);
    if (c == EOF)
      return EOF;
    if (toupper((int) *s) != 'C') {
      ungetc(c, fp->fp);
      break;
    }
    if (fgets(s+1, RDMC_MAXLINE-1, fp->fp) == NULL)
      return EOF;
    fp->fpos++;
  } while (rdmc_a_C(s,&(fp->creator)) == 0);

  *s = c = getc(fp->fp);
  if (c == EOF)
    return EOF;
  if (*s != 'M') {                          /* if there is *not* a 'M' line */
    ungetc(*s, fp->fp);                            /* ungetch the last char */
    return 0;
  } /* if *s != M */
  
  
  if (fgets(s+1,RDMC_MAXLINE-1,fp->fp) == NULL)         /* try to read next line */
    return EOF;
  fp->fpos += 1;

  if (rdmc_a_M(s,fp) != 0)                             /* monte carlo infos */
    return ERRLINE;
  
  return 0;

} /* function rhd_ascii() */

/****************************************************************************/
/* skipevt_ascii() skips the next event from a dumand-like file             */
/****************************************************************************/

int rdmc_skipevt_ascii(mcfile *fp)
{
  char s[RDMC_MAXLINE];                                            /* input line */
  int c;                                                      /* first char */

  *s = c = getc(fp->fp);                                     /* read a char */
  if (c == EOF)
    return EOF;

  fp->fpos += 1;
  if (fgets(s+1,RDMC_MAXLINE-1,fp->fp) == NULL)        /* try to read first line */
    return EOF;

  for (*s =c = getc(fp->fp); (*s != 'E') && (c != EOF); *s = c = getc(fp->fp)){
                                  /* for all lines *not* beginning with 'E' */
    fp->fpos += 1;
    if (fgets(s+1,RDMC_MAXLINE-1,fp->fp) == NULL)               /* read the line */
      return EOF;
  }

  if (c != EOF)
    ungetc(*s, fp->fp);

#if 0
  fp->enr++;                                       /* increase event number */
#endif

  return 0;

} /* function skipevt_ascii() */


/****************************************************************************/
/* functions for reading the dumand-like files                              */
/****************************************************************************/

int rdmc_a_V(const char *s, mcfile *fp) /* scans the V card  */
{ 

  int form = 0;

  if(sscanf(s,"V %i ",&form) != 1) 
    return ERRLINE;
  fp->fmajor = form;
  fp->fminor = 0;

  return 0;

} /* function rdmc_a_V() */

/****************************************************************************/

static int rdmc_a_C(const char *s, char **comment) /* read comment line */
{ 

  while (*s == ' ') 
    s++;
  
  if (toupper(*s != 'C'))
    return 1;

  s += 2;

  /*now append it */
  rdmc_concat_comment(comment,s,DUMAND_ASCII_F);

  return 0;

} /* function rdmc_a_C() */

/****************************************************************************/

static int rdmc_a_G(const char *s, array *ar) /* scans the G line  */
{ 

  int form;

  form = sscanf(s,"G %i %i %i %f %f %f", 
		&(ar->id), &(ar->nch), &(ar->nstr), 
		&(ar->lattitude), &(ar->longitude), &(ar->depth));

  switch(form) {
  case 0: ar->id = 0;
  case 1: ar->nch = 0;
  case 2: ar->nstr = 0;
  case 3: ar->lattitude = 0.0;
  case 4: ar->longitude = 0.0;
  case 5: ar->depth = 0.0;
  default: break;
  } /* switch form */

  return 0;

} /* function rdmc_a_G() */

/****************************************************************************/

static int rdmc_a_M(const char *s, mcfile *fp)
{
  int form;
  char id;
  int irun,igeo,igen,igtrack,daswrun,date;
  char mc_version[RDMC_MAXLINE];

  if (s[0] != 'M' ) 
    return -1;

  if( (form = sscanf(s,"M %c %s %i %i %i %i %i %i ",
	&id,mc_version,&irun,&igeo,&igen,&igtrack,&daswrun,&date) )  != 8) {
    if (form != 7 ){
      return ERRLINE;
    }
    else {
      date= 700101;  
    }
  }
  switch (id) {
  case 'D': /*allowed id's */
  case 'U':
  case 'G':
  case 'M':
    fp->info.dum.mc_id = id;
    break;
  case 'B': /* old data */
    fp->info.dum.mc_id = 'D';
    break;
  case 'J': /* old MC */
  case 'S':
    fp->info.dum.mc_id = 'M';
    break;
  default: /* unknown */
    fp->info.dum.mc_id = 'U';
    break;
  }
  strcpy(fp->info.dum.mc_vers,mc_version);
  fp->info.dum.nrun = irun;
  fp->info.dum.igen =igen;
  fp->info.dum.igtrack = igtrack ;
  fp->info.dum.time = rdmc_o_dateconv(date);
  fp->info.dum.daswrun = daswrun;
  switch (fp->info.dum.mc_id) {
  case 'G':
    fp->info.dum.igeo =0;
    fp->info.dum.igtrack =-1;
    break;
  case 'D': /* data */
    fp->info.dum.igtrack =-1;
    fp->info.dum.igen =0;
    fp->info.dum.daswrun =-1;
    break;
  case 'M':
    fp->info.dum.igtrack =-1;
    break;
  }
  return 0;
  
} /* function rdmc_a_M() */

/****************************************************************************/

static int rdmc_a_P(const char *s, array *ar)
{

  int form;
  long ch, type, serial, str,  iclust;
  float x,y,z;
  const float cs_ori[3] = {0.0, 1.0, -1.0};
  float thresh, rsense, ori;

  form = sscanf(s,"P %li %li %li %li %f %f %f %f %li %f %f"
		,&ch
		,&type
		,&serial
		,&str
		,&x
		,&y
		,&z
		,&ori
		,&iclust
		,&thresh
		,&rsense);
  if (form == 0) return 0;                         /* not a "P" line, ignore */
  if ((form <8) || (form >11)) return ERRLINE;
  ch--;                                     /* I want a 'C' style of indizes */
  if ((ch < 0) ||(ch >= RDMC_MAXCHANNELS) || (ch >= ar->nch))
    return ERRLINE;

  ar->str[ch] = str;
  ar->type[ch] = type;
  ar->x[ch] = x/1000.0;
  ar->y[ch] = y/1000.0;
  ar->z[ch] = z/1000.0;
  if ((ori < 0) || (ori > 2))
    return ERRLINE;

  if (form > 8)                             /* if there was an icluster info */
    ar->clust[ch] = iclust - 1;                 /* fill it as channel number */
  else                                                /* if no icluster info */
    ar->clust[ch] = -1;

  ar->costh[ch] = cs_ori[(int)ori];
  ar->serial[ch] = serial;
  if (form > 9)
    ar->thresh[ch] = thresh;
  else
    ar->thresh[ch] = RDMC_SMALL;
  if (form > 10)
    ar->sensit[ch] = rsense;
  else
    ar->sensit[ch] = 1.0;

  /* now set the geo cal flag */
  ar->is_calib.geo=1;
  return 0;

} /* function rdmc_a_P() */

/****************************************************************************/

static int rdmc_a_K(const char *s, array *ar)
{

  int form;
  array_calib_t cal;
  int ich;

  
  form = sscanf(s,"K %i %f %f %f %f %f %f",
		&ich,
		&(cal.beta_t), &(cal.t_0), &(cal.alpha_t),
		&(cal.ped), &(cal.beta_a), &(cal.kappa));
  switch (form) {
  case 1: cal.beta_t = 0.0;
  case 2: cal.t_0 = 0.0;
  case 3: cal.alpha_t = 0.0;
  case 4: cal.ped = 0.0;
  case 5: cal.beta_a = 0.0;
  case 6: cal.kappa = 0.0;
  case 7: 
    cal.flag = RDMC_CALIB_TDC | RDMC_CALIB_ADC ;
    break;
  default:
    return 0;
  }

  if ((ich > 0) && (ich <= RDMC_MAXCHANNELS))
    ar->cal[ich-1] = cal;

  /* now set the  cal flag */
  ar->is_calib.tdc=1;
  ar->is_calib.adc=1;
  
  return 0;

} /* function rdmc_a_K() */

/****************************************************************************/

static int rdmc_a_Q(const char *s, array *ar)
{
  int form;
  array_hdef_t trig;
  char id;

  rdmc_init_array_hdef(&(trig));

  form = sscanf(s,"Q %c %s %s %s %s %s",
		&(id), trig.tag
		, trig.pars[0], trig.pars[1]
		, trig.pars[2], trig.pars[3] );
  switch (form) {
  case 3: 
  case 4: 
  case 5: 
  case 6: 
  case 2: 
    break;
  default:
    return 0;
  }

  trig.id  = ar->n_trigger+1;
  trig.nwords =  0;
  trig.npars = form -2 ;

  if (ar->n_trigger+1 < RDMC_MAXTRIGID){
    rdmc_add_trigger_def(ar,&trig,ar->n_trigger);
  }

  return 0;

} /* function rdmc_a_Q() */

/****************************************************************************/

static int rdmc_a_E(const char *s, mevt *ev)
{

  char trig[RDMC_MAXLINE];
  int form; 
  unsigned long trigger; 
  int sec;
  unsigned nsec;
  unsigned DUMtLSW;                                             /* dummies */
  int DUMtMSW, rperiodtime;                                      /* dummies */

  form = sscanf(s,"E %i %i %i %u %i %u %i %li %g %i %s ",
                &(ev->enr),&(ev->nhits),
                &sec, &nsec, 
                &DUMtMSW, &DUMtLSW, 
		&rperiodtime,
	        &(trigger), &(ev->t_offset), &(ev->ntrack),
                trig);
  switch( form) {
  case 2:  sec = 0;
  case 3:  nsec = 0;
  case 4:  DUMtMSW = 0;
  case 5:  DUMtLSW = 0;
  case 6:  rperiodtime = 0;
  case 7: trigger = 0;
  case 8:  ev->t_offset=0.;
  case 9:  ev->ntrack = 0;
  case 10: strcpy(trig,"0");
  case 11:
    break;
  default:
    return ERRLINE;
    break;
    }

  if (sec == -1) sec = 0;           /* -1 means: no time was stored */
  if (DUMtMSW == -1) {DUMtMSW = 0; DUMtLSW = 0; }
  if (rperiodtime == -1) rperiodtime = 0;

  ev->mjd = rdmc_unixtime_to_mjd(sec);            /* convert unix time into mjd */
  ev->secs = rdmc_unixtime_to_mjdsecs(sec);           /* and seconds after mjd */
  ev->nsecs = nsec;

  /* now add the trigger */
  {
    mevt_special_t ptrig;
    int itrig;
    rdmc_init_mevt_special(&ptrig,0);
    for (itrig = 0 ; itrig <RDMC_MAXTRIGID ; itrig++){
      if ( (trigger) & (1  << itrig)){
	ptrig.id = itrig;
	rdmc_add_trigger(ev, &ptrig, ev->ntrig, itrig);
      }
    }
  }


  return 0;
  
} /* function rdmc_a_E() */

/****************************************************************************/

static int rdmc_a_H(const char *s, mhit *h)
{

  int form;
  char htrig[RDMC_MAXLINE];

  form = sscanf(s,"H %i %i %g %g %g %s %i %i",
		&(h->str),&(h->ch),
		&(h->tot),
		&(h->amp),
		&(h->t),
		htrig,
		&(h->mt),
		&(h->ma));
  switch (form) {
  case 5:
    strcpy(htrig, "U0");
  case 6:
    h->mt = RDMC_NA;
  case 7:
    h->ma = RDMC_NA;
  case 8:
    break;
  default:
    return ERRLINE;
  }
  h->ch--;                        /* for 'c' style of channel number (0,...) */
#if 0 /* this is stupid */
  h->mt = abs(h->mt);         /* it's possible to have negative track id's */
  h->ma = abs(h->ma);         /* it's possible to have negative track id's */
#endif
  return 0;

} /* function rdmc_a_H() */

/****************************************************************************/

static int rdmc_a_T(const char *s, mtrack *t)
{

  int form;

  rdmc_init_mtrack(t); 
  form  = sscanf(s,"T %i %f %f %f %f %f %f %f %f %i %f "
		 ,&(t->tag)
		 ,&(t->x),&(t->y),&(t->z)
		 ,&(t->px),&(t->py),&(t->pz)
		 ,&(t->e),&(t->t) 
		 ,&(t->id),&(t->length));
  switch (form)
    {
    case 10:
      t->length = RDMC_NA;
    case 11:
      t->length *= 1e-3;       /* dumand format: millimeter, in rdmc: meter */
      if (t->length < 0.)
	t->length = RDMC_NA;  /*10. * t->e/100000. ; */
      t->x *= 1.e-3;
      t->y *= 1.e-3;
      t->z *= 1.e-3;
      rdmc_tauinv_tr(t);                         /* calc cos theta and phi */
      break;
    default:
      return ERRLINE;
    }
  
  return 0;

} /* function rdmc_a_T() */

/****************************************************************************/

static int rdmc_a_F(const char *s, mtrack *t, mevt_special_t *krabi)
{

  int form;

  form  = sscanf(s,"F %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %i %g"
		 ,&(krabi->val[JK_FITID])
		 ,&(t->x),&(t->y),&(t->z)
		 ,&(t->px),&(t->py),&(t->pz)
		 ,&(t->e)
		 ,&(t->t) 
		 ,&(krabi->val[JK_RCHI2])
		 ,&(krabi->val[JK_PROB])
		 ,&(krabi->val[JK_SIGTH])
		 ,&(krabi->val[JK_COVMIN])
		 ,&(krabi->val[JK_COVMAX])
		 ,&(krabi->val[JK_CUTFLAG])
		 ,&(t->id) 
		 ,&(t->length) );


  krabi->val[JK_CHI2]=RDMC_NA;
  switch (form)   /* minimum 10 maximal 17 */
    {
    case 10:
      krabi->val[JK_PROB] =  RDMC_NA;
    case 11:
      krabi->val[JK_SIGTH] = RDMC_NA ;
    case 12:
      krabi->val[JK_COVMIN] =  RDMC_NA;
    case 13:
      krabi->val[JK_COVMAX] = RDMC_NA ;
    case 14:
      krabi->val[JK_CUTFLAG] = RDMC_NA  ;
    case 15:
      t->id = MUON_MINUS ;
    case 16:
      t->length = RDMC_NA;
    case 17:

      t->length *= 1e-3;       /* dumand format: millimeter, in rdmc: meter */
      if (t->length < 0.)
	t->length = RDMC_NA;  /*10. * t->e/100000. ; */
      t->x *= 1.e-3;
      t->y *= 1.e-3;
      t->z *= 1.e-3;
      rdmc_tauinv_tr(t);                               /* calc cos theta and phi */
      break;
    default:
      return ERRLINE;
    }
  
  return 0;

} /* function rdmc_a_F() */

/****************************************************************************/
/* rdmc_a_U() reads optional user fields from a comment line                */
/*  it reallocates the the user field                                       */
/****************************************************************************/

static int rdmc_a_U(const char *s, mevt_special_t **user, int *nuser)
{
  char *t=NULL;
  char *tp;
  mevt_special_t temp_user;
  int iuser,imod;
  int i;
  int r=0;

  if (s[0] != 'U')
    return 0;

#ifndef CRAY
  t = alloca(strlen(s)+1);
#else
  t = malloc(strlen(s)+1);
#endif
  strcpy(t,s);
  
  rdmc_init_mevt_special(&temp_user, RDMC_MAXUSER_PER_LINE);
  
  /* get number of previous these users */
  iuser = 0;
  for (i = 0; i < *nuser; i++){
    iuser += (*user)[i].nval;
  }
  /* iuser is now the index for the new user or the number of previous ones */

  /* now get the last user into temp_user */
  if (iuser > 0){
    rdmc_cp_mevt_special(&temp_user, &((*user)[*nuser-1]));
    r |= rdmc_del_mevt_special(user, nuser, (*nuser)-1);
    imod = temp_user.nval;
  }  else {
    imod = 0;
    temp_user.id = iuser % DUMAND_USER_MAX;
  }

  /* now parse the line and add users */
  for (tp = strtok(t," \t\n"); tp && (iuser < DUMAND_USER_MAX);
       tp = strtok(NULL," \t\n")) {
    iuser++;
    /* check if current user is filled */
    if (imod >= DUMAND_MAXUSER_PER_LINE) {
      r |= rdmc_add_mevt_special(user,nuser,&temp_user,*nuser); 
      rdmc_init_mevt_special(&temp_user,RDMC_MAXUSER_PER_LINE);
      temp_user.id = iuser%DUMAND_USER_MAX;
      imod = 0;
    }
    /* now fill the values */
    temp_user.val[imod] = atof(tp);
    temp_user.nval = imod++;
  }
  /* now commit a last uncommitted block */
  if (iuser > 0) 
    r |= rdmc_add_mevt_special(user,nuser,&temp_user,*nuser); 

#ifdef CRAY
  free(t);
#endif
  return r;                             /* there should not occure an error */

} /* function rdmc_a_U() */

#if 0
/****************************************************************************/
/* a_skip_C() skips all comment lines and points to the next line in ASCII  */
/****************************************************************************/
static int a_skip_C(mcfile *fp)
{
  
  int c;
  char s[RDMC_MAXLINE];

  while ((c = getc(fp->fp)) == 'C'){
    ungetc(c, fp->fp);
    fp->fpos += 1;
    fgets(s, RDMC_MAXLINE-1, fp->fp);
  }

  if (c == EOF) 
    return EOF;
  else {
    ungetc(c, fp->fp);
    return 0;
  }

} /* function a_skip_C() */
#endif

/****************************************************************************/
/* declare_rdmc_user() declares a default user tag               */
/****************************************************************************/

static int  declare_rdmc_dumand_user(array *ar){
  array_hdef_t dumand_user_def;
  int i;
  int r=0;
  int iblock=0;

  rdmc_init_array_hdef(&dumand_user_def);

  for (i=0 ; i < DUMAND_USER_MAX ;i++){
    dumand_user_def.nwords++;
    sprintf(dumand_user_def.words[i%DUMAND_MAXUSER_PER_LINE],"u%i",i+1);
    if (((i+1)%DUMAND_MAXUSER_PER_LINE)==0){
      iblock++;
      dumand_user_def.id=iblock;
      sprintf(dumand_user_def.tag,"dumand%i",iblock);
      r |= rdmc_add_user_def(ar, &dumand_user_def, ar->n_user);
      rdmc_init_array_hdef(&dumand_user_def);
    }
  }
  /* clean up the rest */
  if (((i+1)%DUMAND_MAXUSER_PER_LINE) !=0 ){
    iblock++;
    dumand_user_def.id=iblock;
    sprintf(dumand_user_def.tag,"dumand%i",iblock);
    r |= rdmc_add_user_def(ar, &dumand_user_def, ar->n_user);
  }

  return r;
}

/****************************************************************************/
/* parse_fit_def_ascii() parses the creator line for "recoos"               */
/****************************************************************************/

static int parse_fit_def_ascii(array *ar, const char *s)
{
  char *line=NULL;                  /* temporary copy of the string */
  char *reco;                            /* reco pointer  */
  char *eol;                                         /* end of current line */
  char *bol;                                       /* beginning of the line */
  array_hdef_t dum_fitdef;
  
  if(s){
#ifndef CRAY
    line = alloca(sizeof(char)*(1+strlen(s)));
#else
    line = malloc(sizeof(char)*(1+strlen(s)));
#endif
    strcpy(line,s);                                       /* copy the string */

    for (bol = eol = line; eol; bol = eol+1) {              /* for each line */
      eol = strchr(bol,'\n');                             /* search line end */
      if (eol) *eol = '\0';                          /* replace it with '\0' */
      if (*bol != '\0'){                         /* if it contains something */
	if( (reco=strstr(bol,"recoos (")) != NULL ){ 
	  /* there is a recoos in this history line */ 
	  rdmc_jk_to_fitdef(&(dum_fitdef),ar->n_fit+1);
	  rdmc_add_fit_def(ar,&(dum_fitdef), ar->n_fit);
	}
      }
    }
#ifdef CRAY
    free(line);
#endif
  }
  return 0;

} /* parse_fit_def_ascii() */



#endif /* DUMAND_ASCII_F */
/****************************************************************************/
/********************************** E O F ***********************************/
/****************************************************************************/
/* 
   This is just for EMACS:
   Local Variables:
   compile-command: "cd .. ; make -k rdmc" 
   End: 
*/
