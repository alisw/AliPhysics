
/********implements functions for the mevt structure *********/
/* $Header: /net/local/cvsroot/siegmund/rdmc/rdmc_mevt.c,v 1.27 2004/02/19 17:10:10 wiebusch Exp $ */

#include <string.h>
#include <stdlib.h>

#include "rdmc.h"
#include "rdmc_local.h"

#include "amanda.h"
#include "dumand.h"
#include "uwi.h"
#include "baikal.h"

static int  rdmc_add_mevt_WF(waveform **list, int *count, waveform *new, int ipos);
static int  rdmc_del_mevt_WF(waveform **list, int *count, int ipos);

/****************************************************************************/
/* The function revt() reads an event of a mc file                          */
/*     if it was not successfull, try to skip to the next event             */
/****************************************************************************/

int rdmc_revt(mcfile *fp, array *ar, mevt *ev)
{

  int r;                                                   /* return value */

  if (feof(fp->fp))                                /* test for end of file */
    return EOF;

  switch (fp->format) {
#ifdef DUMAND_ASCII_F
  case DUMAND_ASCII_F:
    r = rdmc_revt_ascii(fp, ev,ar);
    break;
#endif
#ifdef UWI_ASCII_F
  case UWI_ASCII_F:
    r = rdmc_revt_uwi(fp, ev, ar);
    break;
#endif
#ifdef BAIKAL_BIN_F
  case BAIKAL_BIN_F:                      /* if it is a baikal-like format */
    r = rdmc_revt_baikal(fp, ev, ar);
    break;
#endif
#ifdef AMANDA_ASCII_F
  case AMANDA_ASCII_F:
    r = rdmc_revt_amanda(fp, ev, ar);
    break;
#endif
  default: 
    r = RDMC_UNKNOWN_FORMAT;
  } /* switch fp->format */

  if ((r != 0) && (r != EOF)) {               /* if illegal format occured */
#if 0
    long filepos=fp->fpos;              /* shows the current file position */
                                         /* store the old file position */
#endif
    fp->errline = r;              /* set the source code line of the error */
    if (fp->sloppy){
      rdmc_skipevt(fp);                            /* skip to the next event */
#if 0
      fp->fpos = filepos;                       /* restore old file position */
#endif
    }
#if 0
    r = RDMC_ILF;     /* set the return value in case of format error */
#endif
#if 0
    if (r)
      rdmc_err_print(fp,r);
#endif
  }

  return r;

} /* function revt() */


/****************************************************************************/
/* skipevt skips over the next event record in a mc/data file               */
/****************************************************************************/
int rdmc_skipevt(mcfile *fp)
{
  int r;

  switch(fp->format) {
#ifdef DUMAND_ASCII_F
  case DUMAND_ASCII_F:                        /* for the dumand-like format */
    r = rdmc_skipevt_ascii(fp);                                /* skip the event */
    break;
#endif
#ifdef UWI_ASCII_F
  case UWI_ASCII_F:                        /* for the UWI-like format */
    r = rdmc_skipevt_uwi(fp);
    break;
#endif
#ifdef AMANDA_ASCII_F
  case AMANDA_ASCII_F:                        /* for the amanda-like format */
    r = rdmc_skipevt_amanda(fp);                               /* skip the event */
    break;
#endif
#ifdef BAIKAL_BIN_F
  case BAIKAL_BIN_F:                          /* for the baikal-like format */
    r = rdmc_skipevt_baikal_mc(fp);                                   /* skip the event */
    break;
#endif
  default:
    r = RDMC_UNKNOWN_FORMAT;
    break;
  }
  return r;
} /* function skipevt() */


/****************************************************************************/
/* The function wevt() writes an event of a mc file                         */
/****************************************************************************/

int rdmc_wevt(mcfile *fp, const mevt *ev, const array *ar)
{
  int r;

#if (DEBUG == 1)
    fprintf(stderr, "<Event(write): enr=%i nch=%i %s %s (%s)>\n",
	   ev->enr,ev->nch,(ev->gen!=NULL)?"gen":"",
	   (ev->rec!=NULL)?"rec":"",(fp->mc!=0)?"mc":"data");
#endif

  switch(fp->format) {
#ifdef DUMAND_ASCII_F
  case DUMAND_ASCII_F:                            /* for a dumand-like file */
    r = rdmc_wevt_ascii(fp, ev,ar );
    break;
#endif
#ifdef UWI_ASCII_F
  case UWI_ASCII_F:                                       /* for a UWI file */
    r = rdmc_wevt_uwi(fp, ev, ar);
    break;
#endif
#ifdef AMANDA_ASCII_F
  case AMANDA_ASCII_F:                            /* for a amanda-like file */
    r = rdmc_wevt_amanda(fp, ev, ar);
    break;
#endif
#ifdef BAIKAL_BIN_F
  case BAIKAL_BIN_F:                              /* for a baikal-like file */
    if (fp->info.bai.mc)                          /* if it is a monte carlo file */
      r = rdmc_wevt_baikal_mc(fp, ev, ar);
    else                                            /* if it is a data file */
      r = rdmc_wevt_baikal_data(fp, ev, ar);
    fp->fpos = ftell(fp->fp);
    break;
#endif
  default:
    r = RDMC_UNKNOWN_FORMAT;
  }

  return r;

} /* function wevt() */

/*******************************************/
/* Initialization functions ***********/
/*******************************************/
/**** init =  set to default values ****/
/***  clear = free also memory    *****/
/****************************************/


/****************************************************************************/
/* clear_mevt() frees the memory used by the pointers of a event and        */
/*              resets the mevt structure itself                            */
/****************************************************************************/

void rdmc_clear_mevt(mevt *ev)
{
  unsigned long int event_id;

  rdmc_free_mevt(ev);

  event_id=ev->event_id;                         /* save the event-id */
  rdmc_init_mevt(ev);                              /* re-init the event */
  ev->event_id = event_id+1;                 /* increment the event id */

} /* function clear_mevt() */

void rdmc_free_mevt(mevt *ev){
  int i;
  if (ev->h != NULL){                                 /* if there were hits */
    for (i =0 ; i < ev->nhits ; i++)
      rdmc_free_mhit(&(ev->h[i]));
    free(ev->h);                                               /* free them */
    ev->nhits=0;
    ev->h=NULL;
  }
  if (ev->wf != NULL){                            /* if there were waveforms */
    for (i =0 ; i < ev->nwf ; i++)
      rdmc_free_WF(&(ev->wf[i]));

    free(ev->wf);                                               /* free them */
    ev->nwf = 0;
    ev->wf  = NULL;
  }

  if (ev->gen != NULL){                  /* if there were generating tracks */
    for (i =0 ; i < ev->ntrack ; i++)
      rdmc_free_mtrack(&(ev->gen[i]));
    free(ev->gen);                                               /* free it */
    ev->ntrack=0;
    ev->gen=NULL;
  }
  if ((ev->rec != NULL)&&(ev->fresult != NULL)){ /* if there are recos */
    for (i =0 ; i < ev->nfit ; i++){
      rdmc_free_mtrack(&(ev->rec[i]));
      rdmc_free_mevt_special(&(ev->fresult[i]));
    }
    free(ev->rec);                                               /* free it */
    free(ev->fresult);
    ev->nfit=0;
    ev->rec=NULL;
    ev->fresult=NULL;
  }

  if (ev->user != NULL){
    for (i =0 ; i < ev->nuser ; i++)
      rdmc_free_mevt_special(&(ev->user[i]));
    ev->nuser=0;
    free(ev->user);
    ev->user=NULL;
  }

  if (ev->status != NULL){
    for (i =0 ; i < ev->nstat ; i++)
      rdmc_free_mevt_special(&(ev->status[i]));
    ev->nstat=0;
    free(ev->status);
    ev->status=NULL;
  }
  if (ev->ptrig != NULL){
    for (i =0 ; i < ev->ntrig ; i++)
      rdmc_free_mevt_special(&(ev->ptrig[i]));
    ev->ntrig=0;
    free(ev->ptrig);
    ev->ptrig=NULL;
  }
  if (ev->mcinfo != NULL){
    for (i =0 ; i < ev->nmcinfo ; i++)
      rdmc_free_mevt_special(&(ev->mcinfo[i]));
    ev->nmcinfo=0;
    free(ev->mcinfo);
    ev->mcinfo=NULL;
  }

#if !NEW_USES
  if (ev->trig_uses != NULL){
    free(ev->trig_uses);
    ev->ntrig_uses = 0;
  }
  if (ev->fit_uses != NULL){
    free(ev->fit_uses);
    ev->nfit_uses = 0;
  }
#endif

  if (ev->comment != NULL){
    free(ev->comment);
    ev->comment=NULL;
  }

  if (ev->tmp){
    free(ev->tmp);
    ev->tmp=NULL;
  }
} /* function free_mevt() */


/****************************************************************************/
/* init_mevt() resets the event values to start values                      */
/****************************************************************************/

void rdmc_init_mevt(mevt *ev)
{
   
  ev->nrun = RDMC_NA;
  ev->enr = RDMC_NA;
  ev->mjd = 0;
  ev->secs = 0;
  ev->nsecs = 0;

  ev->trigger = 0;
  ev->trigid = 0;
  ev->t_offset = 0.0;
  ev->event_id=0;
  ev->sort_status=RDMC_NO_SORT;
  ev->ee = 0;
  ev->weight = 1.;

  ev->nch = 0;
  ev->nstr = 0;

  ev->nhits     = 0;
  ev->nwf       = 0;
  ev->ntrack    = 0;
  ev->nfit      = 0;
  ev->nuser     = 0;
  ev->nstat     = 0;
  ev->ntrig     = 0;
  ev->nmcinfo   = 0;
  ev->ntrig_uses =0;
  ev->nfit_uses  =0;

  ev->h          = NULL;
  ev->wf         = NULL;
  ev->gen        = NULL;
  ev->rec        = NULL;
  ev->fresult    = NULL;
  ev->user       = NULL;
  ev->status     = NULL;
  ev->ptrig      = NULL;
  ev->mcinfo     = NULL;
  ev->trig_uses  = NULL;
  ev->fit_uses   = NULL;
  ev->comment    = NULL;
  ev->tmp        = NULL;
   
} /* function init_mevt() */



/*************************************************************************/
/* counts number of strings in event e                                   */
/*************************************************************************/

int rdmc_count_nstr(const mevt *ev)
{
  /* too speed up we need to do some dirty memory management ! */
  /* there are usually less hits than MAX_CHANNELS */
  int i;
  int nstr = 0;
  static int do_init=1;
  static int hitted[RDMC_MAXCHANNELS];

  if(do_init){
    for (i = 0; i < RDMC_MAXCHANNELS; i++) /* maximum number of strings */
      hitted[i] = 0;
    do_init=0;
  }

  /* now fill the block */
  for (i = 0; i < ev->nhits; i++){
    if ((ev->h[i].str <= RDMC_MAXCHANNELS) && (ev->h[i].str > 0))
      hitted[ev->h[i].str-1]++;
  }

  /* now evaluate and clean on the fly -> strange loop  */
  for (i = 0; i < ev->nhits; i++){
    if ((ev->h[i].str <= RDMC_MAXCHANNELS) && (ev->h[i].str > 0)){
      if(hitted[ev->h[i].str-1]){
	hitted[ev->h[i].str-1]=0;
	nstr++;
      }
    }
  }

  return nstr;

} /* rdmc_count_nstr() */


/*************************************************************************/
/* counts number of channels in event e                                  */
/*************************************************************************/
int rdmc_count_nch(const mevt *ev)
{
  /* too speed up we need to do some dirty memory management ! */
  /* there are usually less hits than MAX_CHANNELS */
  int i;
  int nch = 0;
  static int do_init=1;
  static int hitted[RDMC_MAXCHANNELS];

  if(do_init){
    for (i = 0; i < RDMC_MAXCHANNELS; i++)
      hitted[i] = 0;
    do_init=0;
  }

  /* now fill the block */
  for (i = 0; i < ev->nhits; i++){
    if ((ev->h[i].ch < RDMC_MAXCHANNELS) && (ev->h[i].ch >= 0))
      hitted[ev->h[i].ch]++;
  }

  /* now evaluate and clean on the fly -> strange loop  */
  for (i = 0; i < ev->nhits; i++){
    if ((ev->h[i].ch < RDMC_MAXCHANNELS) && (ev->h[i].ch >= 0)){
      if(hitted[ev->h[i].ch]){
	hitted[ev->h[i].ch]=0;
	nch++;
      }
    }
  }
  return nch;
} /* rdmc_count_nch() */

/****************************************************************************/
/* ns_evt() calculates the number of strings and fills it in the nstr field */
/*         it also fills the string number into the h[].str fields          */
/****************************************************************************/

int rdmc_fill_mhit_str(mevt *ev, const array *ar)
{

  int i;
  int error_flag = 0;                                 /* indicates an error */
  int geocal;                          /* flag if array is a valid geometry */

  if (ev->nch >= RDMC_MAXCHANNELS) return RDMC_INCONSISTENT_GEOMETRY;
                                          /* if there are too much channels */
  geocal = ((ar->is_calib.geo) && (ar->nch >0)) ? 1 : 0;

  for (i = 0; i < ev->nhits; i++) {              /* for all hits */
    
    if ( (geocal) && 
	 (ev->h[i].ch >= 0) && 
	 (ev->h[i].ch < ar->nch))    /* if ch in range */
      ev->h[i].str = ar->str[ev->h[i].ch];           /* write string number */
    else {                                       /* if channel not in range */
      ev->h[i].str = 0;                          /* errornous string number */
    } /* if not (0 < ch < ar->nch) */
  } /* for i */

  return (error_flag == 0)? ev->nstr:error_flag;/* return number of strings */

} /* function rdmc_fill_mhit_str() */

/****************************************************************************/
/* event_cp                                                               */
/* copy values of channel in into out                                       */
/* returns 0 if OK                                                          */
/****************************************************************************/
int rdmc_cp_mevt(mevt *out,mevt *in)
{
  int r=0;
#if 0
  rdmc_init_mevt(out);
#else
  rdmc_clear_mevt(out);
#endif

  /* just copy everything and overwrite dynamic pointers later */
  memcpy(out,in,sizeof(mevt)); 

  if (in->comment) {
    out->comment = (char *) malloc(strlen(in->comment)+1);
    if (!out->comment) {
      r |= -1;
    }else
      strcpy(out->comment, in->comment);
  }

  if (out->nhits){
    mhit *from,*to;
    int i;
    out->h = (mhit *) malloc(in->nhits*sizeof(mhit));
    if (!(out->h)) {
      out->nhits = 0;
      r |= -1;
    }else {
      from = in->h;
      to = out->h;
      for (i = 0 ; i < in->nhits ; i++){
	rdmc_cp_mhit(to,from);
	to++;
	from++;
      }
    }
  }


  if (out->nwf){                                     /*   copy waveforms  */
    waveform *from,*to;
    int i;
    out->wf = (waveform *) malloc(in->nwf*sizeof(waveform));
    if (!(out->wf)) {
      out->nwf = 0;
      r |= -1;
    }else {
      from = in->wf;
      to = out->wf;
      for (i = 0 ; i < in->nwf ; i++){
	rdmc_cp_WF(to,from);
	to++;
	from++;
      }
    }
  }

  if (out->ntrack){
    int i;
    mtrack *from,*to;
    out->gen = (mtrack *) malloc(in->ntrack*sizeof(mtrack));
    if (!(out->gen)) {
      out->ntrack=0;
       r |= -1;
    }else {
      from = in->gen;
      to = out->gen;
      for (i = 0 ; i < in->ntrack ; i++){
	rdmc_cp_mtrack(to,from);
	to++;
	from++;
      }
    }
  }

  if (out->nfit){
    int i;
    mtrack *ft_from,*ft_to;
    mevt_special_t *fr_from,*fr_to;
    out->rec = (mtrack *) malloc(in->nfit*sizeof(mtrack));
    if (!(out->rec)) {
      out->nfit=0;
      r |= -1;
    }else {
      ft_from = in->rec;
      ft_to = out->rec;
      out->fresult = (mevt_special_t *) malloc(in->nfit*sizeof(mevt_special_t));
      if (!(out->fresult)){
	free(out->rec);
	out->rec=NULL;
	out->nfit=0;
	r |= -1;
      } else {
	fr_from = in->fresult;
	fr_to = out->fresult;
	for (i = 0 ; i < in->nfit ; i++){
	  rdmc_cp_mtrack(ft_to,ft_from);
	  rdmc_cp_mevt_special(fr_to,fr_from);
	  fr_to++ ; fr_from++ ; ft_to++ ; ft_from++;
	}
      }
    }
  }

  if (out->nuser){
    int i;
    mevt_special_t *from,*to;
    out->user = (mevt_special_t *) malloc(in->nuser*sizeof(mevt_special_t));
    if (!(out->user)){
      out->nuser =0;
      r |= -1;
    } else {
      from = in->user;
      to = out->user;
      for (i = 0 ; i < in->nuser ; i++){
	rdmc_cp_mevt_special(to,from);
	to++ ; from++;
      }
    }
  }

  if (out->nstat){
    int i;
    mevt_special_t *from=in->status,*to;
    out->status = (mevt_special_t *) malloc(in->nstat*sizeof(mevt_special_t));
    if (! (to=out->status)){
      out->nstat =0;
      r |= -1;
    }
    for (i = 0 ; i < in->nstat ; i++){
      rdmc_cp_mevt_special(to,from);
      to++ ; from++;
    }
  }

  if (out->ntrig){
    int i;
    mevt_special_t *from=in->ptrig,*to;
    out->ptrig = (mevt_special_t *) malloc(in->ntrig*sizeof(mevt_special_t));
    if (! (to=out->ptrig)){
      r |= -1;
      out->ntrig =0;
    }else {
      for (i = 0 ; i < in->ntrig ; i++){
	rdmc_cp_mevt_special(to,from);
	to++ ; from++;
      }
    }
  }

  if (out->nmcinfo){
    int i;
    mevt_special_t *from=in->mcinfo,*to;
    out->mcinfo = (mevt_special_t *) 
      malloc(in->nmcinfo*sizeof(mevt_special_t));
    if (! (to=out->mcinfo)){
      r |= -1;
      out->nmcinfo =0 ;
    }
    for (i = 0 ; i < in->nmcinfo ; i++){
      rdmc_cp_mevt_special(to,from);
      to++ ; from++;
    }
  }

  if (out->ntrig_uses){
    out->trig_uses = (mevt_uses_t *)
      malloc(in->ntrig_uses*sizeof(mevt_uses_t));
    if (! (out->trig_uses) ){
      r |= -1;
      out->ntrig_uses =0;
    }else
      rdmc_cp_mevt_nuses(out->trig_uses,in->trig_uses,in->ntrig_uses);
  }

  if (out->nfit_uses){
    out->fit_uses = (mevt_uses_t *)
      malloc(in->nfit_uses*sizeof(mevt_uses_t));
    if (! (out->fit_uses) ){
      r |= -1;
      out->nfit_uses =0;
    }else
      rdmc_cp_mevt_nuses(out->fit_uses,in->fit_uses,in->nfit_uses);
  }

  if (out->tmp){ /* if there is a new stack delete it */
    out->tmp=NULL;
  }

  return r;
}


/***************** Functions to add/remove mevt objects *******/

int rdmc_add_gen(mevt *ev, mtrack *tr, int itrack){
  return rdmc_add_mevt_mtrack( &(ev->gen), 
			       &(ev->ntrack),
			       tr,
			       itrack);
}
int rdmc_del_gen(mevt *ev , int itrack){
  return  rdmc_del_mevt_mtrack( &(ev->gen),
				&(ev->ntrack),
				itrack);
}

int rdmc_add_fit(mevt *ev, mtrack *ft, mevt_special_t *fr, int ifit){
  int rf=1,rr=1;
  rf = rdmc_add_mevt_mtrack( &(ev->rec), 
			     &(ev->nfit),
			     ft,
			     ifit);
  if (!rf){
    ev->nfit--; /* set counter back  for fresult  */
    rr=rdmc_add_mevt_special( &(ev->fresult), 
			      &(ev->nfit),
			      fr,
			      ifit);
    if (rr) {
      ev->nfit++; /* increase to get the last track back out */
      rdmc_del_mevt_mtrack( &(ev->rec),
			    &(ev->nfit),
			    ifit);
    }
  }
  return (rr || rf) ;
}

int rdmc_del_fit(mevt *ev , int ifit){
  int rf=1,rr=1;
  rf = rdmc_del_mevt_mtrack( &(ev->rec), 
			     &(ev->nfit),
			     ifit);
  if (!rf){
    ev->nfit++; /* set counter back  for fresult  */
    rr=rdmc_del_mevt_special( &(ev->fresult), 
			     &(ev->nfit),
			     ifit);
    if (rr) {
      ev->nfit--;
      /* cannot add it again */
    }
  }
  return (rr || rf) ;
}

int rdmc_add_user(mevt *ev, mevt_special_t *us, int iu){
  return rdmc_add_mevt_special( &(ev->user), 
			       &(ev->nuser),
			       us,
			       iu);
}

int rdmc_del_user(mevt *ev , int iu){
  return  rdmc_del_mevt_special( &(ev->user),
				&(ev->nuser),
				iu);
}

int rdmc_add_status(mevt *ev, mevt_special_t *st, int is){
  return rdmc_add_mevt_special( &(ev->status), 
				&(ev->nstat),
				st,
				is);
}

int rdmc_del_status(mevt *ev , int is){
  return  rdmc_del_mevt_special( &(ev->status),
				&(ev->nstat),
				is);
}

int rdmc_add_mcinfo(mevt *ev, mevt_special_t *mci, int imci){
  return rdmc_add_mevt_special( &(ev->mcinfo), 
				&(ev->nmcinfo),
				mci,
				imci);
}
int rdmc_del_mcinfo(mevt *ev , int imci){
  return  rdmc_del_mevt_special( &(ev->mcinfo),
				 &(ev->nmcinfo),
				 imci);

}
int rdmc_add_trigger(mevt *ev, mevt_special_t *tr, int itr, int idef){
  int imask;

  imask=tr->id;
  if ((imask < 0)||(idef <0))
    return 1;

  if (imask >= RDMC_MAXTRIGID)
    ev->trigger |= ((unsigned long)1<<(RDMC_MAXTRIGID-1));
  else  
    ev->trigger |= ((unsigned long)1<<(imask));

  if (idef >= RDMC_MAXTRIGID)
    ev->trigid |= ((unsigned long)1<<(RDMC_MAXTRIGID-1));
  else  
    ev->trigid |= ((unsigned long)1<<(idef));

  return rdmc_add_mevt_special( &(ev->ptrig), 
				&(ev->ntrig),
				tr,
				itr);
}

int rdmc_del_trigger(mevt *ev , int itr, int idef){
  int imask;
  if ((itr < 0 ) || (idef < 0) )
    return 1;

  imask=ev->ptrig[itr].id;
  if (imask <0 )
    return 1;
  
  if (imask < RDMC_MAXTRIGID)
    ev->trigger &= (~(1<<(imask)));
  /* in the else part we keep the overflow bit because it is not uniqe */
  if (idef < RDMC_MAXTRIGID)
    ev->trigid &= (~(1<<(idef)));
  /* in the else part we keep the overflow bit because it is not uniqe */

  return  rdmc_del_mevt_special( &(ev->ptrig),
				 &(ev->ntrig),
				 itr);
}


int rdmc_add_trig_uses(mevt *ev, mevt_uses_t *tu, int ituse){
  return rdmc_add_mevt_uses( &(ev->trig_uses), 
			     &(ev->ntrig_uses),
			     tu,
			     ituse);
}

int rdmc_del_trig_uses(mevt *ev , int ituse){
  return rdmc_del_mevt_uses( &(ev->trig_uses), 
			     &(ev->ntrig_uses),
			     ituse);
}

int rdmc_add_fit_uses(mevt *ev, mevt_uses_t *fu, int ifuse){
  return rdmc_add_mevt_uses( &(ev->fit_uses), 
			     &(ev->nfit_uses),
			     fu,
			     ifuse);
}

int rdmc_del_fit_uses(mevt *ev , int ifuse){
  return rdmc_del_mevt_uses( &(ev->fit_uses), 
			     &(ev->nfit_uses),
			     ifuse);
}

int rdmc_add_mhit(mevt *ev, mhit *h, int ihit){
  return rdmc_add_mevt_mhit( &(ev->h), 
			     &(ev->nhits),
			     h,
			     ihit);
}

int rdmc_del_mhit(mevt *ev , int ihit){
  return rdmc_del_mevt_mhit( &(ev->h), 
			     &(ev->nhits),
			     ihit);
}

int rdmc_add_WF(mevt *ev, waveform *wf, int iwf){
  return rdmc_add_mevt_WF( &(ev->wf), &(ev->nwf), wf, iwf);
}

int rdmc_del_WF(mevt *ev , int iwf){
  int rr=1;
  if(iwf >=0 && iwf < ev->nwf){
    rdmc_free_WF(&(ev->wf[iwf]));

    rr = rdmc_del_mevt_WF( &(ev->wf), 
         		 &(ev->nwf),
			   iwf);
  }
  return rr ;
}

int  rdmc_add_mevt_mtrack(mtrack **list, int *count, mtrack *new, int ipos){
  int i;
  if( ( list == NULL )
      || (new == NULL )
      || (count == NULL )
      || (*count < 0)
      || (ipos < 0)
      || (ipos > *count)
      )
    return 1;
  (*count)++;
  *list = (mtrack *)  realloc(*list, sizeof(mtrack)*(*count)); 

  for (i = (*count - 1) ; i > ipos ; i-- ){
    rdmc_cp_mtrack( &((*list)[i]), &((*list)[i-1]) );
  }
  rdmc_cp_mtrack( &((*list)[ipos]), new );

  return 0;
} 

int  rdmc_add_mevt_special(mevt_special_t **list, int *count
			   ,mevt_special_t *new, int ipos){
  int i;
  if( ( list == NULL )
      || (new == NULL )
      || (count == NULL )
      || (*count < 0)
      || (ipos < 0)
      || (ipos > *count)
      )
    return 1;

  (*count)++;
  *list = (mevt_special_t *)  realloc( *list,
				       sizeof(mevt_special_t)*(*count));

  for (i = (*count - 1) ; i > ipos ; i-- ){
    rdmc_cp_mevt_special( &((*list)[i]), &((*list)[i-1]) );
  }
  rdmc_cp_mevt_special( &((*list)[ipos]), new );
  
  return 0;
} 

int  rdmc_add_mevt_mhit(mhit **list, int *count, mhit *new, int ipos){
  int i;
  if( ( list == NULL )
      || (new == NULL )
      || (count == NULL )
      || (*count < 0)
      || (ipos < 0)
      || (ipos > *count)
      )
    return 1;

  (*count)++;
  *list = (mhit *)  realloc( *list, sizeof(mhit)*(*count));

  for (i = (*count - 1) ; i > ipos ; i-- ){
    rdmc_cp_mhit( &((*list)[i]), &((*list)[i-1]) );
  }
  rdmc_cp_mhit( &((*list)[ipos]), new );
  return 0;
} 

static int  rdmc_add_mevt_WF(waveform **list, int *count, waveform *new, int ipos){
  int i;
  if( ( list == NULL )
      || (new == NULL )
      || (count == NULL )
      || (*count < 0)
      || (ipos < 0)
      || (ipos > *count)
      )
    return 1;

  (*count)++;
  *list = (waveform *)  realloc( *list, sizeof(waveform)*(*count));

  for (i = (*count - 1) ; i > ipos ; i-- ){
    rdmc_cp_WF( &((*list)[i]), &((*list)[i-1]) );
  }
  rdmc_cp_WF( &((*list)[ipos]), new );
  return 0;
} 

int  rdmc_add_mevt_uses(mevt_uses_t **list, int *count
			,mevt_uses_t *new, int ipos){
  int i;
  if( ( list == NULL )
      || (new == NULL )
      || (count == NULL )
      || (*count < 0)
      || (ipos < 0)
      || (ipos > *count)
      )
    return 1;

  (*count)++;
  *list = (mevt_uses_t *)  realloc( *list, sizeof(mevt_uses_t)*(*count));

  for (i = (*count - 1) ; i > ipos ; i-- ){
    rdmc_cp_mevt_nuses( &((*list)[i]), &((*list)[i-1]), 1);
  }
  rdmc_cp_mevt_nuses( &((*list)[ipos]), new , 1);
  return 0;
} 


int  rdmc_del_mevt_mtrack(mtrack **list, int *count, int ipos){
  int i;
  if( (list == NULL )
      || (*list == NULL )
      || (*count <= 0)
      || (ipos < 0)
      || (ipos >= *count)
      )
    return 1;
  for (i = ipos ; i < (*count-1) ; i++ ){
    rdmc_cp_mtrack( &((*list)[i]), &((*list)[i+1]) );
  }
  (*count)--;
  if (*count > 0)
    *list = (mtrack *) 
      realloc(*list,sizeof(mtrack)*(*count));
  else {
    free(*list);
    *list=NULL;
  }
  return 0;
} 

int  rdmc_del_mevt_special(mevt_special_t **list, int *count, int ipos){
  int i;
  if( (list == NULL )
      || (*list == NULL )
      || (*count <= 0)
      || (ipos < 0)
      || (ipos >= *count)
      )
    return 1;
  for (i = ipos ; i < (*count-1) ; i++ ){
    rdmc_cp_mevt_special( &((*list)[i]), &((*list)[i+1]) );
  }
  (*count)--;
  if (*count > 0)
    *list = (mevt_special_t *) 
      realloc(*list,sizeof(mevt_special_t)*(*count));
  else {
    free(*list);
    *list=NULL;
  }
  return 0;
} 

int  rdmc_del_mevt_mhit(mhit **list, int *count, int ipos){
  int i;
  if( (list == NULL )
      || (*list == NULL )
      || (*count <= 0)
      || (ipos < 0)
      || (ipos >= *count)
      )
    return 1;
  for (i = ipos ; i < (*count-1) ; i++ ){
    rdmc_cp_mhit( &((*list)[i]), &((*list)[i+1]) );
  }
  (*count)--;
  if (*count > 0)
    *list = (mhit *) 
      realloc(*list,sizeof(mhit)*(*count));
  else {
    free(*list);
    *list=NULL;
  }
  return 0;
} 



static int  rdmc_del_mevt_WF(waveform **list, int *count, int ipos){
  int i;
  if( (list == NULL )
      || (*list == NULL )
      || (*count <= 0)
      || (ipos < 0)
      || (ipos >= *count)
      )
    return 1;
  for (i = ipos ; i < (*count-1) ; i++ ){
    rdmc_cp_WF( &((*list)[i]), &((*list)[i+1]) );
  }
  (*count)--;
  if (*count > 0)
    *list = (waveform *) realloc(*list,(*count)*sizeof(waveform));
  else {
    free(*list);
    *list=NULL;
  }
  return 0;
} 


int  rdmc_del_mevt_uses(mevt_uses_t **list, int *count, int ipos){
  int i;
  if( (list == NULL )
      || (*list == NULL )
      || (*count <= 0)
      || (ipos < 0)
      || (ipos >= *count)
      )
    return 1;
  for (i = ipos ; i < (*count-1) ; i++ ){
    rdmc_cp_mevt_nuses( &((*list)[i]), &((*list)[i+1]), 1 );
  }
  (*count)--;
  if (*count > 0)
    *list = (mevt_uses_t *) 
      realloc(*list,sizeof(mevt_uses_t)*(*count));
  else {
    free(*list);
    *list=NULL;
  }
  return 0;
} 


int rdmc_del_ifit_uses(mevt *e , int fitid){
  int i;
  int r=0;
  for (i=0 ; i < e->nfit_uses ; ){
    if (e->fit_uses[i].useid == fitid )
      r |= rdmc_del_fit_uses(e , i);
    else
      i++;
  }
  return r;
}

int rdmc_del_itrig_uses(mevt *e , int trigid){
  int i;
  int r=0;
  
  for (i=0 ; i < e->ntrig_uses ; ){
    if (e->trig_uses[i].useid == trigid )
      r |= rdmc_del_trig_uses(e , i);
    else
      i++;
  }
  return r;
}

int rdmc_del_ihit_uses(mevt *e , int hitid){
  int i;
  int r=0;
  
  /* remove fit uses */
  for (i=0 ; i < e->nfit_uses ; ){
    if (e->fit_uses[i].hitid == hitid )
      r |= rdmc_del_fit_uses(e , i);
    else
      i++;
  }

  /* remove trig uses */
  for (i=0 ; i < e->ntrig_uses ; ){
    if (e->trig_uses[i].hitid == hitid )
      r |= rdmc_del_trig_uses(e , i);
    else
      i++;
  }
  return r;
}

/* in order to create objects one needs unique id's*/
/* these functions return them */
int rdmc_unique_mhit_id(const mevt *ev){
  int id = ev->nhits-1;
  int is_not_unique;

  do{
    int i;
    id++;
    is_not_unique=0;
    for (i = ev->nhits-1 ; i >=0 ; i-- ){ /* downward to be faster  */
      if (ev->h[i].id == id){
	is_not_unique = 1;
	break;
      }
    } /*for */
  } while (is_not_unique);

  return id;
}

 /* give new uniqe hit ids if hit id <0 */
int  rdmc_repair_mhit_id(mevt *ev){
  int i,j;
  int repaired=0;
  for (i = 0 ; i < ev->nhits ; i++ ){
    if (ev->h[i].id < 0){
      ev->h[i].id = rdmc_unique_mhit_id(ev);
      repaired++;
    }
    else{ /* test for a double count  */
      for (j = 0 ; j  < i ; j++){
	if (ev->h[i].id == ev->h[j].id){
	  ev->h[i].id = rdmc_unique_mhit_id(ev);
	  repaired++;
	}
      } /* for j */
    } /* else */
  }
  if (repaired){ /*  the ids were inconsistent */
    /* -> delete all uses ! */
    if(ev->ntrig_uses){
      ev->ntrig_uses= 0;
      free(ev->trig_uses);
      ev->trig_uses = NULL;
    }
    if(ev->nfit_uses){
      ev->nfit_uses= 0;
      free(ev->fit_uses);
      ev->fit_uses = NULL;
    }
  }
  return repaired;
}

/* in order to create objects one needs unique id's*/
/* these functions return them */
int rdmc_unique_gen_id(const mevt *ev){
  int id = ev->ntrack+1;
  int is_not_unique;

  do{
    int i;
    id++;
    is_not_unique=0;
    for (i = ev->ntrack-1 ; i >=0 ; i-- ){ /* downward to be faster  */
      if (ev->gen[i].tag == id){
	is_not_unique = 1;
	break;
      }
    } /*for */
  } while (is_not_unique);

  return id;
}

 /* give new uniqe hit ids if hit id <0 */
int  rdmc_repair_gen_id(mevt *ev){
  int i,j;
  int repaired=0;
  for (i = 0 ; i < ev->ntrack ; i++ ){
    if (ev->gen[i].tag < 0){
      ev->gen[i].tag = rdmc_unique_gen_id(ev);
      repaired++;
    } else{ /* test for a double count  */
      for (j = 0 ; j  < i ; j++){
	if (ev->gen[i].tag == ev->gen[j].tag){
	  ev->gen[i].tag = rdmc_unique_gen_id(ev);
	  repaired++;
	}
      } /* for j */
    } /* else */
  }
  if (repaired){ /*  the ids were inconsistent */
    /* -> delete all track parent id's */
    for (i = 0 ; i < ev->ntrack ; i++)
      ev->gen[i].parent = RDMC_NA;

    /* -> delete all hit parent id's */
    for (i = 0 ; i < ev->nhits ; i++)
      ev->h[i].mt = ev->h[i].ma =  RDMC_NA;
  }

  return repaired;
}
/* return the index of the hit hitid */
int rdmc_mevt_ihit(mevt *ev, int hitid){
  int i;
  for (i = 0 ; i < ev->nhits ;  i++){
    if (ev->h[i].id == hitid)
      return i;
  }
  return RDMC_NA;
}
