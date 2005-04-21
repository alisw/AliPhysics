
/******** implements functions for the mhit structure *********/

#include <string.h>
#include <stdlib.h>

#include "rdmc.h"


void rdmc_init_mhit(mhit *h){
  h->ch = h->str = h->mt = h->ma = h->id = RDMC_NA ;
  h->t= -RDMC_BIG;
  h->amp = h->tot = RDMC_NA ;
  h->weight = 1.;
  rdmc_init_mhit_stat(&(h->hstat));
  h->nuser=0;
  h->user=NULL;
  h->comment=NULL;
  h->tmp=NULL;
}

void rdmc_init_mhit_stat(mhit_stat_t *s){
  s->n_tdc_edges=0;
  s->tdc_flag=0;
}

void rdmc_clear_mhit(mhit *h){
  rdmc_free_mhit(h);
  rdmc_init_mhit(h);
}

void rdmc_free_mhit(mhit *h){
  int i;
  if (h->user !=NULL){
    for(i=0 ; i < h->nuser ; i++){
      rdmc_free_mevt_special(&(h->user[i]));
    }
    free(h->user);
    h->nuser=0;
    h->user=NULL;
  }
  if (h->tmp !=NULL){
    free(h->tmp);
    h->tmp=NULL;
  }
  if (h->comment !=NULL){
    free(h->comment);
    h->comment=NULL;
  }
}


/****************************************************************************
 * merge two hits.  the amplitude of the merged hit is added only if ch=-1
 * else the amplitude of the channel ch is used.
 * If there is no such channel, a zero amplitude is assumed
 ****************************************************************************/

void rdmc_merge_hits(mhit *in_1, mhit *in_2, int ch)
{

  float ampa,tfirst,tota;
  int idt,ida;

  if (in_1 == in_2) return; /* the hits are the SAME! */

  if ((ch == -1) || ((ch == in_1->ch) && (ch == in_2->ch) )){ /* use both */
    if ((in_1->amp >= 0.) && (in_2->amp >= 0.))
      ampa=  in_1->amp + in_2->amp;  /* sum of amplitudes*/
    else if (in_1->amp >= 0.)
      ampa = in_1->amp;
    else if  (in_2->amp >= 0.)
      ampa = in_2->amp;
    else
      ampa = -1.;
    if (ampa > 0)
      ida = (in_1->amp >= in_2->amp) ? (in_1->ma) : (in_2->ma) ; 
    else
      ida = -1;  /* amp muon id*/
  } else if ((ch == in_1->ch) && (in_1->amp >= 0.)) { /* only first */
    ampa = in_1->amp;
    ida = (ampa > 0)? in_1->ma : -1;
  } else if ((ch == in_2->ch) && (in_2->amp >= 0.)) { /* only second */
    ampa = in_2->amp;
    ida = (ampa > 0)? in_2->ma : -1;
    in_1->ch = ch; /* set correct (for amp) channel number */
  }  else {
    ampa = -1.0;
    ida = -1;
  }

  tfirst = (in_1->t <= in_2->t ) ? (in_1->t) : (in_2->t) ; /* keep first time*/
  idt = (in_1->t <= in_2->t ) ? (in_1->mt) : (in_2->mt) ; /* time muon id*/

  in_1->amp = ampa;  /* new amplitude */
  in_2->amp = 0.;

  if ((in_1->tot >= 0.)&&(in_2->tot >= 0.))
    tota=  in_1->tot + in_2->tot;  /* sum of totlitudes*/
  else
    if (in_1->tot >= 0.)
      tota = in_1->tot;
    else
      if  (in_2->tot >= 0.)
	tota = in_2->tot;
      else
	tota =-1.;
  if (ida <0 )
    if (tota > 0 )
      ida = (in_1->tot >= in_2->tot ) ? (in_1->ma) : (in_2->ma) ; 
  /* amp muon id*/

  in_1->tot = tota;  /* new amplitude */
  in_2->tot = 0.;

  in_1->t = tfirst;
  in_2->t = 0.;
  in_1->mt = idt;
  in_2->mt = 0;
  in_1->ma = ida;
  in_2->ma = 0;

  in_1->hstat.n_tdc_edges = 
    (in_1->hstat.n_tdc_edges >=  in_2->hstat.n_tdc_edges )
    ? in_1->hstat.n_tdc_edges
    : in_2->hstat.n_tdc_edges;

  in_1->hstat.tdc_flag = 
    (in_1->hstat.tdc_flag >=  in_2->hstat.tdc_flag )
    ? in_1->hstat.tdc_flag
    : in_2->hstat.tdc_flag;

  in_1->weight = in_1->weight * in_2->weight ;

  return;

} /* rdmc_merge_hits() */

/****************************************************************************
 * copies the hit *in to *out by overwriting.
 * If in == out, nothing is done.
 ****************************************************************************/
void rdmc_cp_mhit(mhit *out, mhit *in)
{
  if ((in == NULL) || (out == NULL) || (in == out))
    return;
  
  memcpy(out, in, sizeof(mhit));
  if (out->nuser){
    int i;
    mevt_special_t *from,*to;
    out->user = (mevt_special_t *) malloc(in->nuser*sizeof(mevt_special_t));
    if (!(out->user)){
      out->nuser = 0;
    } else {
      from = in->user;
      to = out->user;
      for (i = 0 ; i < in->nuser ; i++){
	rdmc_cp_mevt_special(to,from);
	to++ ; from++;
      }
    }
  }

  if (in->comment) {
    out->comment = (char *) malloc(strlen(in->comment)+1);
    if (!out->comment)
      strcpy(out->comment, in->comment);
  }

  if (out->tmp){ /* if there is a new stack delete it */
    out->tmp=NULL;
  }


} /* rdmc_cp_mhit() */
