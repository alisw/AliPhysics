
/******** implements functions for the mtrack structure *********/

#include <string.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
char *rdmc_mtrack_cvsid = 
"$Header: /net/local/cvsroot/siegmund/rdmc/rdmc_mtrack.c,v 1.13 2001/06/08 17:23:59 wiebusch Exp $";

#include "rdmc.h"

/****************************************************************************/
/* init_mtrack() resets the event values to start values                    */
/****************************************************************************/

void rdmc_init_mtrack(mtrack *t)
{
  t->id = MUON_PLUS;
  t->e = 0.0;
  t->x = t->y = t->z = 0.0;
  t->t = 0.0 ;
  t->costh = 1.0;
  t->phi = 0.0;
  t->px = t->py = 0;
  t->pz = -1.0;
  t->length = RDMC_NA;
  t->nmuon = 1;
  t->parent = RDMC_NA;
  t->tag = RDMC_NA;
  t->weight = 1.;
  t->nuser=0;
  t->user=NULL;
  t->comment=NULL;
  t->tmp=NULL;
} /* function init_mtrack() */

void rdmc_clear_mtrack(mtrack *t){
  rdmc_free_mtrack(t);
  rdmc_init_mtrack(t);
}

void rdmc_free_mtrack(mtrack *t){
  int i;
  if (t->user !=NULL){
    for(i=0 ; i < t->nuser ; i++){
      rdmc_free_mevt_special(&(t->user[i]));
    }
    free(t->user);
    t->nuser=0;
    t->user=NULL;
  }
  if (t->tmp !=NULL){
    free(t->tmp);
    t->tmp=NULL;
  }
  if (t->comment !=NULL){
    free(t->comment);
    t->comment=NULL;
  }

  rdmc_init_mtrack(t);
}

/*****************************************************************************/
/* tau_tr() calcs the direction cosinuus of a given track and fills it into  */
/*         the track                                                         */
/*****************************************************************************/

void rdmc_tau_tr(mtrack *tr)
{

  double sinth;                                               /* sinus theta */
  
  if (fabs(tr->costh) < 1.0){
    sinth = sqrt( (1 - (double) tr->costh ) * (1 + (double) tr->costh) );
  }else{
    sinth = 0.0;
  }

  tr->px = - sinth * cos(tr->phi);
  tr->py = - sinth * sin(tr->phi);
  tr->pz = - tr->costh;

  return;

} /* function tau_tr() */


/****************************************************************************/
/* tauinv_tr() calcs cos theta and phi from the direction cosinus and fills */
/*               them into the track                                        */
/****************************************************************************/

void rdmc_tauinv_tr(mtrack *tr)
{

  double r;
  
  r = tr->px * tr->px + tr->py *tr->py + tr->pz * tr->pz;
  if (r > 0.0) r = sqrt(r);
  if (r != 0) {
    r = 1/r;
    tr->px *= r;
    tr->py *= r;
    tr->pz *= r;
  } /* if r != 0 */

  tr->costh = - tr->pz;

#if 0 /* this is wrong and not needed anyway .. atan2 can handle this */  
  if (tr->px == 0.0)
    tr->phi = M_PI / 2.0;
  else
#endif
    tr->phi = atan2(- tr->py, - tr->px);

  if (tr->phi < 0)                       /* I want positive phis (0,..2*Pi) */
    tr->phi += 2*M_PI;

} /* function tauinv_tr() */

/****************************************************************************
 * copies the track *in to *out by overwriting.
 * If in == out, nothing is done.
 ****************************************************************************/
void rdmc_cp_mtrack(mtrack *out, const mtrack *in)
{
  if ((in == NULL) || (out == NULL) || (in == out))
    return;
  
  memcpy(out, in, sizeof(mtrack));

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
    if (out->comment) 
      strcpy(out->comment, in->comment);
  }

  if (out->tmp){ /* if there is a new stack delete it */
    out->tmp=NULL;
  }

} /* rdmc_cp_mtrack() */
