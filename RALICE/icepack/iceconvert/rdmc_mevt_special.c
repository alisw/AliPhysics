
/********implements functions for the mevt_special structure *********/

#include <string.h>

#include "rdmc.h"

void rdmc_init_mevt_special(mevt_special_t *s, int ntoken){
  register int i,imax=RDMC_MAXTOKEN_PER_LINE;
  s->id = RDMC_NA;
  s->nval = 0;
  if ((ntoken >= 0) &&  (ntoken < RDMC_MAXTOKEN_PER_LINE) )
    imax=ntoken;
  for (i=0 ; i< imax ; i++)
    s->val[i]=RDMC_NA;
  return;
}

void rdmc_clear_mevt_special(mevt_special_t *s, int ntoken){
  rdmc_free_mevt_special(s);
  rdmc_init_mevt_special(s,ntoken);
}

void rdmc_free_mevt_special(mevt_special_t *s){
}


void rdmc_cp_mevt_special(mevt_special_t *out ,mevt_special_t  *in){
  int nv;
  if ((in == NULL) || (out == NULL) || (in == out))
    return;
  out->id = in->id;
  out->nval= in->nval;
  nv = (in->nval <= RDMC_MAXTOKEN_PER_LINE ) 
    ? in->nval : RDMC_MAXTOKEN_PER_LINE;
  if (nv > 0){
    memcpy(out-> val,in->val,nv*sizeof(float));
    out->nval= nv;
  }else {
    out->nval= 0;
  }
  return;
}
