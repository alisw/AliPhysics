
/********implements functions for the mevt_uses structure *********/

#include <string.h>
#include <stdlib.h>
#include <malloc.h>

#include "rdmc.h"

static int comp_uses(const void *h1, const void *h2);


void rdmc_init_mevt_uses(mevt_uses_t *u){
  u->hitid = RDMC_NA;
  u->useid = RDMC_NA;
  return;
}

void rdmc_free_mevt_uses(mevt_uses_t *u){
  return;
}

void rdmc_clear_mevt_uses(mevt_uses_t *u){
  rdmc_free_mevt_uses(u);
  rdmc_init_mevt_uses(u);
  return;
}

void rdmc_cp_mevt_nuses(mevt_uses_t *out, mevt_uses_t *in, int nuses){

  if ((in == NULL) || (out == NULL) || (in == out))
    return;
  memcpy(out,in,nuses*sizeof(mevt_uses_t));
  return;
}
  

/* sorts the uses array for a) uses id and the hit_id */
void rdmc_sort_uses(mevt_uses_t *use, int nuse){
  qsort(use,nuse,sizeof(mevt_uses_t),comp_uses);
}

static int comp_uses(const void *h1, const void *h2){
  int s;
  const mevt_uses_t *u1=h1,*u2=h2;
  s = u1->useid - u2->useid;
  switch (s){
    case 0 : 
      s =  u1->hitid - u2->hitid;
      switch (s){
      case 0 : 
	s= -1;
	break;
      default:
	break;
      }
      break;
  default:
    break;
  }
  return s;
}


int rdmc_create_fit_uses(mevt *e, int ifit){
  mevt_uses_t use;
  int i;

  rdmc_init_mevt_uses(&use);
  use.useid=ifit;
  
  /* loop over all hits in the old event */
  for (i=0 ; i < e->nhits ; i++ ){
    use.hitid=e->h[i].id;
    rdmc_add_mevt_uses(&(e->fit_uses),&(e->nfit_uses)
			,&use,e->nfit_uses); 
  } /* all old hits */
  return 0;
}

int rdmc_create_trig_uses(mevt *e, int itrig){
  mevt_uses_t use;
  int i;

  rdmc_init_mevt_uses(&use);
  use.useid=itrig;
  
  /* loop over all hits in the old event */
  for (i=0 ; i < e->nhits ; i++ ){
    use.hitid=e->h[i].id;
    rdmc_add_mevt_uses(&(e->trig_uses),&(e->ntrig_uses)
			,&use,e->ntrig_uses); 
  } /* all old hits */
  return 0;
}


 /* count the number of channels which are used by fit or trigger i*/
int rdmc_count_fit_nch_uses(mevt *e, array *a, int ifit){
  int nch=0;
  static int *ich=NULL;
  int iuse;

  /* first sort fit_uses data for fit number */
  rdmc_sort_uses(e->fit_uses, e->nfit_uses);

#if 1
  /* get pointer to first element for fit i */
  for (iuse=0 ; iuse < e->nfit_uses ; iuse++){
    if (e->fit_uses[iuse].useid == ifit){
      break; /* ok we found one */
    }
  }
  if (iuse >= e->nfit_uses) return 0; 
#else
  iuse=0;
#endif

  /* get temporary array for the number of channels and init to 0*/
  ich = alloca(a->nch*sizeof(int));
  memset(ich,0,a->nch*sizeof(int));

  /* flag  each channel number which was hit */
  while (iuse < e->nfit_uses) {
    if(e->fit_uses[iuse].useid == ifit){
      int ihit =  rdmc_mevt_ihit(e, e->fit_uses[iuse].hitid );
      if (ihit == RDMC_NA ) {
#if 1 	
	/* consider this hit as non existent for good reason */
#else
	nch++; /* assume a different channel was hit */
#endif
#if 0
	rdmc_warnprintf("ev=%i, fit=%i hit id=%i is not existent ",e->enr,ifit,e->fit_uses[iuse].hitid);
#endif
      } else if ( ich[e->h[ihit].ch] == 0){
	nch++;   /* count all newly hit channels */
	ich[e->h[ihit].ch] += 1;
      }
#if 1
    }else{
      break;
#endif
    } 
    iuse++;
  }

  return nch;
}






