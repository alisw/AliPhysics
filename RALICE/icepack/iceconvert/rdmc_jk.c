
/* this c-file implemtes the usage of a special FResult type */
/* it is call jk in honor of Jaanus Krabi who used the values first */


#include <string.h>

#include "rdmc.h"


const array_hdef_t rdmc_jk_fit_def = 
{
  0, "rdmc-jk" , 8 , 0, { "id" , "rchi2" , "prob", "sigth" 
		       , "covmin", "covmax", "cutflag" , "chi2"}
  , { "" }
};

void rdmc_jk_to_fitdef(array_hdef_t *jk_def, int id){
  rdmc_cp_hdef(jk_def,&rdmc_jk_fit_def);
  jk_def->id = id;
}

void rdmc_init_fit_jk(mevt_special_t *r,int id){

  r->id = id;
  r->nval = 8;
  r->val[JK_FITID]=RDMC_NA;
  r->val[JK_CUTFLAG] = -1;
  r->val[JK_CHI2] 
    = r->val[JK_RCHI2] 
    = r->val[JK_PROB] 
    = r->val[JK_SIGTH] 
    = r->val[JK_COVMIN] 
    = r->val[JK_COVMAX] 
    = RDMC_NA;

}

int rdmc_is_fresult_jk(const array *ar,const mevt *ev, int ifit){
  int idef;
  if ((ifit >= ev->nfit)
      || (ifit < 0))
    return 0;
  idef = ev->fresult[ifit].id;
  if ((idef >= ar->n_fit)
      || (idef < 0))
    return 0;
  return rdmc_is_this_jk(&(ar->def_fit[idef]),&(ev->fresult[ifit]));
}


int rdmc_is_this_jk(const array_hdef_t *rd,const mevt_special_t *res){
  int ret = 0;
  if ( (strcmp(rd->tag,"rdmc-jk") == 0)
       && (res->nval == 8 )){
    ret=1;
  }
  return ret;
      
}

int rdmc_is_fitdef_jk(const array *ar,int idef){
  int ret = 0;
  if ( (idef<0) || (idef>= ar->n_fit) )
    ret=0;
  else{
    if (strcmp(ar->def_fit[idef].tag,"rdmc-jk") == 0)
      ret=1;
  }
  
  return ret;
      
}

