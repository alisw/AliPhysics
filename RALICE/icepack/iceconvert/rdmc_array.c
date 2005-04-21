
/*
 *  functions for the array structure 
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rdmc.h"

#include "amanda.h"
#include "baikal.h"
#include "dumand.h"
#include "uwi.h"


/****************************************************************************/
/* The function rarr() reads the array info of a file                       */
/****************************************************************************/

int rdmc_rarr(mcfile *fp, array *ar)
{
  int r;
  
  switch(fp->format) {
#ifdef DUMAND_ASCII_F
  case DUMAND_ASCII_F:       /* if it is a dumand-like format */
    r = rdmc_rarr_ascii(fp, ar);                                /* read it */
    break;
#endif
#ifdef UWI_ASCII_F
  case UWI_ASCII_F:       /* if it is the UWI format: read default geometry */
    r = rdmc_rarr_uwi(NULL, ar);
    break;
#endif
#ifdef AMANDA_ASCII_F
  case AMANDA_ASCII_F:       /* if it is a amanda-like format */
    r = rdmc_rarr_amanda(fp, ar);                           /* read it */
    break;
#endif
#ifdef BAIKAL_BIN_F
  case BAIKAL_BIN_F:                       /* if it is a baikal-like format */
    r = rdmc_rarr_baikal_mc(fp, ar);                  /* read it here */
    fp->fpos = ftell(fp->fp);
    break;
#endif
  default: 
    r = RDMC_UNKNOWN_FORMAT;
  }
  if ((r != RDMC_IO_OK) && (r != EOF)) {
    fp->errline = r;              /* set the source code line of the error */
#if 0
    r = RDMC_ILF;                /* set the return value in case of format error */
#endif
  }

  return r;

} /* funtion rarr() */

/****************************************************************************/
/* warr() writes the array info to a file                                   */
/****************************************************************************/

int rdmc_warr(mcfile *fp, const array *ar)
{
  int r;

  switch(fp->format) {
#ifdef DUMAND_ASCII_F
  case DUMAND_ASCII_F:                            /* for a dumand-like file */
    r = rdmc_warr_ascii(fp, ar);
    break;
#endif
#ifdef UWI_ASCII_F
  case UWI_ASCII_F:                                       /* for a UWI file */
    rdmc_whd_uwi(fp);                        /* store the file header */
    r = rdmc_warr_uwi(NULL, ar);
    break;
#endif
#ifdef AMANDA_ASCII_F
  case AMANDA_ASCII_F:                            /* for a amanda-like file */
    r = rdmc_warr_amanda(fp, ar);
    break;
#endif
#ifdef BAIKAL_BIN_F
  case BAIKAL_BIN_F:                              /* for a baikal-like file */
    r = rdmc_warr_baikal_mc(fp, ar); 
    break;
#endif
  default:
    r = RDMC_UNKNOWN_FORMAT;
  }

  return r;

} /* function warr() */

/****************************************************************************/
/* init_array() resets the event values to start values                     */
/****************************************************************************/

void rdmc_init_array(array *ar)
{

  int i;

  /* init the generic array structure */
  ar->id = AMANDA;
  ar->nch = 0;
  ar->nstr = 0;
  ar->longitude = -180.0;
  ar->lattitude = -90.0;
  ar->depth =  RDMC_DEFAULT_DEPTH;
  ar->tbegin=RDMC_NA;
  ar->tend=RDMC_NA;
  ar->nrun=0;

  ar->comment = NULL;

  ar->array_id = 0;


  for (i = 0; i < RDMC_MAXCHANNELS; i++) {
    ar->str[i] = ar->clust[i] = ar->serial[i] = 0;
    ar->type[i] = RDMC_DEFAULT_OM;
    ar->x[i] = ar->y[i] = ar->z[i] = 0.0;
    ar->costh[i] = -1.0;
    ar->thresh[i] =  RDMC_SMALL_THRESH;
    ar->sensit[i] = 1.0;
  } /* for i */

  /* now the calib structure */
  rdmc_init_array_calib_stat(&(ar->is_calib));

  for (i = 0; i < RDMC_MAXCHANNELS; i++) {
    rdmc_init_array_calib(&(ar->cal[i]));
  } /* for i */

  rdmc_init_array_calib_utc(&(ar->cal_utc));

  /* now init the trigger structure */
  ar->n_trigger = 0;
  ar->def_trig=NULL;

  /* now the user init */
  ar->n_user=0;
  ar->def_user=NULL;

  /* now the fit_def */
  ar->n_fit=0;
  ar->def_fit = NULL;

  /* now the status init */
  ar->n_stat=0;
  ar->def_stat=NULL;

  /* now the status init */
  ar->n_mcinfo=0;
  ar->def_mcinfo=NULL;

  /* now the stack */
  ar->tmp=NULL;

} /* function init_array() */

/****************************************************/
/* add specific header_def elements to ar */

int rdmc_add_user_def(array *ar, array_hdef_t *user, int iuse){
 return rdmc_add_array_hdef( &(ar->def_user), 
				    &(ar->n_user),
				    user,
				    iuse);
}
int rdmc_del_user_def(array *ar , int iuse){
  return  rdmc_del_array_hdef(
			      &(ar->def_user),
			      &(ar->n_user),
			      iuse);
}
int rdmc_add_stat_def(array *ar, array_hdef_t *stat, int istat){
  return rdmc_add_array_hdef( &(ar->def_stat), 
				    &(ar->n_stat),
				    stat,
				    istat);
}
int rdmc_del_stat_def(array *ar , int istat){
  return  rdmc_del_array_hdef(
			      &(ar->def_stat),
			      &(ar->n_stat),
			      istat);
}
int rdmc_add_fit_def(array *ar, array_hdef_t *fit, int ifit){
  return rdmc_add_array_hdef( &(ar->def_fit), 
				    &(ar->n_fit),
				    fit,
				    ifit);
}
int rdmc_del_fit_def(array *ar , int ifit){
  return  rdmc_del_array_hdef(
			      &(ar->def_fit),
			      &(ar->n_fit),
			      ifit);
}
int rdmc_add_trigger_def(array *ar, array_hdef_t *trig, int itrig){
  /* trig_def is a static array -> workaround */
  return rdmc_add_array_hdef( &(ar->def_trig), 
			      &(ar->n_trigger),
			      trig,
			      itrig);
}


int rdmc_del_trig_def(array *ar , int itrig){
  return  rdmc_del_array_hdef(
			      &(ar->def_trig),
			      &(ar->n_trigger),
			      itrig);
}


void rdmc_clear_array(array *ar){
  unsigned long int array_id;
  /* Free all the dynamically allocated memory you can in the structure ar. */
  rdmc_free_array(ar);

  array_id=ar->array_id;                         /* save the id */
  rdmc_init_array(ar); /* Reset the array structure to initializarion values */
  ar->array_id =  array_id + 1;

}/* function rdmc_clear_array() */

void rdmc_free_array(array *ar){
  int i;
  /* Free all the dynamically allocated memory you can in the structure ar. */
  
  if(ar->def_user != NULL){
    for (i=0 ; i < ar->n_user ; i++){
      rdmc_free_array_hdef(&(ar->def_user[i]));
    }
    ar->n_user=0;
    free(ar->def_user);
    ar->def_user=NULL;
  }

  if(ar->def_trig != NULL){
    for (i=0 ; i < ar->n_trigger ; i++){
      rdmc_free_array_hdef(&(ar->def_trig[i]));
    }
    ar->n_trigger=0;
    free(ar->def_trig);
    ar->def_trig=NULL;
  }

  if(ar->def_fit != NULL){
    for (i=0 ; i < ar->n_fit ; i++){
      rdmc_free_array_hdef(&(ar->def_fit[i]));
    }
    ar->n_fit=0;
    free(ar->def_fit);
    ar->def_fit=NULL;
  }
  if(ar->def_stat != NULL){
    for (i=0 ; i < ar->n_stat ; i++){
      rdmc_free_array_hdef(&(ar->def_stat[i]));
    }
    free(ar->def_stat);
    ar->def_stat=NULL;
    ar->n_stat=0;
  }
  if(ar->def_mcinfo  != NULL){
    for (i=0 ; i < ar->n_mcinfo ; i++){
      rdmc_free_array_hdef(&(ar->def_mcinfo[i]));
    }
    ar->def_mcinfo=NULL;
    free(ar->def_mcinfo);
    ar->n_mcinfo=0;
  }
  

  if (ar->comment != NULL){
    free(ar->comment);
    ar->comment=NULL;
  }
  
  if (ar->tmp){
    free(ar->tmp);
    ar->tmp=NULL;
  }
}/* function rdmc_free_array() */
                   

void rdmc_cp_array(array *out, const array *in){/* copies two geometries */
  int i;

  /* just cp everything and cleanup later */
  memcpy(out,in,sizeof(array));

  if(out->comment != NULL){
    out->comment = malloc(strlen(in->comment)*sizeof(char));
    strcpy(out->comment,in->comment);
  }  

  if (out->tmp){ /* if there is a new stack delete it */
    out->tmp=NULL;
  }

  if(out->n_user >0){
    out->def_user = malloc(sizeof(array_hdef_t)*out->n_user);
    for (i=0 ; i < out->n_user ; i++)  
      rdmc_cp_hdef(&(out->def_user[i]),&(in->def_user[i]));
  }

  if(out->n_trigger >0){
    out->def_trig =  malloc(sizeof(array_hdef_t)*out->n_trigger);
    for (i=0 ; i < out->n_trigger ; i++)  
      rdmc_cp_hdef(&(out->def_trig[i]),&(in->def_trig[i]));
  }

  if(out->n_fit >0){
    out->def_fit =  malloc(sizeof(array_hdef_t)*out->n_fit);
    for (i=0 ; i < out->n_fit ; i++)  
      rdmc_cp_hdef(&(out->def_fit[i]),&(in->def_fit[i]));
  }

  if(out->n_stat >0){
    out->def_stat = malloc(sizeof(array_hdef_t)*out->n_stat);
    for (i=0 ; i < out->n_stat ; i++)  
      rdmc_cp_hdef(&(out->def_stat[i]),&(in->def_stat[i]));
  }

  if(out->n_mcinfo >0){
    out->def_mcinfo =  malloc(sizeof(array_hdef_t)*out->n_mcinfo);
    for (i=0 ; i < out->n_mcinfo ; i++)  
      rdmc_cp_hdef(&(out->def_mcinfo[i]),&(in->def_mcinfo[i]));
  }

}

/****************************************************************************/
/* comp_array() compares two array for geometry etc                         */
/*     returns 0 if they are equal, and 1 if not                            */
/* it looks only for: number of channels, of strings, of the string number  */
/* and the coordinates of each pmt, not for id or sensitivity etc           */
/****************************************************************************/

int rdmc_comp_array(const array *a1, const array *a2)
{

  int i;                                                   /* loop variable */
  static const double maxdiff = 0.01;      /* max geometric deviation of pmt */
  int r=RDMC_ARRAY_COMP_OK;

  if (a1->id != a2->id) r |= RDMC_ARRAY_COMP_HEADER;
  if (a1->longitude != a2->longitude) r |=RDMC_ARRAY_COMP_HEADER;
  if (a1->lattitude != a2->lattitude) r |=RDMC_ARRAY_COMP_HEADER;
  if (a1->depth != a2->depth) r |=RDMC_ARRAY_COMP_HEADER;
  if (a1->array_id != a2->array_id)  r |=RDMC_ARRAY_COMP_HEADER;
  if( a1->nrun != a2->nrun)  r |=RDMC_ARRAY_COMP_HEADER;
#if 0 /* this is ignored */
  if (a1->tmp != a2->tmp)  r |=RDMC_ARRAY_COMP_HEADER;
#endif

  if (a1->nch != a2->nch) r |= RDMC_ARRAY_COMP_GEO;
  if (a1->nstr != a2->nstr) r |= RDMC_ARRAY_COMP_GEO;

  for (i = 0; i < a1->nch; i++) {
    if (a1->str[i] != a2->str[i]) r |= RDMC_ARRAY_COMP_GEO;
    if (fabs(a1->x[i] - a2->x[i]) > maxdiff) r |= RDMC_ARRAY_COMP_GEO;
    if (fabs(a1->y[i] - a2->y[i]) > maxdiff)  r |= RDMC_ARRAY_COMP_GEO;
    if (fabs(a1->z[i] - a2->z[i]) > maxdiff)  r |= RDMC_ARRAY_COMP_GEO;
    if (fabs(a1->costh[i] - a2->costh[i]) > maxdiff)  r |= RDMC_ARRAY_COMP_GEO;
    if (a1->type[i] != a2->type[i] ) r|= RDMC_ARRAY_COMP_OMS;
    if (a1->serial[i] != a2->serial[i] ) r|= RDMC_ARRAY_COMP_OMS;
    if (fabs(a1->thresh[i] - a2->thresh[i] )) r|= RDMC_ARRAY_COMP_OMS;
    if (fabs(a1->sensit[i] - a2->sensit[i] )) r|= RDMC_ARRAY_COMP_OMS;
    
    if (fabs(a1->cal[i].beta_t - a2->cal[i].beta_t )) r|= RDMC_ARRAY_COMP_CALIB;
    if (fabs(a1->cal[i].t_0 - a2->cal[i].t_0 )) r|= RDMC_ARRAY_COMP_CALIB;
    if (fabs(a1->cal[i].alpha_t - a2->cal[i].alpha_t )) r|= RDMC_ARRAY_COMP_CALIB;
    if (fabs(a1->cal[i].ped - a2->cal[i].ped )) r|= RDMC_ARRAY_COMP_CALIB;
    if (fabs(a1->cal[i].beta_a - a2->cal[i].beta_a )) r|= RDMC_ARRAY_COMP_CALIB;
    if (fabs(a1->cal[i].kappa - a2->cal[i].kappa )) r|= RDMC_ARRAY_COMP_CALIB;
    if (fabs(a1->cal[i].ped_tot - a2->cal[i].ped_tot )) r|= RDMC_ARRAY_COMP_CALIB;
    if (fabs(a1->cal[i].beta_tot - a2->cal[i].beta_tot )) r|= RDMC_ARRAY_COMP_CALIB;
    if (fabs(a1->cal[i].kappa_tot - a2->cal[i].kappa_tot )) r|= RDMC_ARRAY_COMP_CALIB;
    if (a1->cal[i].flag != a2->cal[i].flag ) r|= RDMC_ARRAY_COMP_CALIB;
  } /* for i */

  if (memcmp(&(a1->is_calib),&(a2->is_calib),sizeof(array_calib_stat_t))
      != 0)  
    r |= RDMC_ARRAY_COMP_CALIB;


  if (a1->n_trigger != a2->n_trigger)  
    r |= RDMC_ARRAY_COMP_TRIGGER;
  else{
    for (i=0 ; i < a1->n_trigger ; i++){
      if ( rdmc_comp_array_hdef( 
				&(a1->def_trig[i]),
				&(a2->def_trig[i]))  
	   ){
	r |= RDMC_ARRAY_COMP_TRIGGER;
	break;
      }
    }
  }

  if (a1->n_user != a2->n_user)  
    r |= RDMC_ARRAY_COMP_USER;
  else{
    for (i=0 ; i < a1->n_user ; i++){
      if( rdmc_comp_array_hdef( &(a1->def_user[i]),
				 &(a2->def_user[i])) ){
	r |= RDMC_ARRAY_COMP_USER;
	break;
      }
    }
  }

  if (a1->n_fit != a2->n_fit)  
    r |= RDMC_ARRAY_COMP_FIT;
  else{
    for (i=0 ; i < a1->n_fit ; i++){
      if( rdmc_comp_array_hdef( &(a1->def_fit[i]),
				 &(a2->def_fit[i])) ){
	r |= RDMC_ARRAY_COMP_FIT;
	break;
      }
    }
  }


  if (a1->n_stat != a2->n_stat)  
    r |= RDMC_ARRAY_COMP_STATUS;
  else{
    for (i=0 ; i < a1->n_stat ; i++){
      if( rdmc_comp_array_hdef( &(a1->def_stat[i]),
				 &(a2->def_stat[i])) ){
	r |= RDMC_ARRAY_COMP_STATUS;
	break;
      }
    }
  }

  if (a1->n_mcinfo != a2->n_mcinfo)  
    r |= RDMC_ARRAY_COMP_MC;
  else{
    for (i=0 ; i < a1->n_mcinfo ; i++){
      if ( rdmc_comp_array_hdef( &(a1->def_mcinfo[i]),
			      &(a2->def_mcinfo[i]))
	   ){
	r |= RDMC_ARRAY_COMP_MC;
	break;
      }
    }
  }


  return r;

} /* function comp_array() */


/*************************************************************************/
/* Copies the values for OM:                                             */
/* in_i from arry in into position out_i in array out                    */
/*************************************************************************/

void rdmc_channel_cp(array  *out,int out_i, const array *in ,int in_i)
{ 
  out->str[out_i] = in->str[in_i];
  out->clust[out_i] = in->clust[in_i];
  out->x[out_i] = in->x[in_i];
  out->y[out_i] = in->y[in_i];
  out->z[out_i] = in->z[in_i];
  out->costh[out_i] = in->costh[in_i];
  out->type[out_i] = in->type[in_i];
  out->serial[out_i] = in->serial[in_i];
  out->thresh[out_i] = in->thresh[in_i];
  out->sensit[out_i] = in->sensit[in_i];

  /* the calibration for this channel */
  out->cal[out_i] = in->cal[in_i];
  
  return;
}

/* in order to create objects one needs unique id's*/
/* these functions return them */

int rdmc_unique_trigger_id(char *tag, const array *ar,const  char *tag_root){
  rdmc_unique_hdef_tag(tag, ar->def_trig, ar->n_trigger, tag_root);
  return rdmc_unique_hdef_id(ar->def_trig, ar->n_trigger);
}
int rdmc_unique_status_id(char *tag, const array *ar, const char *tag_root){
  rdmc_unique_hdef_tag(tag, ar->def_stat, ar->n_stat, tag_root);
  return rdmc_unique_hdef_id(ar->def_stat, ar->n_stat);
}
int rdmc_unique_user_id(char *tag, const array *ar, const char *tag_root){
  rdmc_unique_hdef_tag(tag, ar->def_user, ar->n_user, tag_root);
  return rdmc_unique_hdef_id(ar->def_user, ar->n_user);
}





