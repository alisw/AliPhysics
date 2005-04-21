
/* implement functions special for the f2k 1.1 format */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>


#include "rdmc.h"
#include "amanda.h"
#include "f2k.h"

#define USE_TOT_CAL 1 /* this is defined in 1.3 but rdmc always used it in 1.1 already */

#if 0
#define ATOI(X) strtol(X, (char **)NULL, 10);
#endif
#define SEC_PER_DAY 86400

  /* this block defines a uses structuere which is hirachically
     build by uses elements ... only n amanda.c so far */
typedef struct {
  mevt_uses_t *u;
  int nu;
  int id;
} mevt_usesblock_t;

static void rdmc_init_usesblock(mevt_usesblock_t *ub){
  ub->u=NULL;
  ub->nu=0;
  ub->id=0;
}
static void rdmc_clear_usesblock(mevt_usesblock_t *ub){
  if(ub->u)
    free(ub->u);
  rdmc_init_usesblock(ub);
}

const f2000_line_t V_line_1x1 = 
  {"V " , V_LINE , COMP_STRINGWISE , rdmc_f2k_dummy_parser };
const f2000_line_t COMMENT_line_1x1 = 
  {"*!/$%&(][?~+-_:,;@|<>^#\"\\", COMMENT_LINE,COMP_CHARPUNCT,rdmc_f2k_dummy_parser };
const f2000_line_t HI_line_1x1 = 
  {"HI ", HI_LINE, COMP_STRINGWISE, rdmc_amanda_HI_1x1 };
const f2000_line_t ARRAY_line_1x1 = 
  {"ARRAY ", ARRAY_LINE, COMP_STRINGWISE, rdmc_amanda_ARRAY_1x1 };
const f2000_line_t TRIG_DEF_line_1x1 = 
  {"TRIG_DEF ", TRIG_DEF_LINE, COMP_STRINGWISE, rdmc_amanda_TRIG_DEF_1x1 };
const f2000_line_t TRIG_PAR_line_1x1 = 
  {"TRIG_PAR ", TRIG_PAR_LINE, COMP_STRINGWISE, rdmc_amanda_TRIG_PAR_1x1 };
const f2000_line_t FIT_DEF_line_1x1 = 
  {"FIT_DEF " , FIT_DEF_LINE , COMP_STRINGWISE, rdmc_amanda_FIT_DEF_1x1 };
const f2000_line_t STAT_DEF_line_1x1 = 
  {"STAT_DEF ", STAT_DEF_LINE, COMP_STRINGWISE, rdmc_amanda_STAT_DEF_1x1 };
const f2000_line_t USER_DEF_line_1x1 = 
  {"USER_DEF ", USER_DEF_LINE, COMP_STRINGWISE, rdmc_amanda_USER_DEF_1x1 };
const f2000_line_t KH_line_1x1 = 
  {"KH ", KH_LINE, COMP_STRINGWISE, rdmc_amanda_KH_1x1 };
const f2000_line_t OM_line_1x1 = 
  {"OM ", OM_LINE, COMP_STRINGWISE, rdmc_amanda_OM_1x1 };
const f2000_line_t KADC_line_1x1 = 
  {"KADC ", KADC_LINE, COMP_STRINGWISE, rdmc_amanda_KADC_1x1 };
const f2000_line_t KTDC_line_1x1 = 
  {"KTDC ", KTDC_LINE, COMP_STRINGWISE, rdmc_amanda_KTDC_1x1 };
const f2000_line_t KTOT_line_1x1 = 
  {"KTOT ", KTOT_LINE, COMP_STRINGWISE, rdmc_amanda_KTOT_1x1 };
const f2000_line_t KUTC_line_1x1 = 
  {"KUTC ", KUTC_LINE, COMP_STRINGWISE, rdmc_amanda_KUTC_1x1 };
const f2000_line_t TBEGIN_line_1x1 = 
  {"TBEGIN ", TBEGIN_LINE, COMP_STRINGWISE, rdmc_amanda_TBEGIN_1x1 };
const f2000_line_t EM_line_1x1 = 
  {"EM ", EM_LINE , COMP_STRINGWISE, rdmc_amanda_EM_1x1 };
const f2000_line_t TR_line_1x1 = 
  {"TR ", TR_LINE, COMP_STRINGWISE, rdmc_amanda_TR_1x1 };
const f2000_line_t CH_line_1x1 = 
  {"CH ", CH_LINE, COMP_STRINGWISE, rdmc_amanda_CH_1x1 };
const f2000_line_t HT_line_1x1 = 
  {"HT ", HT_LINE, COMP_STRINGWISE, rdmc_amanda_HT_1x1 };
const f2000_line_t US_line_1x1 = 
  {"US ", US_LINE, COMP_STRINGWISE, rdmc_amanda_US_1x1 };
const f2000_line_t STATUS_line_1x1 = 
  {"STATUS " , STATUS_LINE , COMP_STRINGWISE, rdmc_amanda_STATUS_1x1 };
const f2000_line_t FIT_block_1x1 = 
  {"FIT ", FIT_LINE , COMP_STRINGWISE, rdmc_amanda_fitblock_1x1 };
const f2000_line_t TRIG_block_1x1 = 
  {"TRIG ", TRIG_LINE, COMP_STRINGWISE, rdmc_amanda_trigblock_1x1 };
const f2000_line_t TRIG_line_1x1 = 
  {"TRIG ", TRIG_LINE , COMP_STRINGWISE, rdmc_amanda_TRIG_1x1 };
const f2000_line_t FIT_line_1x1 = 
  {"FIT ", FIT_LINE, COMP_STRINGWISE, rdmc_amanda_FIT_1x1 };
const f2000_line_t FRESULT_line_1x1 = 
  {"FRESULT " , FRESULT_LINE, COMP_STRINGWISE, rdmc_amanda_FRESULT_1x1 };
const f2000_line_t USES_line_1x1 = 
  {"USES ", USES_LINE, COMP_STRINGWISE, rdmc_amanda_USES_1x1 };
const f2000_line_t EE_line_1x1 = 
  {"EE", EE_LINE, COMP_STRINGWISE_NOTRAIL, rdmc_f2k_dummy_parser };
const f2000_line_t TEND_line_1x1 = 
  {"TEND ", TEND_LINE, COMP_STRINGWISE, rdmc_f2k_dummy_parser };
const f2000_line_t END_line_1x1 = 
  {"END", END_LINE, COMP_STRINGWISE_NOTRAIL, rdmc_f2k_dummy_parser };
/* this line is needed because a 1.1 header may be empty */
#if 0
const f2000_line_t DUMMY_line_1x1 = 
  {"", DUMMY_LINE, COMP_NOTHING, rdmc_f2k_dummy_parser };
#endif

const f2000_event_t f2000_preamble_1x1 =
{ RDMC_EVENT_HEADER_PREF,
  {&(V_line_1x1),NULL},
  {&(HI_line_1x1) , &(COMMENT_line_1x1),  NULL},
  {NULL},
  rdmc_rhd_f2k_1x1, 
  rdmc_f2k_dummy_event_writer
};

const f2000_event_t f2000_mhead_1x1 =
{ RDMC_EVENT_HEADER,
  {    NULL  },
  { &(ARRAY_line_1x1),&(KH_line_1x1),&(FIT_DEF_line_1x1),&(TRIG_DEF_line_1x1)
   ,&(TRIG_PAR_line_1x1),&(STAT_DEF_line_1x1),&(USER_DEF_line_1x1)
   ,&(KADC_line_1x1),&(KTDC_line_1x1),&(KTOT_line_1x1),&(KUTC_line_1x1)
   ,&(OM_line_1x1),&(COMMENT_line_1x1),&(TBEGIN_line_1x1)
   , NULL
  },
  {NULL},
  rdmc_mhead_f2k_1x1,
  rdmc_whead_f2k_1x1_2
};


const f2000_event_t f2000_mevt_1x1 =
{ RDMC_EVENT_MUON,
  { &(EM_line_1x1),
    NULL
  },
  {
    &(TR_line_1x1),
    &(HT_line_1x1),  &(CH_line_1x1),    
    &(STATUS_line_1x1),   &(US_line_1x1),
    &(FIT_block_1x1),  &(FRESULT_line_1x1),
    &(TRIG_block_1x1), &(USES_line_1x1),
    &(COMMENT_line_1x1)
    ,  NULL
  },
  { &(EE_line_1x1)
    , NULL},
  rdmc_mevt_f2k_1x1,
  rdmc_wevt_f2k_1x1_2
};

const f2000_event_t f2000_mfoot_1x1 =
{ RDMC_EVENT_FOOT,
  { &(TEND_line_1x1),
    NULL
  },
  {    &(COMMENT_line_1x1)
    ,  NULL
  },
  { &(END_line_1x1)
    , NULL},
  rdmc_mfoot_f2k_1x1,
  rdmc_wfoot_f2k_1x1_2
};

const f2000_event_t  * f2000_events_1x1[] 
  = { 
    &f2000_mevt_1x1,
    NULL 
  };


/****************************************************************************
 * rhd_amanda() reads the format relevant informations for amanda-like
 *          formats - just the Comments and history lines. 
 ****************************************************************************/
int rdmc_rhd_f2k_1x1(mcfile *fp, array *a, mevt *e){
  int r=RDMC_IO_OK;
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  while( f2k_buff->iline <   f2k_buff->lines ){
    switch(f2k_buff->type_def[f2k_buff->iline]->line_id){
    case COMMENT_LINE:
      rdmc_append_comment(&(fp->comment),f2k_buff->line_pt[f2k_buff->iline]);
      break;
    default: /* HI line, v-line is dummy */
      r=f2k_buff->type_def[f2k_buff->iline]->parser(fp,a,e,NULL);
      break;
    }
    if(r != RDMC_IO_OK){
      rdmc_f2k_errorlog(fp);
      return r;
    }
    f2k_buff->iline++;
  }
  return r;
  
} /* function rhd_amanda() */


/****************************************************************************
 * The function mhead reads the header of a amanda like file
 ****************************************************************************/
int rdmc_mhead_f2k_1x1(mcfile *fp, array *ar, mevt *e){
  int ret=RDMC_IO_OK;  /* the return value */
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;

  rdmc_init_array(ar);                           /* reset the array */

  while (f2k_buff->iline < f2k_buff->lines){
    switch(f2k_buff->type_def[f2k_buff->iline]->line_id){
    case COMMENT_LINE:
      rdmc_append_comment(&(ar->comment),f2k_buff->line_pt[f2k_buff->iline]);
      break;
    default: /* any other line */
      ret = f2k_buff->type_def[f2k_buff->iline]->parser(fp,ar,e,NULL);
      break;
    }
    if(ret != RDMC_IO_OK){
      rdmc_f2k_errorlog(fp);
      return ret;
    }
    f2k_buff->iline++;
  }
  return ret;
} /* function rhd_amanda() */

/****************************************************************************
 * The function mhead reads the header of a amanda like file
 ****************************************************************************/
int rdmc_mevt_f2k_1x1(mcfile *fp, array *ar, mevt *e){
  int ret=RDMC_IO_OK;  /* the return value */
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;

  /*************************************************/
  while (f2k_buff->iline < f2k_buff->lines){
    switch(f2k_buff->type_def[f2k_buff->iline]->line_id){
    case COMMENT_LINE:
      rdmc_append_comment(&(e->comment),f2k_buff->line_pt[f2k_buff->iline]);
      break;
    case EE_LINE:
      e->ee = 1;
      ret = RDMC_IO_OK;
      break;
    default: /* any other line */
      ret = f2k_buff->type_def[f2k_buff->iline]->parser(fp,ar,e,NULL);
      break;
    }
    if(ret != RDMC_IO_OK){
      rdmc_f2k_errorlog(fp);
      return ret;
    }
    f2k_buff->iline++;
  }
  e->nch=rdmc_count_nch(e);         /* calc the number of hit channels */
  rdmc_fill_mhit_str(e,ar);
  e->nstr = rdmc_count_nstr(e);
  if (fp->sloppy)
    rdmc_repair_mhit_id(e);

  return ret;
} /* function revt_amanda() */


/****************************************************************************
 * The function mhead reads the header of a amanda like file
 ****************************************************************************/
int rdmc_mfoot_f2k_1x1(mcfile *fp, array *ar, mevt *e){
  int ret=RDMC_IO_OK;  /* the return value */
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;

  /*************************************************/

  while (f2k_buff->iline < f2k_buff->lines){
    switch(f2k_buff->type_def[f2k_buff->iline]->line_id){
    case COMMENT_LINE:
      rdmc_append_comment(&(e->comment),f2k_buff->line_pt[f2k_buff->iline]);
      break;
    default: /* any other line */
      ret = f2k_buff->type_def[f2k_buff->iline]->parser(fp,ar,e,NULL);
      break;
    }
    if(ret != RDMC_IO_OK){
      rdmc_f2k_errorlog(fp);
      return ret;
    }
    f2k_buff->iline++;
  }
  return ret;
} /* function revt_amanda() */



/****************************************************************************
 * Read a history line
 ****************************************************************************/
int rdmc_amanda_HI_1x1( mcfile *fp , array *a, 
			       mevt *e, void *tmp)
{
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char *p=s+strlen("HI ");

  /* now append it */
  if (*p == '\0')  /* nothing there, -> add a '\n' */
    rdmc_append_comment(&(fp->creator),"\n");
  else
    rdmc_append_comment(&(fp->creator),p);

  return RDMC_IO_OK;

} /* rdmc_amanda_HI() */


/****************************************************************************
 * read the various lines
 ***************************************************************************/

int rdmc_amanda_ARRAY_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  args = rdmc_f2k_tokenize(s,&nargs);
  if (nargs != 7) return RDMC_ILF;

  /* arg[0] is the token */
  a->id = rdmc_amanda_idet(args[1]);
  a->longitude = rdmc_amanda_strtof(args[2],RDMC_LONGI_NA);
  a->lattitude = rdmc_amanda_strtof(args[3],RDMC_LATTI_NA);
  a->depth = rdmc_amanda_strtof(args[4],RDMC_DEPTH_NA);
  a->nstr =  rdmc_amanda_strtoi(args[5],RDMC_NA);
  a->nch =   rdmc_amanda_strtoi(args[6],RDMC_NA);

  if (a->nch > RDMC_MAXCHANNELS)
    return RDMC_TOO_MANY_CHANNELS;

  a->is_calib.geo=1; /* set geo cal flag */

  return RDMC_IO_OK;
} /* rdmc_amanda_ARRAY() */

int rdmc_amanda_TRIG_DEF_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  array_hdef_t def_trig;
  int i;
  int i_trigger;

  args = rdmc_f2k_tokenize(s,&nargs);
  /* token + id + token arguments */
  if (nargs < 3) return RDMC_ILF;
  
  rdmc_init_array_hdef(&def_trig);

  /* now check if this id is existing r */
  i_trigger = rdmc_get_hdef_tag(a->def_trig, a->n_trigger,args[2] );
  if (i_trigger >= 0)                      /* this is not unique */
    return RDMC_INCONSISTENT_GEOMETRY;


  def_trig.id = atoi(args[1]) ;
  strcpy(def_trig.tag,args[2] );
  
  if ( (def_trig.nwords = nargs - 3) >  RDMC_MAXTOKEN_PER_LINE)
    return RDMC_LINE_NOT_PARSED;

  for (i=0 ; i < def_trig.nwords ; i++){
    strcpy(def_trig.words[i],args[i+3]);
  }
  rdmc_add_trigger_def(a,&def_trig,a->n_trigger);

  return RDMC_IO_OK;
} /* rdmc_amanda_TRIG_DEF() */

int rdmc_amanda_TRIG_PAR_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0,npars;
  int i,i_trigger;

  args = rdmc_f2k_tokenize(s,&nargs);
  /* token + id + token arguments */
  if (nargs < 2) return RDMC_ILF;
  
  /* now check if this id is existing r */
  i_trigger = rdmc_get_hdef_tag(a->def_trig, a->n_trigger,args[1] );
  if (i_trigger < 0)
    return RDMC_INCONSISTENT_GEOMETRY;

  if (a->def_trig[i_trigger].npars > 0)
    return RDMC_INCONSISTENT_GEOMETRY;
  else
    npars = nargs - 2; 
  
  if (npars  >  RDMC_MAXTOKEN_PER_LINE){
    return RDMC_LINE_NOT_PARSED;
  }

  a->def_trig[i_trigger].npars=npars;

  for (i=0 ; i < npars ; i++){
    strcpy(a->def_trig[i_trigger].pars[i],args[i+2]);
  }

  return RDMC_IO_OK;
} /* rdmc_amanda_TRIG_PAR() */


int rdmc_amanda_FIT_DEF_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  array_hdef_t def_fit;
  int i;
  int i_fit;

  args = rdmc_f2k_tokenize(s,&nargs);
  /* token + id + token arguments */
  if (nargs < 3) return RDMC_ILF;
  
  rdmc_init_array_hdef(&def_fit);

  def_fit.id = atoi(args[1]) ;
  strcpy(def_fit.tag,args[2] );

  /* now check if this id is existing r */
  i_fit = rdmc_get_hdef_id(a->def_fit, a->n_fit, def_fit.id );
  if (i_fit >= 0)                      /* this is not unique */
    return RDMC_INCONSISTENT_GEOMETRY;
  
  if ( (def_fit.nwords = nargs - 3) >  RDMC_MAXTOKEN_PER_LINE)
    return RDMC_LINE_NOT_PARSED;

  for (i=0 ; i < def_fit.nwords ; i++){
    strcpy(def_fit.words[i],args[i+3]);
  }
  rdmc_add_fit_def(a,&def_fit,a->n_fit);

  return RDMC_IO_OK;
} /* rdmc_amanda_FIT_DEF() */

int rdmc_amanda_STAT_DEF_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  array_hdef_t def_stat;
  int i;
  int i_stat;

  args = rdmc_f2k_tokenize(s,&nargs);
  /* token + tag + token arguments */
  if (nargs < 2) return RDMC_ILF;
  
  rdmc_init_array_hdef(&def_stat);

  def_stat.id = a->n_stat;
  strcpy(def_stat.tag,args[1] );

  /* now check if this id is existing r */
  i_stat = rdmc_get_hdef_id(a->def_stat, a->n_stat, def_stat.id );
  if (i_stat >= 0)                      /* this is not unique */
    return RDMC_INCONSISTENT_GEOMETRY;
  
  if ( (def_stat.nwords = nargs - 2) >  RDMC_MAXTOKEN_PER_LINE)
    return RDMC_LINE_NOT_PARSED;

  for (i=0 ; i < def_stat.nwords ; i++){
    strcpy(def_stat.words[i],args[i+2]);
  }
  rdmc_add_stat_def(a,&def_stat,a->n_stat);

  return RDMC_IO_OK;
} /* rdmc_amanda_STAT_DEF() */


int rdmc_amanda_USER_DEF_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  array_hdef_t def_user;
  int i;
  int i_user;

  args = rdmc_f2k_tokenize(s,&nargs);
  /* token + id + token arguments */
  if (nargs < 2) return RDMC_ILF;
  
  rdmc_init_array_hdef(&def_user);

  strcpy(def_user.tag,args[1] );

  /* now check if this id is existing r */
  i_user = rdmc_get_hdef_tag(a->def_user, a->n_user, def_user.tag );
  if (i_user >= 0)                      /* this is not unique */
    return RDMC_INCONSISTENT_GEOMETRY;
  
  if ( (def_user.nwords = nargs - 2) >  RDMC_MAXTOKEN_PER_LINE)
    return RDMC_LINE_NOT_PARSED;

  for (i=0 ; i < def_user.nwords ; i++){
    strcpy(def_user.words[i],args[i+2]);
  }
  rdmc_add_user_def(a,&def_user,a->n_user);

  return RDMC_IO_OK;
} /* rdmc_amanda_USER_DEF() */

int rdmc_amanda_KH_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){

  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  int i;

  args = rdmc_f2k_tokenize(s,&nargs);
  /* token + id + token arguments */
  if (nargs < 1) return RDMC_ILF;

  for (i = 1 ; i < nargs ; i++ ){
    if ( !strcmp(args[i],"ADC")  ){
      a->is_calib.adc=1;
    }
    else if ( !strcmp(args[i],"TDC")  ){
      a->is_calib.tdc=1;
    }
    else if ( !strcmp(args[i],"UTC")  ){
      a->is_calib.utc=1;
    }
#if 1 /* this is not part of f2000 1.1 but who cares */
    else  if ( !strcmp(args[i],"TOT")  ){
      a->is_calib.tot=1;
    }
#endif    
    else   { 
      if (!(fp->sloppy)) /* no sloppy mode -> return with error */
	return RDMC_LINE_NOT_PARSED;
      else
	continue;
    }
  } /* for */

  return RDMC_IO_OK;
} /* rdmc_amanda_KH() */


int rdmc_amanda_OM_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){

  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  int nr;
  args = rdmc_f2k_tokenize(s,&nargs);
  /* token + id + token arguments */
  if (nargs < 8) return RDMC_ILF;

  nr = atoi(args[1]) - 1;
  if ((nr < 0) || (nr >= a->nch)) return RDMC_INCONSISTENT_GEOMETRY;
  if ( ( a->str[nr] = atoi(args[3])) > a->nstr)
    return RDMC_INCONSISTENT_GEOMETRY;
  a->clust[nr] = atoi(args[2]) - 1;
  if (isnan(a->x[nr] = atof(args[4]))) return RDMC_ILF;
  if (isnan(a->y[nr] = atof(args[5]))) return RDMC_ILF;
  if (isnan(a->z[nr] = atof(args[6]))) return RDMC_ILF;
  a->costh[nr] = rdmc_amanda_fori(args[7]);
  a->type[nr]  = rdmc_amanda_iomid(
				    (nargs > 8) ? args[8] : "unknown" 
				    ,a->id);
  a->serial[nr] = rdmc_amanda_strtoi(
				     (nargs > 9) ? args[9] : "na"
				     ,RDMC_NA);
  a->sensit[nr] = rdmc_amanda_strtof(
				     (nargs > 10) ? args[10] : "1.0"
				     , 1.0);

  if( nargs > 11 ) 
    a->thresh[nr] = rdmc_amanda_strtof(args[11],RDMC_SMALL);
  else{
    char tmp[RDMC_MAXTOKENLENGTH];
    sprintf(tmp,"%g",RDMC_SMALL);
    a->thresh[nr] = rdmc_amanda_strtof(tmp,RDMC_SMALL);
  }

  /* finally at least one channel is geocalibratet -> so set the flag */
  a->is_calib.geo=1;
  return RDMC_IO_OK;
} /* rdmc_amanda_KH() */

/****************************************************************************
 * read the ADC calibration (only nr ped beta lin)
 ****************************************************************************/
int rdmc_amanda_KADC_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  int iom;
  args = rdmc_f2k_tokenize(s,&nargs);
  /* token + id + token arguments */
  if ((nargs < 4) || (nargs > 5)) return RDMC_ILF;

  iom = atoi(args[1]) - 1;
  if ((iom < 0) || (iom >= a->nch)) return RDMC_INCONSISTENT_GEOMETRY;
  a->cal[iom].ped = rdmc_amanda_strtof(args[2],0.0);
  a->cal[iom].beta_a = rdmc_amanda_strtof(args[3],0.0);
  if (nargs == 4)
      a->cal[iom].kappa = rdmc_amanda_strtof("?",0.0);
    else
      a->cal[iom].kappa = rdmc_amanda_strtof(args[4],0.0);
  a->cal[iom].flag |= RDMC_CALIB_ADC;

  /* set the kh flag, even if it was not set before */
  a->is_calib.adc = 1;

  return 0;
} /* rdmc_amanda_KADC() */

/****************************************************************************
 * Read the TDC calibration (only nr beta shift alpha now)
 ****************************************************************************/
int rdmc_amanda_KTDC_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  int iom;
  args = rdmc_f2k_tokenize(s,&nargs);
  /* token + id + token arguments */
  if (nargs != 5) return RDMC_ILF;
  iom = atoi(args[1]) - 1;  
  if ((iom < 0) || (iom >= a->nch)) return RDMC_INCONSISTENT_GEOMETRY;

  a->cal[iom].beta_t =  rdmc_amanda_strtof(args[2],0.0);
  a->cal[iom].t_0 =  rdmc_amanda_strtof(args[3],0.0);
  a->cal[iom].alpha_t =  rdmc_amanda_strtof(args[4],0.0);
  a->cal[iom].flag |= RDMC_CALIB_TDC;

  return 0;

} /* rdmc_amanda_KTDC() */
/****************************************************************************
 * Read the UTC calibration (only nr beta shift alpha now)
 ****************************************************************************/
int rdmc_amanda_KUTC_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){

  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  char *t;

  args = rdmc_f2k_tokenize(s,&nargs);
  /* token + id + token arguments */
  if (nargs != 3) return RDMC_ILF;

  strcpy(a->cal_utc.utc_src,args[1]);

  /* now check sec.nsec, append trailing 0's, if necessary */

  t = strchr(args[2],'.');
  if (t == NULL){
    a->cal_utc.secs =  rdmc_amanda_strtoi(args[2],0);
    a->cal_utc.nsecs =  rdmc_amanda_strtoi("?",0);
  }else{
    *t='\0';
    a->cal_utc.secs =  rdmc_amanda_strtoi(args[2],0);
    a->cal_utc.nsecs =  rdmc_amanda_strtoi(t,0);
  }
  
  return RDMC_IO_OK;
} /* rdmc_amanda_UTC() */

int rdmc_amanda_KTOT_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  int iom;
  args = rdmc_f2k_tokenize(s,&nargs);
  /* token + id + token arguments */
  if ((nargs < 4) || (nargs > 5)) return RDMC_ILF;
  iom = atoi(args[1]) - 1;  

  if ((iom < 0) || (iom >= a->nch)) return RDMC_INCONSISTENT_GEOMETRY;

  a->cal[iom].ped_tot = rdmc_amanda_strtof(args[2],0.0);
  a->cal[iom].beta_tot = rdmc_amanda_strtof(args[3],1.0);
  if (nargs == 5)
    a->cal[iom].kappa_tot = rdmc_amanda_strtof(args[4],0.0);
  else
    a->cal[iom].kappa_tot = rdmc_amanda_strtof("?",0.0);
  a->cal[iom].flag |= RDMC_CALIB_TOT;

  /* set the kh flag, even if it was not set before */
  a->is_calib.tot = 1;

  return 0;
} /* rdmc_amanda_KTOT() */

int rdmc_amanda_TBEGIN_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  int isec,iday,iyear;
  char *csec="?",*cday="?",*cyear="?";

  double fsec;
  struct tm newyear;

  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);
  switch (nargs){
  case 4:
    cyear = args[1];
    cday = args[2];
    csec = args[3];
    break;
  case 3:
    cyear = args[1];
    cday = args[2];
    break;
  case 2:
    cyear = args[1];
    break;
  case 1:
    break;
  default:
   return RDMC_ILF; 
  }
  iyear = rdmc_amanda_strtoi(cyear, 1970);
  iday = rdmc_amanda_strtoi(cday, 1);
  fsec = rdmc_amanda_strtoi(csec, 0);

#if 0
  isec = rdmc_nint(fsec);
#else
  isec = floor(fsec);
#endif
  if (iyear <  0) iyear = 1970;
  else if ((iyear < 70)) iyear += 100;
  if (iyear > 1900) iyear -= 1900;
  
  isec=isec%SEC_PER_DAY ;  /* seconds of begin of day only is mod(sec,86400)*/

  /* get time_t for beginning of that year */
  newyear.tm_sec=newyear.tm_min=newyear.tm_hour
    =newyear.tm_mday=newyear.tm_mon=newyear.tm_year
    =newyear.tm_wday=newyear.tm_yday=newyear.tm_isdst=0;

  newyear.tm_year=iyear;
  newyear.tm_mday = 1;

  a->tbegin = mktime(&newyear) + SEC_PER_DAY*(iday-1) + isec -timezone ;
  return RDMC_IO_OK;
}

/****************************************************************************
 * read an Event header line
 ****************************************************************************/

int rdmc_amanda_EM_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  int iday,iyear;
  char *c_enr="?",*c_nrun="?",*c_iyear="?",*c_iday="?";
  char *c_time="0.000000000",*c_t_offset="?";

  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);


  switch (nargs){
  case 1:
    break;
  case 2:
    c_enr=args[1];
    break;
  case 3:
    c_enr=args[1];
    c_nrun=args[2];
    break;
  case 4:
    c_enr=args[1];
    c_nrun=args[2];
    c_iyear=args[3];
    break;
  case 5:
    c_enr=args[1];
    c_nrun=args[2];
    c_iyear=args[3];
    c_iday=args[4];
    break;
  case 6:
    c_enr=args[1];
    c_nrun=args[2];
    c_iyear=args[3];
    c_iday=args[4];
    c_time=args[5];
    break;
  case 7:
    c_enr=args[1];
    c_nrun=args[2];
    c_iyear=args[3];
    c_iday=args[4];
    c_time=args[5];
    c_t_offset=args[6];
    break;
  default:
   return RDMC_ILF; 
  }
  e->enr = rdmc_amanda_strtoi(c_enr,RDMC_NA);
  e->nrun = rdmc_amanda_strtoi(c_nrun,RDMC_NA);
  iyear = rdmc_amanda_strtoi(c_iyear ,  1970);
  iday = rdmc_amanda_strtoi(c_iday , 1);
  rdmc_amanda_strtimetoi(c_time, &(e->secs), &(e->nsecs));
  e->t_offset = rdmc_amanda_strtof(c_t_offset,0.0);


  /* now patch the year in case of ssome f2k dialects and y2k-type problems */
  iyear = rdmc_f2k_y2k(iyear);
  e->mjd = rdmc_gps_to_mjd(iyear,iday); /* convert GPS date into mjd */
  
  return RDMC_IO_OK;
} /* rdmc_amanda_EM() */

int rdmc_amanda_TR_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  float azi;
  static mtrack gen;

  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);

  if ( nargs != 12) 
    return RDMC_ILF;

  /* get new mem for this track and init index */
  rdmc_init_mtrack(&gen);

  /* now fill it */
  gen.tag = rdmc_amanda_strtoi(args[1], RDMC_NA);;
  gen.parent = rdmc_amanda_strtoi(args[2], RDMC_PARENT_NA);
  gen.id = rdmc_amanda_ipartid(args[3]);

 gen.x = rdmc_amanda_strtof(args[4],RDMC_SPECIAL_NA);
  gen.y = rdmc_amanda_strtof(args[5],RDMC_SPECIAL_NA);
  gen.z = rdmc_amanda_strtof(args[6],RDMC_SPECIAL_NA);

  gen.costh = cos(rdmc_amanda_strtof(args[7],0.)  *PID180); 
  azi = rdmc_amanda_strtof(args[8],0.);
  if ((azi = fmod( azi, 360.)) <0. ) azi += 360;
  gen.phi =  azi * PID180;
  rdmc_tau_tr(&gen); /* calc direction cosinuus from phi and costh */

  gen.length = rdmc_amanda_strtof(args[9],RDMC_NA);
  gen.e = 1000.0 * rdmc_amanda_strtof(args[10],RDMC_NA);/* GeVtoMeV (rdmc)*/
  gen.t = rdmc_amanda_strtof(args[11],RDMC_TDC_NA);
  gen.nmuon = 0;

  rdmc_add_gen(e,&gen,e->ntrack);

#if 0 /* no dynamic allocated stuff */
  rdmc_clear_mtrack(&gen);
#endif

  return RDMC_IO_OK;
} /* rdmc_amanda_TR() */

int rdmc_amanda_CH_1x1(  mcfile *fp , array *a, 
			       mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  int nhits=0;
  mhit h;
  int om;
  int j;

  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);

  if ( nargs < 2) 
    return RDMC_ILF;
  else if ( ((nargs-2)%5) != 0)
    return RDMC_ILF;
  else if ( (nhits = (nargs-2)/5) < 1)
    return RDMC_ILF;
     
  
     /* get new mem for this track and init index */
  rdmc_init_mhit(&h);
  om=atoi(args[1]); /* cp ch token */
  
  for( j=0 ; j < nhits ; j++ ){
    h.str=0;
    h.ch = om-1;
    if( j== 0)
      h.amp = rdmc_amanda_strtof(args[1+j*5+1], RDMC_NA);
    else
      h.amp=RDMC_REPEAT_CH;
    h.id = atoi(args[1+j*5+2]);
    h.mt = rdmc_amanda_strtoi(args[1+j*5+3],RDMC_PARENT_NA);
    h.ma = h.mt;
    h.t = rdmc_amanda_strtof(args[1+j*5+4], RDMC_TDC_NA);
    h.tot = rdmc_amanda_strtof(args[1+j*5+5], RDMC_NA);

    rdmc_add_mhit(e,&h,e->nhits);

    rdmc_clear_mhit(&h);

  } /* for nhits */

  return RDMC_IO_OK;
} /* rdmc_amanda_CH() */


int rdmc_amanda_HT_1x1( mcfile *fp , array *a, 
			mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  mhit h;

  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);

  if ( nargs != 7) 
    return RDMC_ILF;

  rdmc_init_mhit(&h);
  h.str = 0; /* this is filled later at the end of the event read */
  h.ch = atoi(args[1])-1;
  h.amp = rdmc_amanda_strtof(args[2] , RDMC_NA);
  h.id = rdmc_amanda_strtoi(args[3] ,RDMC_NA);
  h.ma = h.mt = rdmc_amanda_strtoi(args[4] ,RDMC_PARENT_NA);
  h.t = rdmc_amanda_strtof(args[5] , RDMC_TDC_NA);
  h.tot = rdmc_amanda_strtof(args[6], RDMC_NA);

  rdmc_add_mhit(e,&h,e->nhits);
#if 1 /* no dynamic allocated stuff */
  rdmc_clear_mhit(&h);
#endif

  return RDMC_IO_OK;
} /* rdmc_amanda_HT() */



int rdmc_amanda_trigblock_1x1( mcfile *fp, array *ar, 
				     mevt *e, void *tmp){
  int r,i; /*ret value */
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  mevt_special_t trig;
  mevt_usesblock_t uses;

  rdmc_init_mevt_special(&trig,0);
  rdmc_init_usesblock(&uses);

  /*  this should parse TRIG  */
  if ( f2k_buff->type_def[f2k_buff->iline]->line_id == TRIG_LINE ){
    r = TRIG_line_1x1.parser(fp,ar,e,&trig);
    if (r != RDMC_IO_OK)
      return RDMC_LINE_NOT_PARSED;
  }else 
    return RDMC_LINE_NOT_PARSED;

  uses.id= e->ntrig;
  while( ++(f2k_buff->iline) < f2k_buff->lines ){/* scan while it is possible to parse */
    switch( f2k_buff->type_def[f2k_buff->iline]->line_id ){
    default: 
      r = RDMC_EVENT_NOT_RECOGNIZED;
      --(f2k_buff->iline);
      break;
    case USES_LINE:
      r = USES_line_1x1.parser(fp,ar,e,&uses);
      break;
    case COMMENT_LINE:
      rdmc_append_comment(&(e->comment),f2k_buff->line_pt[f2k_buff->iline]);
      break;
    }
    if (r != RDMC_IO_OK){
      break;
    }
  } /* while */
  
  if ((r == RDMC_IO_OK) || (r ==  RDMC_EVENT_NOT_RECOGNIZED)){
    rdmc_add_trigger(e, &trig, e->ntrig, trig.id);

    for ( i=0 ; i<uses.nu ; i++) /* this is slow but will change anyway */
      rdmc_add_trig_uses(e,&(uses.u[i]),e->ntrig_uses);
    rdmc_clear_mevt_special(&trig,0);
    rdmc_clear_usesblock(&uses);
    return RDMC_IO_OK;
  } else{
    rdmc_f2k_errorlog(fp);
    return r;
  }
}


int rdmc_amanda_TRIG_1x1(mcfile *fp, array *a, mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  mevt_special_t *trig = tmp;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;

  int itrig,itrig_def,icount,itrig_defid;
  char *tag_s; /* index of parent track  */
  char *tp;

  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);
  if ( nargs < 2) 
    return RDMC_ILF;
  tag_s = args[1];

  /* now get the trigger number */
  itrig_def = rdmc_get_hdef_tag(a->def_trig ,a->n_trigger,tag_s);
  if (itrig_def == RDMC_NA ){
    /* special patch for old amanda files */
    if (strstr(tag_s, "unknown") != NULL) /*this is old rdmc ignore the line */
      return RDMC_IO_OK;
    else
      return RDMC_HEADERTAG_NOTFOUND;
  }
  itrig_defid = a->def_trig[itrig_def].id;

  rdmc_clear_mevt_special(trig,a->def_trig[itrig_def].nwords);
  itrig = e->ntrig;

  trig->id = itrig_def;
  trig->nval = 0;
  icount=0;
  while ( (icount < a->def_trig[itrig_def].nwords) 
	  &&  (icount < RDMC_MAXTOKEN_PER_LINE )
	  &&   ((icount+2) < nargs) ){
    tp = args[icount+2];
    trig->val[icount] =  rdmc_amanda_sptof(tp);
    icount++;
    trig->nval=icount;
  }
  return RDMC_IO_OK;
}


int rdmc_amanda_fitblock_1x1( mcfile *fp, array *ar, 
				     mevt *e, void *tmp){
  int i,r; /*ret value */
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;

  mtrack rec;
  mevt_special_t fresult;
  mevt_usesblock_t uses;

  rdmc_init_mtrack(&rec);  
  rdmc_init_mevt_special(&fresult,0);
  rdmc_init_usesblock(&uses);

  /*  this should parse FIT  */
  if ( f2k_buff->type_def[f2k_buff->iline]->line_id != FIT_LINE )
    return RDMC_LINE_NOT_PARSED;

  r = FIT_line_1x1.parser(fp,ar,e,&rec);
  if (r != RDMC_IO_OK){
    return RDMC_LINE_NOT_PARSED;
  }

  uses.id= e->nfit;

  while( ++(f2k_buff->iline) < f2k_buff->lines ){ /* scan while it is possible to parse */
    switch( f2k_buff->type_def[f2k_buff->iline]->line_id ){
    default: 
      r = RDMC_EVENT_NOT_RECOGNIZED;
      --(f2k_buff->iline);
      break;
    case FRESULT_LINE:
      r = FRESULT_line_1x1.parser(fp,ar,e,&fresult);
      break;
    case USES_LINE:
      r = USES_line_1x1.parser(fp,ar,e,&uses);
      break;
    case COMMENT_LINE:
      rdmc_append_comment(&(e->comment),f2k_buff->line_pt[f2k_buff->iline]);
      break;
    }
    if (r != RDMC_IO_OK){
      break;
    }
  } /* while */
  
  if ((r == RDMC_IO_OK) || (r ==  RDMC_EVENT_NOT_RECOGNIZED)){
    rdmc_add_fit(e,&rec,&fresult,e->nfit);
    for (i=0 ; i<uses.nu ; i++) /* this is slow but will change anyway */
      rdmc_add_fit_uses(e,&(uses.u[i]),e->nfit_uses);
    rdmc_clear_mtrack(&rec);  
    rdmc_clear_mevt_special(&fresult,0);
    rdmc_clear_usesblock(&uses);
    return RDMC_IO_OK;
  }else{
    rdmc_f2k_errorlog(fp);
    return r;
  }
}


int rdmc_amanda_US_1x1(mcfile *fp, array *a, mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;

  int ius_def,icount,ius_defid;
  char *tag_s; /* index of parent track  */
  char *tp;
  static mevt_special_t us;
   

  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);
  if ( nargs < 2) 
    return RDMC_ILF;
  tag_s = args[1];

  /* now get the user number */
  ius_def = rdmc_get_hdef_tag(a->def_user ,a->n_user,tag_s);
  if (ius_def == RDMC_NA ){
    return RDMC_HEADERTAG_NOTFOUND;
  }
  ius_defid = a->def_user[ius_def].id;

  rdmc_init_mevt_special(&us,a->def_user[ius_def].nwords);

  us.id = ius_def;
  us.nval = 0;
  icount=0;
  while ( (icount < a->def_user[ius_def].nwords) 
	  &&  (icount < RDMC_MAXTOKEN_PER_LINE )
	  &&   ((icount+2) < nargs) ){
    tp = args[icount+2];
    us.val[icount] =  rdmc_amanda_sptof(tp);
    icount++;
    us.nval=icount;
  }
  rdmc_add_user(e,&us,e->nuser);

  rdmc_clear_mevt_special(&us,0);
  return RDMC_IO_OK;
}
int rdmc_amanda_STATUS_1x1(mcfile *fp, array *a, mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;

  int istatus_def,icount,istatus_defid;
  char *tag_s; /* index of parent track  */
  char *tp;
  static mevt_special_t status;
   
  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);
  if ( nargs < 2) 
    return RDMC_ILF;
  tag_s = args[1];

  /* now get the user number */
  istatus_def = rdmc_get_hdef_tag(a->def_stat ,a->n_stat,tag_s);
  if (istatus_def == RDMC_NA ){
    return RDMC_HEADERTAG_NOTFOUND;
  }
  istatus_defid = a->def_stat[istatus_def].id;

  rdmc_init_mevt_special(&status,a->def_stat[istatus_def].nwords);

  status.id = istatus_def;
  status.nval = 0;
  icount=0;
  while ( (icount < a->def_stat[istatus_def].nwords) 
	  &&  (icount < RDMC_MAXTOKEN_PER_LINE )
	  &&   ((icount+2) < nargs) ){
    tp = args[icount+2];
    status.val[icount] =  rdmc_amanda_sptof(tp);
    icount++;
    status.nval=icount;
  }
  rdmc_add_status(e,&status,e->nstat);

  rdmc_clear_mevt_special(&status,0);
  return RDMC_IO_OK;
}

int rdmc_amanda_FIT_1x1(mcfile *fp, array *a, mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  mtrack *fit = tmp;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  float azi;
  char *length_s,*energy_s;

  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);

  switch(nargs){
  case 9:
    length_s = energy_s="?";
  case 10:
    length_s = args[9];
    energy_s = "?";
    break;
  case 11:
    length_s = args[9]; 
    energy_s = args[10];
    break;
  default:
    return RDMC_ILF;
  }
  /* get new mem for this track and init index */
  rdmc_clear_mtrack(fit);

  /* now fill it */
  fit->tag = rdmc_amanda_strtoi(args[1], RDMC_NA);;
  fit->id = rdmc_amanda_ipartid(args[2]);
  fit->x = rdmc_amanda_strtof(args[3],RDMC_SPECIAL_NA);
  fit->y = rdmc_amanda_strtof(args[4],RDMC_SPECIAL_NA);
  fit->z = rdmc_amanda_strtof(args[5],RDMC_SPECIAL_NA);


  fit->costh = cos(rdmc_amanda_strtof(args[6],0.)  *PID180); 
  azi = rdmc_amanda_strtof(args[7],0.);
  if ((azi = fmod( azi, 360.)) <0. ) azi += 360;
  fit->phi =  azi * PID180;
  rdmc_tau_tr(fit); /* calc direction cosinus from phi and costh */

  fit->t = rdmc_amanda_strtof(args[8],RDMC_TDC_NA);
  fit->length = rdmc_amanda_strtof(length_s,RDMC_BIG);
  fit->e = 1000.0 * rdmc_amanda_strtof(energy_s,0.);/* GeVtoMeV (rdmc)*/

  return RDMC_IO_OK;
}

int rdmc_amanda_FRESULT_1x1(mcfile *fp, array *a, mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  mevt_special_t *fresult = tmp;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  int ifit_id;

  int ifit_def,icount,ifit_defid;
  char *tp;

  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);
  if ( nargs < 2) 
    return RDMC_ILF;
  /* now get the trigger number */
  ifit_id=atoi(args[1]);
  ifit_def = rdmc_get_hdef_id(a->def_fit ,a->n_fit,ifit_id);
  if (ifit_def == RDMC_NA ){
    return RDMC_HEADERTAG_NOTFOUND;
  }
  ifit_defid = a->def_fit[ifit_def].id;

  rdmc_clear_mevt_special(fresult,a->def_fit[ifit_def].nwords);

  fresult->id = ifit_def;
  fresult->nval = 0;
  icount=0;
  while ( (icount < a->def_fit[ifit_def].nwords) 
	  &&  (icount < RDMC_MAXTOKEN_PER_LINE )
	  &&   ((icount+2) < nargs) ){
    tp = args[icount+2];
    fresult->val[icount] =  rdmc_amanda_sptof(tp);
    icount++;
    fresult->nval=icount;
  }
  return RDMC_IO_OK;
}

int rdmc_amanda_USES_1x1(mcfile *fp, array *a, mevt *e, void *tmp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  char *s=f2k_buff->line_pt[f2k_buff->iline];
  char **args=NULL;
  int nargs=0;
  mevt_usesblock_t  *uses = tmp; /* uses HAS to be initilized and uses.id 
				    has to be set ! */
  int n_allocated=0;
  int i;
  char *c;

  /* start parsing*/
  args = rdmc_f2k_tokenize(s,&nargs);
  if ( nargs < 1) 
    return RDMC_ILF;

  /* nasty but cannot be done differently */
  n_allocated = uses->nu + RDMC_MAXTOKEN_PER_LINE;
  uses->u = realloc(uses->u ,n_allocated * sizeof(mevt_uses_t) );

  for ( i=1 ; i < nargs ; i++){
    if (uses->nu >= n_allocated-1){
      n_allocated += RDMC_MAXTOKEN_PER_LINE;
      uses->u = realloc( uses->u , n_allocated*sizeof(mevt_uses_t));
    }
    if (strcmp(args[i],"all") == 0 ){
      uses->u[uses->nu].hitid = -1; /* all hits ! */
      uses->u[uses->nu].useid = uses->id;
      uses->nu++;
    } else if ( (c = strchr(args[i],'-')) != NULL){
      int ilow,iup,nnew,j;
      *c='\0';
      ilow=atoi(args[i]);
      iup=atoi(c+1);
      nnew =  1 + iup - ilow;
      for (j = 0 ; j <  nnew ; j++){
	if (uses->nu >= n_allocated-1){
	  n_allocated += RDMC_MAXTOKEN_PER_LINE;
	  uses->u = realloc( uses->u , n_allocated*sizeof(mevt_uses_t));
	}
	uses->u[uses->nu].hitid = ilow+j; /* all hits ! */
	uses->u[uses->nu].useid = uses->id;
	uses->nu++;
      }
    }else{
      uses->u[uses->nu].hitid = atoi(args[i]); /* all hits ! */
      uses->u[uses->nu].useid = uses->id;
      uses->nu++;
    }
  } /* for */

  return RDMC_IO_OK;
}

/************************************************/
/***********write functions ********************/
/*************************************************/


int rdmc_wfoot_f2k_1x1_2(const mcfile *fp,const array *a, const mevt *e){
  /* now write the dummy tend line */
  fputs("TEND ? ? ?\n",fp->fp); 
  fputs("END\n",fp->fp); 
  return 0;
}

/****************************************************************************
 * wrhist_amanda() writes history lines to an ASCII file
 ****************************************************************************/

int rdmc_wrhist_f2k_1x1(const mcfile *fp, const char *s, const char *pre)
{
  const char *bol;                               /* beginning of the line */

  int pre_len=0;
  int max_step=0;

  if ((s== NULL)|| (pre == NULL))
    return RDMC_IO_OK;

  /* init */
  pre_len = strlen(pre);
  bol=s;

  /* maximum to write */
  max_step = F2000_MAXLINE - pre_len - RDMC_MAXTOKENLENGTH ; 
  if (max_step <= 0) 
    return RDMC_ILF;
  
  while (1) {
    int do_exit=0;
    int icount=0;
    char line[RDMC_MAXLINE+1]="";          /* temporary copy of the string */

    while (do_exit == 0){ /* print this line */
      if ((*bol)=='\0'){ /* end of string */
	return RDMC_IO_OK;
      }else if ((*bol) == '\n'){ /* end of line */
	if(icount > 0){
	  line[icount]='\0';
	  fputs(pre,fp->fp);fputs(line,fp->fp); putc('\n',fp->fp);
	  /*fprintf(fp->fp, "%s%s\n", pre, line);*/
	}
	++(bol);
	do_exit =1;
      } else if ((*bol) == '\\'){ /* cont line */
	line[icount] = '\0';
	strcat(line,"\\");
	fputs(pre,fp->fp);fputs(line,fp->fp); putc('\n',fp->fp);
	/*fprintf(fp->fp, "%s%s\n", pre, line); */
	++(bol);
	do_exit =1;

      } else if ((icount) >=  (max_step+1)){ /* cont line */
	line[icount]='\0';
	strcat(line,"\\");
	fputs(pre,fp->fp);fputs(line,fp->fp); putc('\n',fp->fp);
	/*fprintf(fp->fp, "%s%s\n", pre, line);*/
	do_exit =1;
      } else { /* just copy */
	line[icount] = *bol;
	++icount;
	++(bol);
	do_exit = 0 ;
      }  
    } /* print this line */
  }

} /* wrhist_amanda() */

/****************************************************************************
 * wrcomment_amanda() writes a comment line to an ASCII file
 ****************************************************************************/

int rdmc_wrcomment_f2k_1x1(const mcfile *fp, const char *s)
{      
  static char *ts=NULL;
  char line[F2000_MAXLINE+1];              /* temporary copy of the string */
  char *eol;                                    /* end of current line */
  char *bol;                                  /* beginning of the line */
  if (s != NULL){
    /* get a local copy due to const decl */
#ifndef CRAY
    ts = alloca(sizeof(char)*(strlen(s)+1));
#else
    ts = realloc(ts,sizeof(char)*(strlen(s)+1));
#endif
    strcpy(ts,s);

    bol=ts;
    while(*bol){
      int thislen;
      eol = strchr(bol,'\n');
      if (eol){
	thislen = 1 + (int) (eol-bol) ;
	if (thislen <=  1){
	  putc('\n',fp->fp);
	/*fprintf(fp->fp, "\n"); */
	}else if (thislen <= F2000_MAXLINE) {
	  *eol='\0';
	  fputs(bol,fp->fp); putc('\n',fp->fp);
	  /*fprintf(fp->fp, "%s\n", bol);*/
	}else {
	  strncpy(line,bol,F2000_MAXLINE-2);
	  line[F2000_MAXLINE-2]='\\';
	  line[F2000_MAXLINE-1]='\0';
	  fputs(line,fp->fp); putc('\n',fp->fp);
	  /* fprintf(fp->fp, "%s\n",line); */
	  eol = bol + (F2000_MAXLINE-4);
	  *(eol + 1) = '!';
	}
	bol=eol+1;
      }else
	break;
    }
  }

  return RDMC_IO_OK;
} /* wrcomment_amanda() */


int rdmc_whead_f2k_1x1_2(const mcfile *fp,const array *geo, const mevt *e)
{
  int iom;                                              /* om index in loop */
  int itrigger;                                    /* trigger index in loop */
  int istat;
  int ifit,iuser;
  int array_defined=0;  /* is the array line to be written) */

  fprintf(fp->fp, "V 2000.%i.%i\n",                     /* write the V flag */
	  AMANDA_ASCII_VERSION, AMANDA_ASCII_MINOR);       /* version number */

  if (fp->creator !=NULL )
    rdmc_wrhist_f2k_1x1(fp,fp->creator,"HI ");           /* write history lines */
  if (fp->comment != NULL)                 /* if there is a common comment */
    rdmc_wrcomment_f2k_1x1(fp, fp->comment);

  if (geo->comment != NULL)                 /* if there is a common comment */
    rdmc_wrcomment_f2k_1x1(fp, geo->comment);

  /* thest if array is to be written */
  if ((geo->is_calib.geo)
      || ( geo->is_calib.tdc)
      || (geo->is_calib.adc)
#if USES_TOT_CAL /* not in f2000 1.1 */
      || (geo->is_calib.tot)
#endif
      || (geo->nch >0 )
      ){
    array_defined =1;
  }

  /* first write the geometry Block header */
  if (array_defined) {     /* if there was a detector geometry */
    fprintf(fp->fp,"ARRAY %s %g %g %g %i %i\n",
	    rdmc_amanda_sdet(geo->id), 
	    geo->longitude, geo->lattitude,
	    geo->depth, geo->nstr, geo->nch);         /* 'ARRAY' flag */
  }  
  


  /**************** now the trigger block ***********************/
  
  for (itrigger = 0; itrigger < geo->n_trigger; itrigger++) { 
    int i;
    fprintf(fp->fp, "TRIG_DEF %i %s", /* write trigger lines */
	    geo->def_trig[itrigger].id, 
	    geo->def_trig[itrigger].tag);
    for (i=0 ; i< geo->def_trig[itrigger].nwords ; i++){
      putc(' ',fp->fp); fputs(geo->def_trig[itrigger].words[i],fp->fp);
      /*fprintf(fp->fp," %s",geo->def_trig[itrigger].words[i]); */
    }
    putc('\n',fp->fp); /*fprintf(fp->fp,"\n"); */
    if (geo->def_trig[itrigger].npars >0 ){ /* write trigger par lines */

      fputs("TRIG_PAR ",fp->fp);fputs(geo->def_trig[itrigger].tag,fp->fp);
      /* fprintf(fp->fp, "TRIG_PAR %s"
	 ,geo->def_trig[itrigger].tag); */
      for (i=0 ; i< geo->def_trig[itrigger].npars ; i++){
	putc(' ',fp->fp); fputs(geo->def_trig[itrigger].pars[i],fp->fp);
	/*fprintf(fp->fp," %s",geo->def_trig[itrigger].pars[i]);*/
      }
      putc('\n',fp->fp);/* fprintf(fp->fp,"\n");  */
    }
  } /* for itrigger */

  
  /** write now the calibration ('K???') of all channels, if there is one *****/
  /* first theKH line */
  if ( geo->is_calib.tdc
       || geo->is_calib.adc
#if USES_TOT_CAL /* not in f2000 1.1 */
       || geo->is_calib.tot
#endif
       || geo->is_calib.utc){
    fputs("KH",fp->fp); /*fprintf(fp->fp,"KH"); */
    if ( geo->is_calib.adc)
      fputs(" ADC",fp->fp); /*fprintf(fp->fp," ADC"); */
    if ( geo->is_calib.tdc)
      fputs(" TDC",fp->fp); /* fprintf(fp->fp," TDC"); */
    if ( geo->is_calib.utc)
      fputs(" UTC",fp->fp);  /* fprintf(fp->fp," UTC");*/
#if USES_TOT_CAL 
    if ( geo->is_calib.tot)
      fputs(" TOT",fp->fp);  /* fprintf(fp->fp," TOT"); */
#endif
    putc('\n',fp->fp); /* fprintf(fp->fp,"\n"); */
  }
  
  /* now the stat_def -> these are only strings ! */
  for (istat=0 ; istat < geo->n_stat ; istat++){
    int i;
    fputs("STAT_DEF ",fp->fp); fputs(geo->def_stat[istat].tag,fp->fp);
    /* fprintf(fp->fp, "STAT_DEF"); fprintf(fp->fp, " %s",geo->def_stat[istat].tag); */
    for (i=0 ; i< geo->def_stat[istat].nwords ; i++){
      putc(' ',fp->fp); fputs(geo->def_stat[istat].words[i],fp->fp);
      /* fprintf(fp->fp," %s",geo->def_stat[istat].words[i]); */
    }
    putc('\n',fp->fp);/* fprintf(fp->fp,"\n"); */
    /* write STAT_DEF lines */
  }
  
  /* now the FIT_DEF -> these are  ! */
  for (ifit=0 ; ifit < geo->n_fit ; ifit++){
    int i;
    fprintf(fp->fp, "FIT_DEF %i %s", /* write trigger lines */
	    geo->def_fit[ifit].id, 
	    geo->def_fit[ifit].tag);
    for (i=0 ; i< geo->def_fit[ifit].nwords ; i++){
      putc(' ',fp->fp); fputs(geo->def_fit[ifit].words[i],fp->fp);
      /* fprintf(fp->fp," %s",geo->def_fit[ifit].words[i]);*/
    }
    putc('\n',fp->fp); /*fprintf(fp->fp,"\n"); */
    /* write FIT_DEF lines */
  }
  
  /********* write now the positions ('OM') of all channels ***************/
  
  if (geo->is_calib.geo){
    for (iom = 0 ; iom < geo->nch ; iom++) {
      fprintf(fp->fp,"OM %i %i %i %g %g %g %2s %s %s %g %g\n",
	      iom+1, 
	      geo->clust[iom]+1,
	      geo->str[iom],
	      geo->x[iom], geo->y[iom], geo->z[iom],
	      rdmc_amanda_sori(geo->costh[iom]),
	      rdmc_amanda_spmtid(geo->type[iom]),
	      rdmc_amanda_itostr(geo->serial[iom],RDMC_NA),
	      geo->sensit[iom],
	      geo->thresh[iom]);
    } /* for iom */  
  }
  
  /*
   * Write ADC/TDC/UTC calibration entries
   */
  if (geo->is_calib.adc)
    for (iom = 0; iom < geo->nch; iom++)
      fprintf(fp->fp,"KADC %i %g %g %g\n",
	      iom+1,
	      geo->cal[iom].ped,
	      geo->cal[iom].beta_a,
	      geo->cal[iom].kappa);
  
  if (geo->is_calib.tdc)
    for (iom = 0; iom < geo->nch; iom++)
      fprintf(fp->fp,"KTDC %i %g %g %g\n",
	      iom+1,
	      geo->cal[iom].beta_t,
	      geo->cal[iom].t_0,
	      geo->cal[iom].alpha_t);
  if (geo->is_calib.tot)
#if USE_TOT_CAL
    for (iom = 0; iom < geo->nch; iom++)
      fprintf(fp->fp,"KTOT %i %g %g %g\n",
	      iom+1,
	      geo->cal[iom].ped_tot,
	      geo->cal[iom].beta_tot,
	      geo->cal[iom].kappa_tot);
#endif
  if (geo->is_calib.utc){
    fprintf(fp->fp, "KUTC %s %i.%09i\n",geo->cal_utc.utc_src
	    ,geo->cal_utc.secs,geo->cal_utc.nsecs);
  }
  

  /* now the USER_DEF -> these are  ! */
  for (iuser=0 ; iuser < geo->n_user ; iuser++){
    int i;
     fputs("USER_DEF ",fp->fp); fputs(geo->def_user[iuser].tag,fp->fp);
     /* fprintf(fp->fp, "USER_DEF %s", geo->def_user[iuser].tag); */
     for (i=0 ; i< geo->def_user[iuser].nwords ; i++){
       putc(' ',fp->fp); fputs(geo->def_user[iuser].words[i],fp->fp);
       /*fprintf(fp->fp," %s",geo->def_user[iuser].words[i]); */
     }
     putc('\n',fp->fp);/* fprintf(fp->fp,"\n");*/
    /* write USER_DEF lines */
  }

  /* now write the tbegin line */
  if (geo->tbegin > 0){ 
    struct tm *gmt;

    gmt = gmtime(&(geo->tbegin));
       fprintf(fp->fp, "TBEGIN %i %i %i\n", /* write trigger lines */
	       1900+gmt->tm_year,gmt->tm_yday+1,
	       gmt->tm_hour*3600+gmt->tm_min*60+gmt->tm_sec
	       );
  } else {
    fputs( "TBEGIN ? ? ?\n",fp->fp); /*fprintf(fp->fp, "TBEGIN ? ? ?\n"); */
  }

  return 0;
} /* warr_amanda() */


/****************************************************************************
 * function wevt_amanda() writes an event to a amanda-like file
 ****************************************************************************/

int rdmc_wevt_f2k_1x1_2(const mcfile *fp,const array *ar, const mevt *event)
{
  int itok, i;
  char *upt;
  float channel;
  
  unsigned long format_version= 100*fp->fmajor + fp->fminor;
   
  /*
   * write Event header and comments
   */

  /* i have an event 0 at every begin ??? throw it away */
  /*  if(event->enr == 0)return 0 ; */

  {
    int gpsyear,gpsday;
    rdmc_mjd_to_gps(event->mjd, &gpsyear, &gpsday);/* GPS/UT year and day */
    fprintf(fp->fp, "EM %i %i %s", event->enr,event->nrun
	    ,rdmc_amanda_itostr(gpsyear,RDMC_NA));

    fprintf(fp->fp, " %s %i.%09i %g\n"
	    ,rdmc_amanda_itostr(gpsday,RDMC_NA)
	    ,event->secs,(event->nsecs>=0)?event->nsecs:0,
	    event->t_offset);
  }

  if (event->comment != NULL)                /* if there is an event comment */
    rdmc_wrcomment_f2k_1x1(fp, event->comment);


  /*
   * write tracks
   */
  for (itok=0; itok<event->ntrack; itok++){ /* write  the 'TR' lines */
    fprintf(fp->fp,"TR %i %s %s %g %g %g %g %g %s"
	    , event->gen[itok].tag
	    ,rdmc_amanda_itostr(event->gen[itok].parent,RDMC_PARENT_NA)
	    ,rdmc_amanda_spartid(event->gen[itok].id)
	    ,event->gen[itok].x,
	    event->gen[itok].y,
	    event->gen[itok].z,
	    acos(event->gen[itok].costh)/PID180,
	    event->gen[itok].phi/PID180
	    ,rdmc_amanda_ftostr(event->gen[itok].length,RDMC_BIG));
    fprintf(fp->fp, " %s",rdmc_amanda_ftostr(event->gen[itok].e/1000.,RDMC_NA));
    fprintf(fp->fp, " %s\n",rdmc_amanda_ftostr(event->gen[itok].t,RDMC_TDC_NA));
  }

  /*
   * write hits
   */
  for (itok = 0; itok < event->nhits; itok++){
    fprintf(fp->fp, "HT %i",event->h[itok].ch+1);

    if (event->h[itok].amp == RDMC_REPEAT_CH)
      fputs(" *",fp->fp);
    else{
      putc(' ',fp->fp); fputs(rdmc_amanda_ftostr(event->h[itok].amp,RDMC_NA),fp->fp); }
    /* fprintf(fp->fp, " %s",rdmc_amanda_ftostr(event->h[itok].amp,RDMC_NA));*/

    fprintf(fp->fp, " %i",event->h[itok].id);
    if (event->h[itok].mt == RDMC_PARENT_NOISE)
      fputs(" N",fp->fp);
    else if(event->h[itok].mt == RDMC_PARENT_AFPULS)  
      fputs(" A",fp->fp);
    else{
      putc(' ',fp->fp); 
      fputs(rdmc_amanda_itostr(event->h[itok].mt,RDMC_PARENT_NA),fp->fp);
    }
    putc(' ',fp->fp); fputs(rdmc_amanda_ftostr(event->h[itok].t,RDMC_TDC_NA),fp->fp);
    putc(' ',fp->fp); fputs(rdmc_amanda_ftostr(event->h[itok].tot,RDMC_NA),fp->fp);
    if ( (format_version ) >= 102 ){ /* only in case of minor f2k.1.2 and further version */
	if (event->h[itok].hstat.n_tdc_edges || event->h[itok].hstat.tdc_flag){
	    putc(' ',fp->fp);
	    if (event->h[itok].hstat.tdc_flag) putc('>',fp->fp);
	    fprintf(fp->fp, "%i",event->h[itok].hstat.n_tdc_edges);
	}
    }

    putc('\n',fp->fp);
    /*      fprintf(fp->fp, " %s",rdmc_amanda_itostr(event->h[itok].mt,RDMC_PARENT_NA));
    fprintf(fp->fp, " %s",
	    rdmc_amanda_ftostr(event->h[itok].t,RDMC_TDC_NA));
    fprintf(fp->fp, " %s\n", 
	    rdmc_amanda_ftostr(event->h[itok].tot,RDMC_NA));
    */
  }


  /*
   * write waveforms
   */
  /*  float channel; */
  
  if ( (format_version ) >= 102 ){ /* only in case of minor f2k.1.2 and further version */
      for (itok = 0; itok < event->nwf; itok++){
	  
	  channel = event->wf[itok].om + event->wf[itok].ch/100.;
	  
	  fprintf(fp->fp, "WF %4.2f", channel); 
	  putc(' ',fp->fp);  
	  fputs(rdmc_amanda_itostr(event->wf[itok].id, RDMC_NA),fp->fp); 
	  putc(' ',fp->fp);  
	  fputs(rdmc_amanda_itostr(event->wf[itok].pa, RDMC_PARENT_NA),fp->fp);
	  fprintf(fp->fp, " %i",event->wf[itok].ndigi);
	  putc(' ',fp->fp);  
	  fputs(rdmc_amanda_ftostr(event->wf[itok].t_start,RDMC_NA),fp->fp); 
	  putc(' ',fp->fp);  
	  fputs(rdmc_amanda_ftostr(event->wf[itok].t_bin,RDMC_NA),fp->fp); 
	  if ( (format_version ) >= 200401 ){ /* only in case of minor f2k.2004.1 and further version */
	      putc(' ',fp->fp);  
	      fputs(rdmc_amanda_ftostr(event->wf[itok].baseline,RDMC_WF_BASELINE_NA),fp->fp); 
	  }
	  for ( i=0; i<event->wf[itok].ndigi; i++){
	      /*      if(i!=0&&i%10 == 0) fprintf(fp->fp, "\n&"); // continuation line  */
	      /* doesn't work with reading */
	      /*     fprintf(fp->fp, " %6.1f",event->wf[itok].digi[i]);*/
	      putc(' ',fp->fp); 
	      fputs(rdmc_amanda_ftostr((int)((event->wf[itok].digi[i])*10)/10.,RDMC_TDC_NA),fp->fp);
	  }
	  putc('\n',fp->fp);
      }
  }
	
  /*
   * write trigger lines
   */
  for (itok = 0; itok < event->ntrig ; itok++){
    int id;
    int i;
    id = event->ptrig[itok].id;
    if ( (id < ar->n_trigger ) ){
      fputs("TRIG ",fp->fp); fputs(ar->def_trig[id].tag,fp->fp);
      /*fprintf(fp->fp, "TRIG %s", ar->def_trig[id].tag); */
      for ( i=0  ;  i < event->ptrig[itok].nval ; i++){
	fprintf(fp->fp, " %g" ,event->ptrig[itok].val[i]);
      }
      putc('\n',fp->fp); /*fprintf(fp->fp, "\n"); */
      /* 
	 write trigger uses lines
      */
      upt = rdmc_amanda_uses_to_str(event->ntrig_uses,event->trig_uses,itok);
      if (upt[0] != '\0'){ 
	fputs(upt,fp->fp); /*fprintf(fp->fp, "%s",upt); */
      }
    }
  }


  /* 
     write status  lines 
  */
  for (itok = 0; itok < event->nstat ; itok++){
    int id;
    int i;
    id = event->status[itok].id;
    if (id < ar->n_stat  ){
       fputs("STATUS ",fp->fp); fputs(ar->def_stat[id].tag,fp->fp);
       /* fprintf(fp->fp, "STATUS %s", ar->def_stat[id].tag); */
      for ( i=0  ;  i < event->status[itok].nval ; i++){
	fprintf(fp->fp, " %g" ,event->status[itok].val[i]);
      }
      putc('\n',fp->fp); /* fprintf(fp->fp, "\n"); */
    }
  }

  
  /*
   * write fits
   */
  for (itok = 0; itok < event->nfit; itok++){
    fprintf(fp->fp, "FIT %i %s %g %g %g %g %g %g",
	    event->rec[itok].tag,
	    rdmc_amanda_spartid(event->rec[itok].id),
	    event->rec[itok].x,
	    event->rec[itok].y,
	    event->rec[itok].z,
	    acos(event->rec[itok].costh)/PID180,
	    event->rec[itok].phi/PID180,
	    event->rec[itok].t);
    putc(' ',fp->fp);  fputs(rdmc_amanda_ftostr(event->rec[itok].length,RDMC_BIG),fp->fp);
    putc(' ',fp->fp);  fputs(rdmc_amanda_ftostr(event->rec[itok].e/1000.,0.),fp->fp);
    putc('\n',fp->fp); 
    /*
    fprintf(fp->fp, " %s"
	    ,rdmc_amanda_ftostr(event->rec[itok].length,RDMC_BIG));
    fprintf(fp->fp, " %s\n",rdmc_amanda_ftostr(event->rec[itok].e/1000.,0.));
    */
    
    {/*      write fresult   */
      int id;
      int i;
      id = event->fresult[itok].id;
      if ((id>=0 ) && (id < ar->n_fit )){
	fprintf(fp->fp, "FRESULT %i", 
		ar->def_fit[id].id);
	for ( i=0  ;  i < event->fresult[itok].nval ; i++){
	  putc(' ',fp->fp);  
	  fputs(rdmc_amanda_ftostr(event->fresult[itok].val[i], 
				   RDMC_SPECIAL_NA),fp->fp);
	  /* fprintf(fp->fp, " %s",
	     rdmc_amanda_ftostr(event->fresult[itok].val[i],
	     RDMC_SPECIAL_NA)); */
	}
	putc('\n',fp->fp); /* fprintf(fp->fp, "\n"); */
      }
    }

    /* 
       fuses lines
    */
    upt = rdmc_amanda_uses_to_str(event->nfit_uses,event->fit_uses,itok);
    if (upt[0] != '\0'){ 
      fputs(upt,fp->fp); /* fprintf(fp->fp, "%s",upt); */
    }
  }


  /*
   * write user defined lines
   */
  /*  fprintf(stderr,"*** rdmc %i\n",event->nuser); */
  for (itok = 0; itok < event->nuser ; itok++){
    int id;
    int i;
    id = event->user[itok].id;
    /*    fprintf(stderr,"*** rdmc %i %i %i\n",id,ar->n_user,ar->def_user[0].nwords); */
    if (id < ar->n_user  ){
      fputs("US ",fp->fp); fputs(ar->def_user[id].tag,fp->fp);
      /* fprintf(fp->fp, "US %s", ar->def_user[id].tag); */
      for ( i=0  ;  i < event->user[itok].nval ; i++){
	fprintf(fp->fp, " %g" ,event->user[itok].val[i]);
      }
      putc('\n',fp->fp);  /*fprintf(fp->fp, "\n"); */
    }
  }
  fputs("EE\n",fp->fp); /* fprintf(fp->fp, "EE\n"); */

  return 0;
} /* function wevt_amanda() */





/****************************************************************************
 ********************************** E O F ***********************************
 ****************************************************************************/



