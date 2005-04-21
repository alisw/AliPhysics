
/* implement functions that help with the f2k format */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "rdmc.h"
#include "amanda.h"
#include "f2k.h"

/* copies the string src to dest and removes inline coments, leading blanks */
/* returns the number of characters copied */
static int rdmc_f2k_strcleancpy(char *dest, char *src);



/****************************************************************************/
int rdmc_f2k_dummy_parser(mcfile *fp, array *a, mevt *e, void *misc){ 
  return RDMC_IO_OK;
}
/****************************************************************************/
int rdmc_f2k_dummy_event_writer(const mcfile *fp, const array *a, const mevt *e){ 
  return RDMC_IO_OK;
}
/****************************************************************************/

/****************************************************************************/
/* returns 1 if the format version is supportet */
int rdmc_is_f2000_supported(int major, int minor){
  static struct {
    int major ;
    int minor ;
  } supportet[] = 
    { {1,1} ,{1,2} , {2004,1},  {0,0} };
  int i=0;

  while (( supportet[i].major != 0) && ( supportet[i].minor != 0)){
    if (supportet[i].major == major){
      if (supportet[i].minor == minor){
	return 1;
      }
    }
    i++;
  }
  return 0;
}

/* stores specific error information int f2k file pointer for print */
void rdmc_f2k_errorlog(mcfile *fp){
  rdmc_f2k_buffer_t *f2k_buff = fp->info.f2000.f2k_buffer;
  int iline= f2k_buff->iline;
  char *token_pt = f2k_buff->line_pt[iline];
  char *target_pt = fp->info.f2000.parse_line;
  char *end_pt=NULL;

  if (*target_pt)
    return; /* already one error */
  *target_pt='\0';
  fp->info.f2000.errline  = fp->fpos - f2k_buff->lines + iline +1;
  /* try to recover the tokenized line  arghhhh */
  if (iline+1 >= f2k_buff->lines) 
    end_pt = &(f2k_buff->buff[f2k_buff->used-1]);
  else
    end_pt = f2k_buff->line_pt[iline+1];
  while (token_pt < end_pt ){
    if(*token_pt == '\0')
      *target_pt=' ';
    else
      *target_pt = *token_pt;
    target_pt++;      
    token_pt++;
  }
  *target_pt = '\0';
  return;
}

/****************************************************************************/
/****************************************************************************/
/* text buffer functions *********************/
void rdmc_push_f2k_buffer(rdmc_f2k_buffer_t *b, char *s, 
			  const f2000_line_t * type_def){
  int nnew;

  /* check the size and increase it if necessary */
  while ( (b->used + RDMC_MAXLINE + 1) >= b->ntot){
    register int i;
    long int pedestal;
    char* bsave = b->buff; /* is needed to save the former location */

#if 0
    fprintf(stderr,"bused=%i lused=%i btot=%i ltot=%i %.10s\n",b->used, b->lines, b->ntot,b->lines_tot,s); 
#endif

    b->ntot += F2K_BUFFERSIZE;
    b->buff = realloc(b->buff,b->ntot*sizeof(char));

    /* now correct the line_pointer pointers */
    pedestal = (long int ) (b->buff - bsave); /* amount buffer was mooved */
#if 0
    fprintf(stderr," %p %p %p %i %p\n",bsave,b->buff,b->line_pt[0],pedestal,b->line_pt[0]+pedestal);
#endif
    for(i=0 ; i< b-> lines; i++){
       b->line_pt[i] += pedestal; /* only line pointers have to be cotrected */

    }
  }

  while( (b->lines+1) >= b->lines_tot){
    b->lines_tot += F2K_LINE_BUFFERSIZE;
    b->line_pt = realloc(b->line_pt,b->lines_tot*sizeof(char *));
    b->type_def = realloc(b->type_def
			  ,b->lines_tot*sizeof(const f2000_line_t *));
  }

  /* copy the string */
    nnew = rdmc_f2k_strcleancpy(&(b->buff[b->used]),s);
    if(nnew){
      b->type_def[b->lines] = type_def;
      b->line_pt[b->lines] = &(b->buff[b->used]);
      b->used += nnew;
      b->lines += 1;
      b->line_pt[b->lines] = NULL; /* flag the last line to NULL */
    }
}

void rdmc_init_f2k_buffer(rdmc_f2k_buffer_t *b){
  b->buff = NULL; b->line_pt = NULL ; b->type_def = NULL;
  b->ntot = b->lines_tot = 0;
  rdmc_reset_f2k_buffer(b);
}

void rdmc_reset_f2k_buffer(rdmc_f2k_buffer_t *b){
  b->used = b->lines = b->iline = 0;

  b->ntot=F2K_BUFFERSIZE;
  b->buff = realloc(b->buff,b->ntot*sizeof(char));

  b->lines_tot = F2K_LINE_BUFFERSIZE;
  b->line_pt = realloc(b->line_pt,b->lines_tot*sizeof(char *));
  b->type_def = realloc(b->type_def,b->lines_tot*sizeof(f2000_line_t *));

}

void rdmc_clear_f2k_buffer(rdmc_f2k_buffer_t *b){
  b->used = b->lines = b->iline = 0;
}

void rdmc_unlink_f2k_buffer(rdmc_f2k_buffer_t *b){
  free(b->buff);
  free(b->type_def);
  free(b->line_pt);
  b->ntot = b->lines_tot = 0;
  b->buff = NULL; b->line_pt = NULL ; b->type_def = NULL;
}

/*********************************************/
static int rdmc_f2k_strcleancpy(char *dest, char *src){
  int ncp=0;                           /* number of chars copied to dest */
  char *spt=src;
  char *dpt=dest;
  int leading_space=1;              /* flag to mark the begining of parsing */

  /* copy the string */
  while (*spt) {
    if (leading_space){ /* scip leading spaces */ 
      if (isspace(*spt)){
	++(spt);
	continue;
      }else{
	leading_space = 0;
      }
    }
    /* now we are pointing to the first non-white */
    /* kill inline comments */
    *dpt = *spt;
    ++(dpt);
    ++(spt);
    ++ncp;
  }
  
  /* append a \0 terminator */
  *dpt='\0';
  ++ncp;
  return ncp;
}


/***************************************************************************/

/**********************************************************************
 * scans a nonempty line into tokens
 **********************************************************************/
char ** rdmc_f2k_tokenize(char *s, int *nargs) {
  static char ** f2k_tokens = NULL; /* buffer for tokens */
  static int f2k_maxtoken =0;       /* buffer for tokens */
  const char delim[]=" \t\n";
  char *tok_pt;
  
  *nargs=0;
  tok_pt = strtok(s,delim); /* find first non-whitespace char */
  while (tok_pt){
    if (*nargs >= f2k_maxtoken ){
      f2k_maxtoken += RDMC_MAXTOKEN_PER_LINE + 1; 
      f2k_tokens = realloc(f2k_tokens,f2k_maxtoken*sizeof(char *));
    }
    f2k_tokens[*nargs] = tok_pt;
    *nargs += 1;
    tok_pt = strtok(NULL,delim); /* find first non-whitespace char */
  }
  f2k_tokens[*nargs + 1] = NULL;

  return f2k_tokens;
} /* amanda_tokenize() */



/****************************************************************************/
/****************************************************************************/
/** string conversion functions                    */

/****************************************************************************
 * return the OM id, with the experiment prefix (see siegmund format 
 *  description about these ids)
 ****************************************************************************/
int rdmc_amanda_iomid(const char *str,int array_id){
  return rdmc_amanda_ipmtid(str) 
    + 100000*(array_id / 100);
}

const char * rdmc_amanda_somid(int id){
  return rdmc_amanda_spmtid(id);
}

/****************************************************************************
 * return the OM id, with the experiment prefix (see siegmund format 
 *  description about these ids)
 ****************************************************************************/
const char *rdmc_amanda_spmtid(int type)
{
  int pmtid, sphereid, dataid;
  int i;
  static char typestring[RDMC_MAXTOKENLENGTH];
  char *sphere, *pmt, *data;

  pmtid = type%100;
  dataid = (type/100)%100;
  sphereid = (type/10000)%10;

  pmt = "std";
  for (i = 0; rdmc_pmt_idtable[i].name != NULL; i++)
    if (rdmc_pmt_idtable[i].id == pmtid) {
      pmt =  rdmc_pmt_idtable[i].name;
      break;
    }

  sphere = "std";
  for (i = 0; rdmc_sphere_idtable[i].name != NULL; i++)
    if (rdmc_sphere_idtable[i].id == sphereid) {
      sphere =  rdmc_sphere_idtable[i].name;
      break;
    }

  data = "std";
  for (i = 0; rdmc_datatrans_idtable[i].name != NULL; i++)
    if (rdmc_datatrans_idtable[i].id == dataid) {
      data =  rdmc_datatrans_idtable[i].name;
      break;
    }

  sprintf(typestring, "%s-%s-%s", pmt, sphere, data);

  return typestring;

} /* rdmc_amanda_spmtid */


/****************************************************************************
 * return the OM id, without the experiment prefix (see siegmund format 
 *  description about these ids)
 ****************************************************************************/
int rdmc_amanda_ipmtid(const char *str)
{
  char *tok;
  char *c;
  int pmtid, dataid, sphereid;
  int i;
  char s1[RDMC_MAXTOKENLENGTH];

  strcpy(s1,str);

  tok = strtok(s1, "-");
  if (tok == NULL) return 0;
  for (c=tok ; *c ; c++)
    *c=tolower((int) *c);

  pmtid = 0;
  for (i = 0; rdmc_pmt_idtable[i].name != NULL; i++)
    if (strcmp(rdmc_pmt_idtable[i].name, tok) == 0) {
      pmtid =  rdmc_pmt_idtable[i].id;
      break;
    }
  
  tok = strtok(NULL, "-");
  if (tok == NULL) return pmtid;
  for (c=tok ; *c ; c++)
    *c=tolower((int) *c);
  
  sphereid = 0;
  for (i = 0; rdmc_sphere_idtable[i].name != NULL; i++)
    if (strcmp(rdmc_sphere_idtable[i].name, tok) == 0) {
      sphereid =  rdmc_sphere_idtable[i].id;
      break;
    }

  tok = strtok(NULL, "-");
  if (tok == NULL) return 10000*sphereid + pmtid;
  for (c=tok ; *c ; c++)
    *c=tolower((int) *c);

  dataid = 0;
  for (i = 0; rdmc_datatrans_idtable[i].name != NULL; i++)
    if (strcmp(rdmc_datatrans_idtable[i].name, tok) == 0) {
      dataid =  rdmc_datatrans_idtable[i].id;
      break;
    }
  
  return 10000*sphereid + 100*dataid + pmtid;

} /* rdmc_amanda_ipmtid() */

/****************************************************************************
 * build a string containing the Detector name, or "Neutrino_telescope" if
 * unknown
 ****************************************************************************/

const char *rdmc_amanda_sdet(int geoid)
{
  int i;
  for (i = 0; rdmc_detector_idtable[i].name != NULL; i++)
    if (rdmc_detector_idtable[i].id == geoid) 
      return  rdmc_detector_idtable[i].name;
  return rdmc_detector_idtable[0].name;
} /* rdmc_amanda_sdet() */

int rdmc_amanda_idet(const char *detnam)
{
  int i;
  char tmp_detnam[RDMC_MAXTOKENLENGTH];
  char *c;
  
  strncpy(tmp_detnam,detnam,RDMC_MAXTOKENLENGTH-1);
  tmp_detnam[RDMC_MAXTOKENLENGTH-1]='\0';

  c=tmp_detnam;
  while(*c){
    *c = tolower((int) *c);
    c++;
  }

  for (i = 0; rdmc_detector_idtable[i].name != NULL; i++)
    if (strcmp(rdmc_detector_idtable[i].name, tmp_detnam) == 0)
      return  rdmc_detector_idtable[i].id;

  /* maybe, we find a substring? */
  for (i = 0; rdmc_detector_idtable[i].name != NULL; i++)
    if (strstr(rdmc_detector_idtable[i].name, tmp_detnam) != NULL)
      return  rdmc_detector_idtable[i].id;
  return rdmc_detector_idtable[0].id;

} /* rdmc_amanda_detid() */


/****************************************************************************
 * build a string containing the particle name
 ****************************************************************************/
const char *rdmc_amanda_spartid(int id)
{
  int i;
  static char a_primary[]="Axxx";
  static char z_primary[]="Zxxx";

  if ((id>=Z_PRIMARY) && (id < (Z_PRIMARY+200))){ /* primary Zxx */
    int idd=id-Z_PRIMARY;
#if 0
    if(idd <10)
      sprintf(z_primary,"Z00%1i",idd);
    else if(idd <100)
      sprintf(z_primary,"Z0%2i",idd);
    else
      sprintf(z_primary,"Z%3i",idd);
#else
      sprintf(z_primary,"Z%i",idd);
#endif

    return z_primary;
  }else if ((id>=A_PRIMARY) && (id < A_PRIMARY+500)){ /* primary Axxx undocumented !! */ 
    int idd=id-A_PRIMARY;
#if 0
    if(idd <10)
      sprintf(a_primary,"A00%1i",idd);
    else if(idd <100)
      sprintf(a_primary,"A0%2i",idd);
    else
      sprintf(a_primary,"A%3i",idd);
#else
    sprintf(a_primary,"A%i",idd);
#endif
    return a_primary;
  }

  for (i = 0; rdmc_particle_idtable[i].name != NULL; i++)
    if (rdmc_particle_idtable[i].id == id) 
      return  rdmc_particle_idtable[i].name;

  return rdmc_amanda_spartid(0); /* return unknown */

} /* rdmc_amanda_spartid */


int rdmc_amanda_ipartid(const char *s)
{
  int i;
  int idd;
  char name[RDMC_MAXTOKENLENGTH];

  {
    char *c=name;
#if 1
    const char *sp=s;
    while ( *sp ){
      *c=tolower((int) *sp);  
      sp++; 
      c++;
    }
    *(c++)='\0';
#else
    strcpy(name,s);
    for ( ; *c ; c++)
      *c=tolower(*c);
#endif
  }

  for (i = 0; rdmc_particle_idtable[i].name != NULL; i++)
    if (strcmp(rdmc_particle_idtable[i].name, name) == 0)
      return  rdmc_particle_idtable[i].id;

  if (name[0] == 'z'){ /* CR primary charge Z */
    idd=0;
    if( 1 != sscanf(name,"z%d",&idd))
      return 0;
    if (isnan(idd))
      return 0;
    idd += Z_PRIMARY;
    if ( (idd <  Z_PRIMARY) && (idd >=  Z_PRIMARY+200) )
      return 0;
    else
      return idd;
  }else if (name[0] == 'a'){ /* CR primary mass A */
    idd=0;
    if( 1 != sscanf(name,"a%d",&idd))
      return 0;
    if (isnan(idd))
      return 0;
    idd += A_PRIMARY;
    if ( (idd <  A_PRIMARY) && (idd >=  A_PRIMARY+500) )
      return 0;
    else
      return idd;
  }
  return 0;

} /* rdmc_amanda_ipartid() */

int rdmc_amanda_mhit_stat(mhit_stat_t *hstat, char *stat_s){
  char *cpt=stat_s;
  long int iedges;

  if (*cpt == '\0')
    return RDMC_IO_OK;
  if( *cpt == '>'){
    hstat->tdc_flag = 1;
    cpt++;
  }
  iedges = strtol(cpt, (char **)NULL, 10);
  if ((iedges == LONG_MIN )|| (iedges == LONG_MAX ))
    return RDMC_ILF;
  else
    hstat->n_tdc_edges = iedges;

  return RDMC_IO_OK;
}


/****************************************************************************/
/****************************************************************************/
/** value conversion functions                    */
/****************************************************************************
 * make an integer from a string. If "na" or "?" return RDMC_NA
 ****************************************************************************/

int rdmc_amanda_strtoi(const char *str, int default_na)
{
  int r;
  /* f2000 writes  not na, but ? for NA */   
  if (strstr(str,"?") != NULL)  
    return default_na;
  if (strstr(str,"na") != NULL) /* old rdmc specific value */
    return default_na;
  if (strstr(str,"*") !=NULL)
    return  RDMC_REPEAT_CH;
  if  (strstr(str,"-inf") !=NULL)
    return -RDMC_IINFTY;
  if  (strstr(str,"inf") !=NULL)
    return RDMC_IINFTY;
  if  (strstr(str,"N") !=NULL)
    return RDMC_PARENT_NOISE;
  if  (strstr(str,"A") !=NULL)
    return RDMC_PARENT_AFPULS;


  r = atoi(str);
  if (isnan( (float) r ))
    return  default_na;
  else
    return r;

} /* rdmc_amanda_strtoi() */

/****************************************************************************
 * make a float from a string. If "na" return default_na
 ****************************************************************************/
double  rdmc_amanda_strtof(char *str, double default_na)
{
  double r;
  /* f2000 writes  not na, but ? for NA */   
  if (strstr(str,"?") != NULL)  
    return  default_na;
  if (strstr(str,"na") != NULL) /* old rdmc specific value */
    return default_na;
  if (strstr(str,"*") !=NULL)
    return  RDMC_REPEAT_CH;
  if  (strstr(str,"-inf") !=NULL)
    return -RDMC_FINFTY;
  if  (strstr(str,"inf") !=NULL)
    return RDMC_FINFTY;
#if 0 /* parent is never a float but an int */
  if  (strstr(str,"N") !=NULL)
    return RDMC_PARENT_NOISE;
  if  (strstr(str,"A") !=NULL)
    return RDMC_PARENT_AFPULS;
#endif
//  r=strtof(str,NULL);
  r = atof(str);
  if (isnan( (float) r ))
    return  default_na;
  else
    return r;

} /* rdmc_amanda_strtof() */

/****************************************************************************
 * make a float to a string.
 ****************************************************************************/
char * rdmc_amanda_ftostr(double f, double default_na)
{
  static char str[RDMC_MAXTOKENLENGTH];
  /* f2000 writes  not na, but ? for NA */   
  if (isnan(f)) return  "NaN";
  if ( fabs(f - default_na) < 1.e-5 )  
    return  "?";
  if (f >= RDMC_FINFTY/2)
    return "inf";
  if (f <= -RDMC_FINFTY/2)
    return "-inf";
  sprintf(str,"%g",f);
  return str;
} /* rdmc_amanda_ftostr() */

/****************************************************************************
 * make a int to a string. 
 ****************************************************************************/
char * rdmc_amanda_itostr(int i, int default_na)
{
  static char str[RDMC_MAXTOKENLENGTH];
  /* f2000 writes  not na, but ? for NA */   

  if (isnan(i)) return  "NaN";

  if (i == default_na)  
    return  "?";
  if (i >= RDMC_IINFTY)
    return "inf";
  if (i <= -RDMC_IINFTY)
    return "-inf";
  sprintf(str,"%i",i);
  return str;

} /* rdmc_amanda_itostr() */

/****************************************************************************
 * make a time string sec.nsec. 
 ****************************************************************************/
char * rdmc_amanda_itimetostr(int sec, int nsec){
  static char str[RDMC_MAXTOKENLENGTH];
  static char s1[RDMC_MAXTOKENLENGTH];
  static char s2[RDMC_MAXTOKENLENGTH];

  if (isnan(sec)  || (sec<0)   ) 
    strcpy(s1,"?.");
  else
    sprintf(s1,"%i.",sec);

  if (isnan(nsec) || (nsec<0) || (nsec>= 1e9 )   ) 
    strcpy(s2,"?");
  else
    sprintf(s2,"%09i",nsec);

  strcpy(str,s1);
  strcat(str,s2);
  return str;
} /* rdmc_amanda_timestr() */

/* returns 0 on succes */
int rdmc_amanda_strtimetoi(const char *stime, int *sec, int *nsec){
  char *spt;
  int slen;
  if ( (spt = strchr(stime,'.')) != NULL){
    
    *sec = rdmc_amanda_strtoi(stime,0); /* test only if a valid number */
    *nsec = rdmc_amanda_strtoi(spt+1,0); /* test only if a valid number */
    if (*nsec >0 ){
      slen = strlen(spt+1) - 9;
      if (slen < 0){
	while (slen){
	  *nsec *= 10;
	  slen++;
	}
      }
    }
  }else{
    *sec = rdmc_amanda_strtoi(stime,0); /* test only if a valid number */
    *sec %= 86400 ;  /* seconds of begin of day only is mod(sec,86400)*/
    *nsec = rdmc_amanda_strtoi("?",0); /* test only if a valid number */
  }
  return 0;
}

int rdmc_f2k_y2k(int year){
  year = (year > 0 ) ? year : 1970;
  if (year < 70 )
    year += 2000;
  if (year < 100 )
    year += 1900;
  return year;
}

/****************************************************************************
 * return the orientation (cosinus theta) as a float
 ****************************************************************************/
static const struct {
  char* str;
  float fori;
  float upper;
  float lower;
} rdmc_updown_table[] = {
  {"hr", 0.0 , 0.5 , -0.5 },
  {"HR", 0.0 , 0.5 , -0.5 },
  {"Hr", 0.0 , 0.5 , -0.5 },
  {"up", 1.0 , 1.1 , 0.5 },
  {"UP", 1.0 , 1.1 , 0.5 },
  {"Up", 1.0 , 1.1 , 0.5 },
  {"dn", -1.0, -0.5 , -1.1 },
  {"DN", -1.0, -0.5 , -1.1 },
  {"Dn", -1.0, -0.5 , -1.1 },
  {"down", -1.0, -0.5 , -1.1 },
  {"DOWN", -1.0, -0.5 , -1.1 },
  {"Down", -1.0, -0.5 , -1.1 },
  {NULL , 0.0 , 0.0 , 0.0 }
};

char *rdmc_amanda_sori(float ori){
  int i;
  for (i=0 ; rdmc_updown_table[i].str != NULL ; i++ ){
    if ( (ori  > rdmc_updown_table[i].lower)
	 &&  ( ori  < rdmc_updown_table[i].upper) ){
      return rdmc_updown_table[i].str;
    }
  }
  return rdmc_updown_table[0].str;
}

float rdmc_amanda_fori(char *ori)
{
  int i;
  for (i=0 ; rdmc_updown_table[i].str != NULL ; i++ ){
    if (strcmp(ori,rdmc_updown_table[i].str) == 0) 
      return rdmc_updown_table[i].fori;
  }
  return rdmc_updown_table[0].fori;
} /* rdmc_amanda_ori() */


static const struct {
  int ival;
  double val;
  char *name;
} rdmc_sp_ftable[] = {
  { 1, 1.0, "yes" },
  { 0, 0.0, "no" },
  { 1, 1.0, "true" },
  { 0, 0.0, "false" },
  { RDMC_NA, RDMC_NA, "nan" },
  { RDMC_NA, RDMC_NA, "na" },
  { -1,  -1.0 , NULL }
};

/* catch some special values */
double rdmc_amanda_sptof(char *str){
  double r=RDMC_NA;
  int ir=RDMC_NA,i,slen;

  /* now check for hex numbers */
  if (strstr(str,"0x") == str){
    if ( 1 == sscanf(str,"%i",&ir)){
      r = ir;
      return r;
    }
  }

  /* try to decode directly first */
  if (sscanf(str,"%lg",&r) == 1){
    /*    rdmc_msgprintf("%.18g %s %i",r,str,DBL_DIG); */
    return r;
  }

  /* make a local copy and tronspose to lower */
  {  
    char t[RDMC_MAXTOKENLENGTH];
    slen = strlen(str);
    slen = (slen < (RDMC_MAXTOKENLENGTH-1)) ? slen : RDMC_MAXTOKENLENGTH-1;
    for (i=0 ; i <= slen ; i++){
      t[i] = tolower((int) str[i]);
    } 
    t[slen]='\0';
    
    /* check special values */
    for (i=0; rdmc_sp_ftable[i].name != NULL  ; i++){
      if( strcmp(rdmc_sp_ftable[i].name,t) == 0){
	return rdmc_sp_ftable[i].val;
      }
    }
  } /* check for special values */

  /* now try a last decode: */
  
  return RDMC_SPECIAL_NA;
}


/* create a string the USES lines for itoken */
char * rdmc_amanda_uses_to_str(int n_uses, mevt_uses_t *uses, int id){
  static char *buffer =  NULL;
  static int bsize = 0;
  int bused = 0;
  int i_uses;
  int char_in_line,token_in_line,last_hid;
  char tbuff[RDMC_MAXTOKENLENGTH];
  int cont_flag;
  int t_len;

#if 1
  rdmc_sort_uses(uses, n_uses);
#endif

  /* init the buffer */
  if (!bsize) {
    bsize += RDMC_MAXLINE;
    buffer = malloc(bsize*sizeof(char));
  }
  buffer[0] = '\0';
  bused = char_in_line = token_in_line = 0 ;
  tbuff[0] = '\0';
  last_hid =  RDMC_NA; /* the last seen OM */
  cont_flag=0;

  /* now loop and look */
  for (i_uses = 0 ; i_uses < n_uses ; i_uses++){
    if ( uses[i_uses].useid == id ){ /*** ahh there is one ***/

      /* check buffersize  large enough */
      /* 2* -> save space for USES and \n */
      while ( (bused + 2*RDMC_MAXTOKENLENGTH) > bsize){ 
	bsize += RDMC_MAXLINE;
	buffer =  (char *) realloc(buffer,bsize*sizeof(char));
      }

      /* check the current lenght */
      if ( ( (char_in_line+RDMC_MAXTOKENLENGTH) >= F2000_MAXLINE )
	   || (token_in_line >= RDMC_MAXTOKEN_PER_LINE ) ){
	if (cont_flag){
	  bused += sprintf(&(buffer[bused]),"-%i",last_hid);
	  cont_flag=0;
	}
	bused += sprintf(&(buffer[bused]),"\n");
	char_in_line = token_in_line = 0 ;
	last_hid =  RDMC_NA;
	cont_flag=0;
      }
      if (char_in_line == 0 ){
	char_in_line = strlen("USES");
	bused += sprintf(&(buffer[bused]),"USES");
      }

      /* now print this channel */
      if ((last_hid == RDMC_NA ) ||
	  (uses[i_uses].hitid > last_hid+1 )){
	if (cont_flag){
	  t_len = sprintf(&(buffer[bused]),"-%i",last_hid);
	  last_hid=RDMC_NA;
	  bused += t_len;
	  char_in_line += t_len;
	  token_in_line++;
	  cont_flag=0;
	}
	t_len = sprintf(&(buffer[bused])," %i",uses[i_uses].hitid);
	last_hid=uses[i_uses].hitid;
	bused += t_len;
	char_in_line += t_len;
	token_in_line++;
	cont_flag=0;
      }else{
	last_hid=uses[i_uses].hitid;
	cont_flag=1;
      }
	
    } else if ( uses[i_uses].useid > id ){
#if 1
      break; /* break to improve speed  since the uses array is sorted */
#endif
    }
  }

  if (bused){
    if(buffer[bused-1] != '\n'){
      if (cont_flag){
	bused += sprintf(&(buffer[bused]),"-%i",last_hid);
	cont_flag=0;
      }
      strcpy(&(buffer[bused]),"\n");
      bused++;
    }
  }

#if 0 /* desireble later */
  if (buffer[0] == '\0')
    strcpy(buffer,"USES all\n");
#endif
  return buffer;
}


/****************************************************************************
 ********************************** E O F ***********************************
 ****************************************************************************/









