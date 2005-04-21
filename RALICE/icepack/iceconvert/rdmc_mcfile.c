/*
 *  functions for the mcfile structure 
 */

#include <stdlib.h>
#include <string.h>

#include "rdmc.h"
#include "rdmc_local.h"

#include "f2k.h"

void rdmc_init_mcfile_f2000(f2000_fp_t *t);
void rdmc_free_mcfile_f2000(f2000_fp_t *t);
void rdmc_init_mcfile_dumand(dumand_fp_t  *t);
void rdmc_init_mcfile_uwi(uwi_fp_t   *t);
void rdmc_init_mcfile_baikal(baikal_fp_t   *t);


/****************************************************************************/
/* mccpfp() copies the relevant fields from one mcfile* to another          */
/* (this is useful for file transforming)                                   */
/****************************************************************************/

int rdmc_mccpfp(const mcfile *src, mcfile *dest)
{
  if ((src == NULL) || (dest == NULL))
    return -1;
  /* do not copy 
     the file pointer, src->fp 
     the mode, src->mode
     the format, src->format 
  also not: fpos,fmajor,fminor,zipped,errline,sloppy, last_line h_buff,e_buff
  But: creator,comment,nrun,time
  The format specific stuff is converted to comments */
#if 0
  dest->time = src->time;
  dest->nrun = src->nrun;
#endif
  memcpy(dest->last_line,src->last_line,RDMC_MAXLINE*sizeof(char));

  /* copy creator + comment structure */
  if( src->creator)
    rdmc_append_comment(&(dest->creator),src->creator);
  if (src->comment)
    rdmc_append_comment(&(dest->comment),src->comment);

  /* now tryy to convert the special info */
  switch (src->format){
#ifdef AMANDA_ASCII_F 
  case AMANDA_ASCII_F :
    /* nothing is to be copied here */
    break;
#endif
#ifdef DUMAND_ASCII_F
  case DUMAND_ASCII_F:       /* if it is a dumand-like format */
    switch (dest->format){
    case DUMAND_ASCII_F:
      memcpy(&(dest->info),&(src->info),sizeof(rdmc_special_format_info_t));
      break;
#ifdef AMANDA_ASCII_F 
    case AMANDA_ASCII_F :
      {     /* copy everything to a comment line */
	char tmp[RDMC_MAXLINE];
	sprintf(tmp,"! DUM-Form M-header: id=%c mc_vers=%s \n",src->info.dum.mc_id,src->info.dum.mc_vers);
	rdmc_append_comment(&(dest->comment),tmp);
	sprintf(tmp,
		"! DUM-Form M-header: igen=%i, igeo=%i, igtrack=%i, daswrun=%i\n"
		,src->info.dum.igen,src->info.dum.igeo
		,src->info.dum.igtrack,src->info.dum.daswrun);
	rdmc_append_comment(&(dest->comment),tmp);
      }
      break;
#endif
#ifdef UWI_ASCII_F 
    case UWI_ASCII_F :
      dest->info.uwi.mc_id = src->info.dum.mc_id;
      strcpy(dest->info.uwi.mc_vers, src->info.dum.mc_vers);
      break;
#endif
#ifdef BAIKAL_BIN_F 
    case BAIKAL_BIN_F :
      dest->info.bai.mc_id = src->info.dum.mc_id;
      dest->info.bai.igen = src->info.dum.igen;
      dest->info.bai.igtrack = src->info.dum.igtrack;
      dest->info.bai.time = src->info.dum.time;
      dest->info.bai.nrun = src->info.dum.nrun;
      strcpy(dest->info.bai.mc_vers, src->info.dum.mc_vers);
      break;
#endif
    default: 
      return RDMC_UNKNOWN_FORMAT;
    }
    break;
#endif
#ifdef UWI_ASCII_F
  case UWI_ASCII_F:       /* if it is the UWI format: read default geometry */
    switch (dest->format){
    case UWI_ASCII_F:
      memcpy(&(dest->info),&(src->info),sizeof(rdmc_special_format_info_t));
      break;
#ifdef AMANDA_ASCII_F 
    case AMANDA_ASCII_F :
      {     /* copy everything to a comment line */
	char tmp[RDMC_MAXLINE];
	sprintf(tmp,"! UWI-Format header: id=%c mc_vers=%s \n"
		,src->info.uwi.mc_id,src->info.uwi.mc_vers);
	rdmc_append_comment(&(dest->comment),tmp);
      }
      break;
#endif
#ifdef DUMAND_ASCII_F 
    case DUMAND_ASCII_F :
      dest->info.dum.mc_id = src->info.uwi.mc_id;
      strcpy(dest->info.dum.mc_vers, src->info.uwi.mc_vers);
      break;
#endif
#ifdef BAIKAL_BIN_F 
    case BAIKAL_BIN_F :
      dest->info.bai.mc_id = src->info.uwi.mc_id;
      strcpy(dest->info.bai.mc_vers, src->info.uwi.mc_vers);
      break;
#endif
    default: 
      return RDMC_UNKNOWN_FORMAT;
    }
    break;
#endif
#ifdef BAIKAL_BIN_F
  case BAIKAL_BIN_F:                       /* if it is a baikal-like format */
    switch (dest->format){
    case BAIKAL_BIN_F:
      memcpy(&(dest->info),&(src->info),sizeof(rdmc_special_format_info_t));
      break;
#ifdef AMANDA_ASCII_F 
    case AMANDA_ASCII_F :
      {     /* copy everything to a comment line */
	char tmp[RDMC_MAXLINE];
	sprintf(tmp,"! BAI-Form header: id=%c mc_vers=%s \n",src->info.bai.mc_id,src->info.bai.mc_vers);
	rdmc_append_comment(&(dest->comment),tmp);
	sprintf(tmp,"! BAI-Form header: igen=%i, igtrack=%i\n"
		,src->info.bai.igen,src->info.dum.igtrack);
	rdmc_append_comment(&(dest->comment),tmp);
	sprintf(tmp,"! BAI-Form header: rec_vers=%i, min_vers=%i mc=%i\n"
		,src->info.bai.rec_vers,src->info.bai.min_vers
		,src->info.bai.mc);
	rdmc_append_comment(&(dest->comment),tmp);
      }
      break;
#endif
#ifdef DUMAND_ASCII_F 
    case DUMAND_ASCII_F :
      dest->info.dum.mc_id = src->info.bai.mc_id;
      dest->info.dum.igen = src->info.bai.igen;
      dest->info.dum.igtrack = src->info.bai.igtrack;
      dest->info.dum.time = src->info.bai.time;
      dest->info.dum.nrun = src->info.bai.nrun;
      strcpy(dest->info.dum.mc_vers, src->info.bai.mc_vers);
      break;
#endif
#ifdef UWI_ASCII_F 
    case UWI_ASCII_F :
      dest->info.dum.mc_id = src->info.bai.mc_id;
      strcpy(dest->info.dum.mc_vers, src->info.bai.mc_vers);
      break;
#endif
    default: 
      return RDMC_UNKNOWN_FORMAT;
    }
    break;
#endif
  default: 
    return RDMC_UNKNOWN_FORMAT;
  }

  return 0; 
} /* mccpfp() */




void rdmc_init_mcfile(mcfile *fp, int format, int mode, FILE *tf)
{
  fp->fp = tf;
  fp->mode = mode;
  fp->format = format;
  fp->fmajor = fp->fminor = fp->zipped  = 0;
  fp->creator = fp->comment = NULL;            /* default: no file creator */
  fp->errline = 0;
  fp->fpos = 0;
#if 0
  fp->time = 0;
  fp->nrun = 0;
#endif
  fp->sloppy = 0;

  switch (format) {
#ifdef DUMAND_ASCII_F
  case DUMAND_ASCII_F:       /* if it is a dumand-like format */
    rdmc_init_mcfile_dumand(&(fp->info.dum));
    fp->fmajor = DUMAND_ASCII_VERSION; 
    break;
#endif
#ifdef UWI_ASCII_F
  case UWI_ASCII_F:       /* if it is the UWI format: read default geometry */
    rdmc_init_mcfile_uwi(&(fp->info.uwi));
    fp->fmajor = UWI_ASCII_VERSION;
    break;
#endif
#ifdef AMANDA_ASCII_F
  case AMANDA_ASCII_F:       /* if it is a amanda-like format */
    rdmc_init_mcfile_f2000(&(fp->info.f2000));
    fp->fmajor = AMANDA_ASCII_VERSION; 
    fp->fminor = AMANDA_ASCII_MINOR; 
    break;
#endif
#ifdef BAIKAL_BIN_F
  case BAIKAL_BIN_F:                       /* if it is a baikal-like format */
    rdmc_init_mcfile_baikal(&(fp->info.bai));
    fp->fmajor = BAIKAL_NEW;                         /* default: new format */
    break;
#endif
  default: 
    break;
  }
}

void rdmc_init_mcfile_f2000(f2000_fp_t *t){
  t->unread=0;
  t->parse_line[0] = '\0' ;
  t->errline=0;
  t->nolex=0;
  t->f2k_buffer = malloc(sizeof(rdmc_f2k_buffer_t));
  rdmc_init_f2k_buffer(t->f2k_buffer);
}

void rdmc_free_mcfile_f2000(f2000_fp_t *t){
  rdmc_unlink_f2k_buffer(t->f2k_buffer);
  free(t->f2k_buffer);
  t->f2k_buffer=NULL;
}

void rdmc_init_mcfile_dumand(dumand_fp_t  *t){
  t->mc_id = 'U';
  t->igtrack = -1;
  t->daswrun = 0;
  t->igeo = 0;
  t->igen = 0;
  t->time = 0;
  t->nrun = 0;
}

void rdmc_init_mcfile_uwi(uwi_fp_t   *t){
  t->mc_id = 'U';
  strcpy(t->mc_vers,"0.");
}

void rdmc_init_mcfile_baikal(baikal_fp_t   *t){
  t->mc = 1;                                /* default: mc file */
  t->swap = 0;                                /* default: no bytes swapped */
  t->fortran_recs = 1;                         /* default: fortran file */
  t->nw_rec_arr = 3;
  t->nw_rec_evt = 9;
  t->rec_vers = 0;
  t->min_vers = 0;
  t->enr = 0;                                        /* first event number */
  t->mc_id = 'U';
  strcpy(t->mc_vers,"0.");
  t->igtrack = -1;
  t->igen = 0;
  t->time = 0;
  t->nrun = 0;
}

void rdmc_free_mcfile(mcfile *fp)
{
  if (fp->format == AMANDA_ASCII_F)
    rdmc_free_mcfile_f2000(&(fp->info.f2000));

  if(fp->creator){
    free(fp->creator);
    fp->creator=NULL;
  }
  if(fp->comment){
    free(fp->comment);
    fp->comment=NULL;
  }
  free(fp);
}


