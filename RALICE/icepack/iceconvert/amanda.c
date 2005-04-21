/*
 * Read/write files in the oncoming Standard AMANDA format - "format 2000"
 */
#if 0 /******* TODO *********/
/* skip event gets array as parameter   */
/* level 4 event types in mevt  */
/*   end of file, tend          */
/******************************/        
#endif


char *amanda_rdmc_cvsid = 
"$Header: /net/local/cvsroot/siegmund/rdmc/amanda.c,v 1.81 2004/02/19 17:10:08 wiebusch Exp $";

#include <stdlib.h>
#include <ctype.h>

#include "rdmc.h"
#ifdef AMANDA_ASCII_F

#include "rdmc_local.h"

#include "amanda.h"
#include "f2k.h"

/****************************************************************************/
/*       rdmc_nint(double a)                                                 */
/* converts a to the nearest integer                                        */
/* but also takes care on the sign                                          */
/****************************************************************************/
long rdmc_nint(double a)
{ long i;
  i =  (a >= 0.) ?  (long) (a + 0.5 ) : (long) (a - 0.5);
 return  i;
}

/*******************************/
/* Work around missing isnan() */
/*******************************/
int isnan(float r)
{
 return 0;
}

/******************************************************/
/*** Definitions of f2k blocks struture in memory   ***/
/******************************************************/


/* check if line is in  buffer or read a new one */
/* returns 0 if OK */
static int amanda_readline(mcfile *fp);

/* dump a line back into the buffer if it is empty */
/* returns 0 if OK */
static int amanda_unreadline(mcfile *fp);

/* read the next block according to a list of allowed lines into the buffer */
/* returns 0 if OK */
static int amanda_read_block(mcfile *fp ,const f2000_event_t *def);

/* try to decode a given line according to line_def */
/* store if OK */
/* returns RDMC_IO_OK or RDMC_LINE_EMPTY or RDMC_EVENT_NOT_RECOGNIZED */ 
static int amanda_lex_line(mcfile *fp, const f2000_line_t * const line_def[F2K_MAX_LINETYPES] );

/* compares each byte of key to the first no white exam until lngth of key*/
/* -1 if line is empty, 0 on no and 1 on succes */
static int tokencmp(const char *key, char *exam);

/* same with no trailing chars */ 
static int tokencmp_notrail(const char *key, char *exam);

/* checks if the first non-white char in exam is in keys */
/* if a line is only whitespace ot empty -1 is returned */
/* 0: if the char is not found */
/* 1: if the char is not found */
static int tokenchr(const char *keys, char *exam);
/* same as above but chars are not alnum (ispunct()) */
static int tokenchrpunct(const char *keys, char *exam);

/***************************************************************************/
/***************************************************************************/
/***                                                                     ***/
/*** Read  functions follow                                              ***/
/***                                                                     ***/
/***************************************************************************/

  
/****************************************************************************
 * read the F2000 "V" line (the first one) and put it into fp
 ****************************************************************************/
/******* This is the only parswer routine, which is different !!! */
/* it reads directly from fp, instead of the readline mode */
int rdmc_amanda_V(const char *s, mcfile *fp)
{
  int major = 0, minor = 0;
  if(sscanf(s,"V 2000.%i.%i",&major, &minor) != 2) 
    return RDMC_LINE_NOT_PARSED;
  if (rdmc_is_f2000_supported(major, minor) == 0)
    return  RDMC_UNKNOWN_FORMAT_VERSION;
  fp->fmajor = major;
  fp->fminor = minor;
  return RDMC_IO_OK;

} /* amanda_rdmc_V() */

int rdmc_rhd_amanda(mcfile *fp){ 
  const  f2000_event_t *ev_def;
  int ret;

  /* get the config def for the right format*/
  switch( 100*fp->fmajor + fp->fminor ){
  case 101:
  case 102:
  case 200401:
    ev_def=&f2000_preamble_1x1;
    break;
  default:
    return RDMC_UNKNOWN_FORMAT_VERSION;
  }

  /* read the block acordng to the def */
  if ((ret = amanda_read_block(fp, ev_def)) != RDMC_IO_OK)
    return ret;

  /* use the pointer from  the event definition and scan*/
  return ev_def->reader(fp, NULL, NULL);
}


int rdmc_rarr_amanda(mcfile *fp, array *ar)
{ 
  const  f2000_event_t *ev_def;
  int ret;

  /* get the config def for the right format*/
  switch( 100*fp->fmajor + fp->fminor ){
  case 101:
  case 102:
  case 200401:
    ev_def=&f2000_mhead_1x1;
    break;
  default:
    return RDMC_UNKNOWN_FORMAT_VERSION;
  }

  /* read the block acordng to the def */
  if ((ret = amanda_read_block(fp, ev_def)) != RDMC_IO_OK)
    return ret;

  /* use the pointer from  the event definition and scan*/
  return ev_def->reader(fp, ar, NULL);

}

int rdmc_revt_amanda(mcfile *fp, mevt *ev, array *ar){ 
  const  f2000_event_t **ev_def,*foot_def;

  int ret = RDMC_ILF;
  int ret2 = RDMC_ILF;

  /* get the config def for the right format*/
  /* get the config def for the right format*/
  switch( 100*fp->fmajor + fp->fminor ){
  case 101:
    ev_def = f2000_events_1x1;
    foot_def =  &f2000_mfoot_1x1 ;
    break;
  case 102:
    ev_def =f2000_events_1x2;
    foot_def =  &f2000_mfoot_1x1 ;
    break;
  case 200401:
    ev_def =f2000_events_2004x1;
    foot_def =  &f2000_mfoot_1x1 ;
    break;
  default:
    return RDMC_UNKNOWN_FORMAT_VERSION;
  }

  while( *ev_def !=NULL ){
    /* read the block acordng to the def */
    ret = amanda_read_block(fp, *ev_def);
    if ( ret == RDMC_IO_OK)
      break;
    ++ev_def;
  }
  if ( ret == RDMC_IO_OK){
    rdmc_clear_mevt(ev);                           /* reset the event */
    /* use the pointer from  the event definition and scan*/
    if (fp->info.f2000.nolex)
      return RDMC_IO_OK;
    else
      return (*ev_def)->reader(fp, ar, ev);
  }else{  /* test footer */
    if ( ret == RDMC_EVENT_NOT_RECOGNIZED) {
      ret2 = amanda_read_block(fp, foot_def);
      if (ret2 == RDMC_IO_OK)
	return RDMC_EOF;
      else
	return ret2;
    }else{
      return ret;
    }
  }
}

int rdmc_skipevt_amanda(mcfile *fp)
{ 
  int ret;
#if 0
  array ar;
#endif
  mevt ev;
#if 0 /* this is the proper, presently slow version */
  rdmc_init_array(&ar); 
  rdmc_init_mevt(&ev);
  fp->info.f2000.nolex=1;
  ret = rdmc_revt_amanda(fp, &ev, &ar);
  fp->info.f2000.nolex=0;
  rdmc_clear_array(&ar);
  rdmc_clear_mevt(&ev);
#else /* presently no array is needed -> dirty patch ! */
  rdmc_init_mevt(&ev);
  fp->info.f2000.nolex=1;
  ret = rdmc_revt_amanda(fp, &ev, NULL);
  fp->info.f2000.nolex=0;
  rdmc_clear_mevt(&ev);
#endif
  return ret;
}



/****************************************************************************
 * Reads a nonempty line 
 ****************************************************************************/
static int amanda_readline(mcfile *fp){

  char *s = fp->last_line;
  do {
    if (fp->info.f2000.unread != 0 ) { /* line in buffer */
      fp->info.f2000.unread = 0;
    } else{ /* no line in buffer - read from file */
      if (fgets(s, RDMC_MAXLINE-1, fp->fp) == NULL)/*fgets includes the '\n'*/
	return EOF;
    }
  } while (*s == '\0');      /* look for a not-empty line */
  /* increment the counter line */
  fp->fpos++;
  return RDMC_IO_OK;

} /* amanda_readline() */

/****************************************************************************
 * "unread" an amanda F2000 format line
 ****************************************************************************/

static int amanda_unreadline(mcfile *fp){
  if (fp->info.f2000.unread != 0) 
    return EOF; /* only one line allowed */
  else
    fp->info.f2000.unread = 1;
  fp->fpos--;
  return RDMC_IO_OK;
} /* amanda_unreadline() */


static int amanda_read_block(mcfile *fp ,const f2000_event_t *def){
  int r;
  int ret=RDMC_ILF;
  int do_reading;
  rdmc_clear_f2k_buffer(fp->info.f2000.f2k_buffer);
  fp->errline = 0;

  /* test if the block starts with the right line */
  do_reading=1;
  do {
    if (def->opener[0]){  /* yes there is one needed ! */
      if ( (r=amanda_readline(fp)) != RDMC_IO_OK)
	ret=r;
      else /* analyse if this line matches and store it then */
	ret = amanda_lex_line(fp, def->opener);
    } else{ /* no opener needed  -> OK */
      ret = RDMC_IO_OK;
      break;
    }
    if( ret != RDMC_LINE_EMPTY ){ /* something happende */
      do_reading=0;
      if ( ret == RDMC_EVENT_NOT_RECOGNIZED){
	amanda_unreadline(fp); /* we have read one line to much  */
	return ret;
      }
    }
  } while(do_reading);

  /* keep on reading body lines  */
  do_reading=1;
  do {
    /* analyse if this line matches and store it then */
    if (def->inner[0]){  /* yes there are some  ! */
      if ( (r=amanda_readline(fp)) != RDMC_IO_OK)
	ret =r;
      else
	ret = amanda_lex_line(fp, def->inner);
    } else{ 
      ret = RDMC_IO_OK;
      do_reading=0;
      break;
    }

    if( ret != RDMC_LINE_EMPTY ){ /* something happend */
      if (ret == RDMC_IO_OK) /* line OK -> next one */
	continue;
      else{ /* line is not parsed -> check if an end marker appeared */
	amanda_unreadline(fp); /* we have read one line to much  */
	do_reading = 0; /* this is no error since the body may be empty */
	ret = RDMC_IO_OK; /* this is no error but maybe  the next event */
	break; /* ok we end here */
      }
    }
  } while(do_reading);

  /* check the end marker   */
  do_reading=1;
  do {
    /* analyse if this line matches and store it then */
    if (def->closer[0]){  /* yes there are some needed ! */
      if ( (r=amanda_readline(fp)) != RDMC_IO_OK)
	ret=r;
      else
	ret = amanda_lex_line(fp, def->closer);
    } else{ 
      ret = RDMC_IO_OK;
      do_reading=0;
      break;
    }

    if( ret != RDMC_LINE_EMPTY ){ /* something happend */
      if (ret == RDMC_IO_OK){ /* line OK -> finish */
	do_reading=0;
	break;
      }
      else{ /* line is not parsed  */
	amanda_unreadline(fp); /* we have read one line to much  */
	do_reading = 0; /* this is no error since the body may be empty */
	break; /* ok we end here */
      }
    }
  } while(do_reading);

  return ret;
}


static int amanda_lex_line(mcfile *fp, const f2000_line_t * const line_def[F2K_MAX_LINETYPES] ){
  int ptindex , cres = 0;
  const f2000_line_t *opt;

  if ( !(line_def[0]) ) /* no match actually needed */
    return RDMC_EVENT_NOT_RECOGNIZED;

  for(ptindex=0 , cres=0, opt = line_def[0] 
	; opt 
	; opt = line_def[++ptindex] 
      ){

    switch (opt->searchtype){
    case COMP_STRINGWISE:
      cres = tokencmp(opt->tag,fp->last_line);
      break;
    case COMP_STRINGWISE_NOTRAIL:
      cres = tokencmp_notrail(opt->tag,fp->last_line);
      break;
    case COMP_CHARWISE:
      cres = tokenchr(opt->tag,fp->last_line);
      break;
    case COMP_CHARPUNCT:
      cres = tokenchrpunct(opt->tag,fp->last_line);
      break;
#if 0
    case COMP_DUMMY: /* no token is needed so stop parsing */
      /* trap to finish here */
      return RDMC_EVENT_NOT_RECOGNIZED;
      break;
#endif
    default:
      return RDMC_LIBRARY_ERROR;
    } /* switch */
    if (cres < 0){ /* empty line */
      return RDMC_LINE_EMPTY;
    } else if(cres == 0){ /* no match found */
      /* needed but not found -> try the next */
      continue;
    }else{ /* ok tag is found */
      rdmc_push_f2k_buffer(fp->info.f2000.f2k_buffer, fp->last_line,line_def[ptindex]); 
      return RDMC_IO_OK;
    } /* check cres */
  } /*for */
  return RDMC_EVENT_NOT_RECOGNIZED;
}

/* compares each byte of key to exam until th lngth of key */
/* leading whitespaces are ignored */
static int tokencmp(const char *key, char *exam){
  int r = -1;
  const char *k = key;
  const char *e = exam;
  while(*e){
    if (isspace(*e))
      ++e;
    else { /* this is the first non white */
      r=0;
      while (*k){ /* now check key */
	if (!(*e)) /* if e ends before key */ 
	  return r=0;
	if( *e == *k ) /* stil agreement ? */
	  r=1;
	else
	  return r=0; /* break if not */
	++e; /* next char */
	++k;
      } /* while k */
      return r; /*  exit here likely with r==1  */
    } 
  } /* while e */

  return r;
}
/* compares each byte of key to exam until the lngth of key */
/* trailing chars are ignored */
static int tokencmp_notrail(const char *key, char *exam){
  int r = -1;
  const char *k = key;
  const char *e = exam;
  while(*e){
    if (isspace(*e))
      ++e;
    else { /* this is the first non white */
      r=0;
      while (*k){ /* now check key */
	if (!(*e)) /* if e ends before key */ 
	  return r=0;
	if( *e == *k ) /* stil agreement ? */
	  r=1;
	else
	  return r=0; /* break if not */
	++e; /* next char */
	++k;
      } /* while k */
      return r; /*  exit here likely with r==1  */
    } 
  } /* while e */
  return r;
}

/* checks if the first non-white char in exam is in keys */
/* if a line is only whitespace ot empty -1 is returned */
/* 0: if the char is not found */
/* 1: if the char is found */
static int tokenchr(const char *keys, char *exam){
  int r = -1;
  const char *k=keys;
  const char *e=exam;
  while(*e){
    if (isspace(*e))
      ++e;
    else { /* this is the first non white */
      r=0;
      while (*k){ /* now check key */
	if( *e == *k ) /* agreement ? */
	  return r=1;
	else
	  ++k;
      } /* while k */
      return r=0; /*  exit here likely with r==1  */
    } 
  } /* while e */
  return r;
}

static int tokenchrpunct(const char *keys, char *exam){
  int r = -1;
  const char *k=keys;
  const char *e=exam;
  while(*e){
    if (isspace(*e))
      ++e;
    else if (!ispunct(*e))
      return r=0;
    else { /* this is the first non white */
      r=0;
      while (*k){ /* now check key */
	if( *e == *k ) /* agreement ? */
	  return r=1;
	else
	  ++k;
      } /* while k */
      return r=0; /*  exit here likely  */
    } 
  } /* while e */
  return r;
}


/***************************************************************************/
/***                                                                     ***/
/*** Write functions follow                                              ***/
/***                                                                     ***/
/***************************************************************************/
/***************************************************************************/




/***************************************************************************/
/* write end of file (F2000 format)                                        */
int rdmc_wrend_amanda(const mcfile *fp)
{
  const  f2000_event_t *ev_def;

  /* get the config def for the right format*/
  switch( 100*fp->fmajor + fp->fminor ){
  case 101:
  case 102:
  case 200401:
    ev_def=&f2000_mfoot_1x1;
    break;
  default:
    return RDMC_UNKNOWN_FORMAT_VERSION;
  }

  /* use the pointer from  the event definition and scan*/
  return ev_def->writer(fp, NULL, NULL);
} /* wrend_amanda() */

/****************************************************************************
 * function warr_amanda() writes the array info to a amanda-like file
 * This function writes out the head of a AMANDA ascii file
 * opposite to reading the input file it writes not only
 * the Geometry banks ('G', 'P') , but also the ('V' and 'M' flags)
 * so the function whd_amanda does not exist
 ****************************************************************************/

int rdmc_warr_amanda(const mcfile *fp,const array *geo)
{
  const  f2000_event_t *ev_def;

  /* get the config def for the right format*/
  switch( 100*fp->fmajor + fp->fminor ){
  default:
    ev_def=&f2000_mhead_1x1;
    break;
  }

  /* use the pointer from  the event definition and scan*/
  return ev_def->writer(fp, geo, NULL);
}

/****************************************************************************
 * function wevt_amanda() writes an event to a amanda-like file
 ****************************************************************************/

int rdmc_wevt_amanda(const mcfile *fp,const mevt *event, const array *ar)
{
  const  f2000_event_t *ev_def;

  /* get the config def for the right format*/
  switch( 100*fp->fmajor + fp->fminor ){
  case 101:
    ev_def=f2000_events_1x1[0]; /* take first array element */
    break;
  case 102:
    ev_def=f2000_events_1x2[0]; /* take first array element */
    break;
  case 200401:
    ev_def=f2000_events_2004x1[0]; /* take first array element */
    break;
  default:
    return RDMC_UNKNOWN_FORMAT_VERSION;
  }

  /* use the pointer from  the event definition and scan*/
  return ev_def->writer(fp, ar, event);

}

#if 0
/****************************************************************************
 * wrhist_amanda() writes history lines to an ASCII file
 ****************************************************************************/

int rdmc_wrhist_amanda(const mcfile *fp, const char *s, const char *pre)
{
  return rdmc_wrhist_f2k_1x1(fp, s, pre);
}

/****************************************************************************
 * wrcomment_amanda() writes a comment line to an ASCII file
 ****************************************************************************/

int rdmc_wrcomment_amanda(const mcfile *fp, const char *s)
{
  return rdmc_wrcomment_f2k_1x1(fp, s);
}
#endif

#endif /* AMANDA_ASCII_F */

/****************************************************************************
 ********************************** E O F ***********************************
 ****************************************************************************/
/* 
   This is just for EMACS:
   Local Variables:
   compile-command: "cd .. ; make -k rdmc" 
   End: 
*/












