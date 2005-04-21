/****************************************************************************/
/* rdmc_mcopen() opens a mc/data file                                       */
/* usage is like fopen(), format is BAIKAL_BIN_F for the baikal like binary */
/* format, and DUMAND_ASCII_F for the DUMAND-like ascii format              */
/****************************************************************************/

#undef RDMC_MCOPEN_DEBUG

#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
///////////////#include <unistd.h>
#include <sys/stat.h>

#if defined(OSF1) || defined (SunOS) || defined(IRIX)
#include <alloca.h>
#endif


#include "rdmc.h"
#include "rdmc_local.h"


#include "baikal.h"
#include "amanda.h"
#include "uwi.h"
#include "dumand.h"


FILE *rdmc_gzipopen(const char *name, const char *mode);                  /* open a gzip file */
FILE *rdmc_bzipopen(const char *name, const char *mode);                  /* open a bzip file */
FILE *rdmc_bzip2open(const char *name, const char *mode);                /* open a bzip2 file */
int rdmc_strsuf(const char *str, const char *suf);           /* test for a suffix */
FILE *rdmc_myopen(const char *name, const char *mode, int *piped);/*open a (compressed) file */
char *rdmc_remote_pipe(const char *name, const char *mode);  /* builds pipe command for rsh  */
char *rdmc_is_in_path(const char *program); /* is certain program in the path. */



mcfile *rdmc_mcopen(const char *name, const char *mode, int format)
{
  mcfile *fp;
  FILE *tf; /* temporary file pointer */
  int r;
  int fmode,zipped; /* read or write, and detected format */
  char s[RDMC_MAXLINE];

  /* test if format is allowed */
  if ( 1
#ifdef DUMAND_ASCII_F
      && (format != DUMAND_ASCII_F)
#endif
#ifdef UWI_ASCII_F
      && (format != UWI_ASCII_F)
#endif
#ifdef AMANDA_ASCII_F
      && (format != AMANDA_ASCII_F)
#endif
#ifdef BAIKAL_BIN_F
      && (format != BAIKAL_BIN_F ) 
#endif
       ){
#ifdef RDMC_MCOPEN_DEBUG
     rdmc_warnprintf("DEBUG: Unknown format %i", format);
#endif
     return (NULL);
  }

  /* now get the mode to open the file */
  if (strstr(mode,"r") != NULL){ /* this is reading */
    fmode = RDMC_READ_MODE;
  }  else if (strstr(mode,"w") != NULL){ /* this is writing */
    fmode = RDMC_WRITE_MODE;
  }else{ /* this is nonesensce */
#ifdef RDMC_MCOPEN_DEBUG
     rdmc_warnprintf("DEBUG: Unknown mode %s (not r or w)", mode);
#endif
    return (NULL);   
  }    

  /* allocate memory for the mcfile structure and initilize to 0*/
  if ((fp = (mcfile *)calloc(1,sizeof(mcfile))) == NULL){
#ifdef RDMC_MCOPEN_DEBUG
     rdmc_warnprintf("DEBUG: Could not allocate mcfile structure");
#endif
    return (NULL);   
  }

  if (name[0] != '\0') {                                 /* if no stdin/out */
    tf = rdmc_myopen(name, mode, &(zipped));
    if (tf == NULL) {                                 /* if open failed */
      free(fp);
#ifdef RDMC_MCOPEN_DEBUG
      rdmc_warnprintf("DEBUG: could not open %s (zipped=%i)",name, zipped);
#endif
      return (NULL);
    } /* if fopen == 0 */
  } /* if name != "" */
  else {                                           /* no name --> stdin/out */
    if (tolower((int) mode[0]) == 'w')                  /* I will write to stdout */
      tf = stdout;
    else {        /* I will read from stdin (currently only *_ASCII_F) */
      if (0
#ifdef DUMAND_ASCII_F
	  || (format == DUMAND_ASCII_F) 
#endif
#ifdef AMANDA_ASCII_F
	  || (format == AMANDA_ASCII_F) 
#endif
#ifdef UWI_ASCII_F
	  || (format == UWI_ASCII_F)
#endif
	  )
	tf = stdin;
      else {
	free(fp);
#ifdef RDMC_MCOPEN_DEBUG
        rdmc_warnprintf("DEBUG: Unknown format %i to read from stdin", format);
#endif
	return(NULL);
      } /* if no 'w' and no ASCII_F */
    } /* in no 'w' */
  } /* if name == "" */


  /* before we can do something with the file 
     we need to track down the exact format for the ASCII format, 
     we try automatic detection in case of file-reading
      dumand starts with "V"
      while UWI starts with "F", 
      F2000 starts with "V 2000"
   */
  if (fmode == RDMC_READ_MODE){
    if (0
#ifdef DUMAND_ASCII_F
	|| (format == DUMAND_ASCII_F)
#endif
#ifdef UWI_ASCII_F
	|| (format == UWI_ASCII_F) 
#endif
#ifdef AMANDA_ASCII_F
	|| (format == AMANDA_ASCII_F)
#endif
	){
      if (fgets(s,RDMC_MAXLINE,tf) == NULL) {  /* try to read first line */
	fclose(tf);
	free(fp);                                         /* free the memory */
#ifdef RDMC_MCOPEN_DEBUG
        rdmc_warnprintf("DEBUG: Cannot read first line");
#endif
	return NULL;
      } else {
	r=1;
	/* ok we have the line in s and have to decode it now */
#ifdef AMANDA_ASCII_F
	if (r) {
	  rdmc_init_mcfile(fp,AMANDA_ASCII_F,fmode,tf);
	  r=rdmc_amanda_V(s,fp);
	  if(r == 0){
	    fp->format = AMANDA_ASCII_F;
#if OLD
	    rdmc_push_line_buffer(&(fp->h_buff),s);
#endif
	    fp->fpos++;
	    strncpy(fp->last_line,s,RDMC_MAXLINE-1); 
	    fp->last_line[RDMC_MAXLINE]='\0';
	    fp->info.f2000.unread=1;
	  }
	}
#endif
#ifdef DUMAND_ASCII_F
	if (r) {
	  rdmc_init_mcfile(fp,DUMAND_ASCII_F,fmode,tf);
	  r= rdmc_a_V(s,fp);
	  if ( r == 0){             /* scan for version number */
	    fp->format = DUMAND_ASCII_F;
	    fp->fpos++;
	    strncpy(fp->last_line,s,RDMC_MAXLINE-1); 
	    fp->last_line[RDMC_MAXLINE]='\0';
	  }
	}
#endif
#ifdef UWI_ASCII_F
	if (r) {
	  rdmc_init_mcfile(fp,UWI_ASCII_F,fmode,tf);
	  r = rdmc_uwi_FH(s,fp);
	  if ( r == 0){
	      fp->format = UWI_ASCII_F;
	      fp->fpos++;
	      strncpy(fp->last_line,s,RDMC_MAXLINE-1); 
	      fp->last_line[RDMC_MAXLINE]='\0';
	  }
	}
#endif
	if (r) { /* file could not be deytected */
	  fclose(tf);
	  free(fp);                                   /* free the memory */
#ifdef RDMC_MCOPEN_DEBUG
          rdmc_warnprintf("DEBUG: Cannot recognize format of the following line:\n%s", s);
#endif
	  return NULL;
	}
      }
    } /* format defined */
    else{ /* e.g. Baikal format */
      rdmc_init_mcfile(fp,format,fmode,tf);
    }
  }else{ /* opened for writing */
    rdmc_init_mcfile(fp,format,fmode,tf);
  }/*atodetect in case of reading */
  
  /* update the zip status */
  fp->zipped=zipped;
   
  if ((strstr(mode,"r") == NULL) && (strstr(mode,"a") == NULL)) 
    return fp;                                  /* only write mode: no test */
  
  /* Now read the history and comment lines until the header starts */
  switch(fp->format) {
#ifdef UWI_ASCII_F
  case UWI_ASCII_F:
    r = rdmc_rhd_uwi(fp);                            /* read the UWI-like header */
    break;
#endif
#ifdef DUMAND_ASCII_F
  case DUMAND_ASCII_F:
    r = rdmc_rhd_ascii(fp);                       /* read the dumand-like header */
    break;
#endif
#ifdef AMANDA_ASCII_F
  case AMANDA_ASCII_F:
    r = rdmc_rhd_amanda(fp);                       /* read the amanda-like header */
    break;
#endif
#ifdef BAIKAL_BIN_F
  case BAIKAL_BIN_F:                               /* for the binary format */
    r = rdmc_rhd_baikal_mc(fp);               /* read the baikal like header */
    break;
#endif
  default:
    r = RDMC_UNKNOWN_FORMAT;
  } /* switch fp->format */

  if (r != 0) {                    /* if header reading was not successfull */
    free(fp);                                            /* free the memory */
#ifdef RDMC_MCOPEN_DEBUG
     rdmc_warnprintf("DEBUG: Cannot read header (%i)",r);
#endif
    fp = NULL;                             /* and reset the pointer to NULL */
  }

  return (fp);
  
} /* function rdmc_mcopen() */

/****************************************************************************/
/* mcclose() closes a mc/data file                                          */
/* usage is like fclose()                                                   */
/****************************************************************************/

int rdmc_mcclose(mcfile *fp)
{
  int c=1;

#ifdef AMANDA_ASCII_F
  /* write end of file mark */
  if (fp->format == AMANDA_ASCII_F){
    rdmc_wrend_amanda(fp);
  }
#endif
  
  /* close the file */
  if (fp->zipped == 0)
    c = fclose(fp->fp);
  else
//    c = pclose(fp->fp);

  rdmc_free_mcfile(fp);   /* free/unallocate file structure */

  return(c);

} /* function mcclose() */


/*****************************************************************************/
/* rdmc_myopen() opens a file, maybe zipped. It returns the file pointer and      */
/*          a flag if it has opened a pipe                                   */
/*****************************************************************************/
FILE *rdmc_myopen(const char *name, const char *mode, int *piped)
{
  typedef FILE*(fopen_like_t)(const char *name, const char *mode);
  struct open_table_s {
    char *suffix;
    fopen_like_t *open_function;
  } open_table[] = {
    {".gz", rdmc_gzipopen},
    {"-gz", rdmc_gzipopen},
    {".GZ", rdmc_gzipopen},
    {"-GZ", rdmc_gzipopen},
    {".bz", rdmc_bzipopen},
    {"-bz", rdmc_bzipopen},
    {".bz2", rdmc_bzip2open},
    {"-bz2", rdmc_bzip2open},
    {NULL, NULL}
  };
  struct open_table_s *entry;
   char *r;
   
  for (entry = open_table; entry->suffix != NULL; entry++)
    if (rdmc_strsuf(name,entry->suffix)) {
      *piped = 1;
//      return entry->open_function(name,mode);
    }

  if ((r = rdmc_remote_pipe(name, mode))) {
    *piped = 1;
//    return popen(r, mode);
  } /* if   file name */

  *piped = 0;
  return fopen(name, mode); /* default: normal open */

} /* rdmc_myopen() */

/*****************************************************************************/
/* strsuffix() tests if suf is a file suffix of str                          */
/*****************************************************************************/
int rdmc_strsuf(const char *str, const char *suf)
{
  const char *str_c, *suf_c;

  if (strlen(suf) > strlen(str))
    return 0;
  for (str_c = str+strlen(str), suf_c = suf+strlen(suf);
       (suf_c >= suf) && (str_c >= str); 
       str_c--, suf_c--)
    if (*suf_c != *str_c) return 0;

  return 1;

} /* strsuffix() */

/*****************************************************************************/
/* gzipopen() opens a gzipped file for reading or writing (no seek!)         */
/*    the file can be used with the usual fprintf()..pclose() functions      */
/*    the function needs the gzip program to work, and will not work under   */
/*    MS-DOS o.a.                                                            */
/*****************************************************************************/

FILE *rdmc_gzipopen(const char *name, const char *mode)
{
  char command[RDMC_MAXLINE];
  char *r;

  r = rdmc_remote_pipe(name, mode);
  if (strchr(mode,'w')) {           /* if there was a 'w' in the mode string */
    if (r)
      sprintf(command, "gzip | %s", r);
    else
      sprintf(command, "gzip > %s", name);
//    return popen(command,"w");
  }
  if (strchr(mode,'r')) {                        /* but if there was an 'r' */
    if (r)
      sprintf(command, "%s | gzip -dc", r);
    else
      sprintf(command, "gzip -dc %s", name);
//    return popen(command,"r");
  }

  return NULL;

} /* rdmc_gzipopen() */

/*****************************************************************************/
/* bzipopen() opens a zipped file for reading or writing (no seek!)          */
/*    the file can be used with the usual fprintf()..pclose() functions      */
/*    the function needs the bzip program to work, and will not work under   */
/*    MS-DOS o.a.                                                            */
/*****************************************************************************/

FILE *rdmc_bzipopen(const char *name, const char *mode)
{ 
  char command[RDMC_MAXLINE];
  char *r;

  r = rdmc_remote_pipe(name, mode);
  if (strchr(mode,'w')) {           /* if there was a 'w' in the mode string */
    if (r)
      sprintf(command, "bzip | %s", r);
    else
      sprintf(command, "bzip > %s", name);
//    return popen(command,"w");
  }
  if (strchr(mode,'r')) {                        /* but if there was an 'r' */
    if (r)
      sprintf(command, "%s | bzip -dc", r);
    else
      sprintf(command, "bzip -dc %s", name);
//    return popen(command,"r");
  }

  return NULL;

} /* bzipopen() */

/*****************************************************************************/
/* bzip2open() opens a zipped file for reading or writing (no seek!)         */
/*    the file can be used with the usual fprintf()..pclose() functions      */
/*    the function needs the bzip2 program to work, and will not work under  */
/*    MS-DOS o.a.                                                            */
/* (bzip2 is a replacement for bzip with better+faster compression and       */
/*  without patent problems)                                                 */
/*****************************************************************************/

FILE *rdmc_bzip2open(const char *name, const char *mode)
{
  char command[RDMC_MAXLINE];
  char *r;

  r = rdmc_remote_pipe(name, mode);
  if (strchr(mode,'w')) {           /* if there was a 'w' in the mode string */
    if (r)
      sprintf(command, "bzip2 | %s", r);
    else
      sprintf(command, "bzip2 > %s", name);
//    return popen(command,"w");
  }
  if (strchr(mode,'r')) {                        /* but if there was an 'r' */
    if (r)
      sprintf(command, "%s | bzip2 -dc", r);
    else
      sprintf(command, "bzip2 -dc %s", name);
//    return popen(command,"r");
  }

  return NULL;

} /* bzip2open() */

/*
 * is_in_path() checks if a certain program can be found in the path.
 * return value is the full path (with program name), or NULL
 */
char *rdmc_is_in_path(const char *program)
{
  char *path ,*p1;
  static char *p=NULL;
  static char progpath[RDMC_MAXLINE];
  struct stat statbuf;

  path = getenv("PATH");
#ifndef CRAY  
  p = alloca(strlen(path)+1);
#else
  p = realloc(p,strlen(path)+1);
#endif
  strcpy(p, path);

  for (p1 = strtok(p,":"); p1 != NULL; p1 = strtok(NULL,":")) {
    sprintf(progpath,"%s/%s",p1,program);
    if (!stat(progpath, &statbuf)) return progpath;
  } /* for p1 */

  return NULL;

} /* is_in_path() */

/*
 * remote_pipe() builds the part of a pipe command for
 * remote file opening. 
 * Just add this string to the beginning or end of your openng command.
 * Actually, it looks just for a ':' in the file
 * name and interprets the part before it as host (to be reached via RSH_CMD)
 * if no remote part is found, NULL will be returned.
 */

char *rdmc_remote_pipe(const char *name, const char *mode)
{
#ifdef ALLOW_REMOTE
  static char command[RDMC_MAXLINE];
  char host[RDMC_MAXLINE];
  char *file;
  static char *rsh_cmd = NULL;

  file = strchr(name,':');
  if (file == NULL) return NULL;

  if (!rsh_cmd) { /* init rsh command */
    rsh_cmd = getenv(RSH_CMD_ENV);
    if (!rsh_cmd) {
      if (rdmc_is_in_path("ssh"))
	rsh_cmd = "ssh";
      else if (rdmc_is_in_path("rsh"))
	rsh_cmd = "rsh";
      else if (rdmc_is_in_path("remsh"))
	rsh_cmd = "remsh";
      else rsh_cmd = NULL;
    }
  }
  if (!rsh_cmd) return NULL;

  strncpy(host, name, file-name);
  if (strchr(mode,'r'))                         /* but if there was an 'r' */
    sprintf(command, "%s %s cat %s", 
	    rsh_cmd, host, file+1);
  else if (strchr(mode,'w'))      /* if there was a 'w' in the mode string */
    sprintf(command, "%s %s \"cat > %s\"", 
	    rsh_cmd, host, file+1);
  else return NULL;
  return command;
#else
  return NULL
#endif
   
} /* remote_pipe() */
