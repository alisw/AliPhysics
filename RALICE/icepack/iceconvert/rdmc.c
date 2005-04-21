/****************************************************************************/
/* This is rdmc.c - some routines for reading mc and preprocessed data      */
/* files in the formats described in format.txt of jk/rw (10.10.93)         */
/* and in the DUMAND ascii format used by cw                                */
/****************************************************************************/

char *rdmc_cvsid = 
"$Header: /net/local/cvsroot/siegmund/rdmc/rdmc.c,v 1.59 2004/04/15 14:45:52 wiebusch Exp $";

#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/////////////////#include <unistd.h>

#include <fcntl.h>

#if defined(OSF1) || defined(AIX) || defined (SunOS)
#include <errno.h>
#endif

#include "rdmc.h"

#include "rdmc_local.h"

#ifdef BAIKAL_BIN_F
#include "baikal.h"
#endif
#ifdef AMANDA_ASCII_F
#include "amanda.h"
#endif
#ifdef DUMAND_ASCII_F
#include "dumand.h"
#endif
#ifdef UWI_ASCII_F
#include "uwi.h"
#endif

#ifndef DEBUG
#define DEBUG 0                       /* DEBUG 1 writes some info to stderr */
#endif

#define PID180 (M_PI/180.)

/****************************************************************************/
/* write_parameters() simply adds the program name and the parameters to    */
/* the fp->creator field                                                    */
/****************************************************************************/

void rdmc_write_parameters(mcfile *fp, int argc,  char **argv, const char *version)
{
  int i;
  char *progname;
  int linelength = 0;

  progname = strrchr(argv[0],'/');
  if (progname == NULL) 
    progname = argv[0];
  else
    progname++;

  rdmc_append_comment(&(fp->creator)," ");
  rdmc_append_comment(&(fp->creator),progname);
  
  if (version != NULL) {
    rdmc_append_comment(&(fp->creator)," (");
    rdmc_append_comment(&(fp->creator),version);
    rdmc_append_comment(&(fp->creator),") ");
  }

  linelength = strlen(progname) + strlen(version) + 5;

  for (i = 1; i < argc; i++) {
    linelength += strlen(argv[i]) + 1;
    if (linelength > F2000_MAXLINE - RDMC_MAXTOKENLENGTH) {
      rdmc_append_comment(&(fp->creator), "\\\n");
      linelength = 1 + strlen(argv[i]);
    }
    rdmc_append_comment(&(fp->creator), " ");
    rdmc_append_comment(&(fp->creator),argv[i]);
  }

  rdmc_append_comment(&(fp->creator),"\n");

} /* write_parameters() */

/**********************************************/
/* appends a comment to a  character field */
/**********************************************/

void rdmc_append_comment(char **comment, const char *s){
  if (comment != NULL){ /* it is there at all */
    if (*comment == NULL) {
      *comment = (char*)malloc(strlen(s) + 1);
      strcpy(*comment, s);
    } else {
      *comment = (char*)realloc(*comment, strlen(*comment) + strlen(s) + 1);
      strcat(*comment, s);
    }
  }
}

void rdmc_concat_comment(char **comment, const char *s, const int form){
  switch(form){
#if 0     /* not needed */
  case   DUMAND_ASCII_F  :
        /* not needed */
    break;
#endif
  default:
    rdmc_append_comment(comment,"!");
    break;
  }
  rdmc_append_comment(comment,s);
}
