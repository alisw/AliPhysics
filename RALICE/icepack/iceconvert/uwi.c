/*
 * Read/write files in the format used at UWI Madison (see J.Jacobsen's PhD 
 * thesis)
 */

char *uwi_rdmc_cvsid = 
"$Header: /net/local/cvsroot/siegmund/rdmc/uwi.c,v 1.33 1999/12/23 12:23:15 wiebusch Exp $";

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

///////////////#include <unistd.h>


#include <fcntl.h>
#include <sys/stat.h>

#if defined(OSF1) || defined(AIX) || defined (SunOS)
#include <errno.h>
#endif

#include "rdmc.h"

#include "rdmc_local.h"

#include "uwi.h"

#ifndef DEBUG
#define DEBUG 0                       /* DEBUG 1 writes some info to stderr */
#endif

#define ERRLINE (-__LINE__)

#define PID180 (M_PI/180.)

#ifdef UWI_ASCII_F

int uwi_rd_EH(char *s, mevt *ev);   /* read the several lines for UWI files */
int uwi_rd_DH(char *s, mevt *ev);
int uwi_rd_MU(char *s, mtrack **tr, int *ntrack);
int uwi_rd_SH(char *s, mtrack **tr, int *ntrack);
int uwi_rd_HT(char *s, mhit **h, int *nhits);
int uwi_rd_FT(char *s, mtrack **fitr, mevt_special_t **fitres, int *nfit);
int uwi_rd_SP(char *s, mtrack **fitr, mevt_special_t **fitres, int *nfit);



/****************************************************************************/
/* The function rarr_uwi() reads the geometry of a UWI like file            */
/****************************************************************************/
int rdmc_rarr_uwi(char *geofile, array *ar)
{
  FILE *fp;
  char s[RDMC_MAXLINE];
  char omfile[RDMC_MAXLINE];
  int form;
  int ich = 0;
  int ichstr = 0;
  int istr;
  float x, y, z;
  int nchstr;          /* temp. x, y, z coor, number of channels per string */

  rdmc_init_array(ar);                               /* reset the array */

  if (geofile == NULL)      /* use the default file (environement variable) */
    geofile = getenv(DEFAULT_UWI_GEO_FILE_ENV);
  if (geofile == NULL)                              /* use the default file */
    geofile = DEFAULT_UWI_GEO_FILE;                            /* (builtin) */
  
  fp = fopen(geofile,"r");
  
  if (fp == NULL) return ERRLINE; /* no geometry file found */

  do {
    if (fgets(s,RDMC_MAXLINE,fp) == NULL)
      return EOF;
    form = sscanf(s,"%f",&(ar->depth));          /* 1. line: detector depth */
  } while (form < 1);
  do {
    if (fgets(s,RDMC_MAXLINE,fp) == NULL)
      return EOF;
    form = sscanf(s,"%i",&(ar->nstr));        /* 2. line: Number of strings */
  } while (form < 1);
  for (istr = 0; istr < ar->nstr; istr++) {             /* read all strings */
    do {
      if (fgets(s,RDMC_MAXLINE,fp) == NULL)
	return EOF;
      form = sscanf(s,"%f %f",&x,&y);         /* 3. line: x,y of 1st string */
    } while (form < 2);
    do {
      if (fgets(s,RDMC_MAXLINE,fp) == NULL)
	return EOF;
      form = sscanf(s,"%i",&nchstr);          /* 4. line: nch of 1st string */
    } while (form < 1);
    ar->nch += nchstr;        /* add the number of channels to the total nr */
    for (ichstr = 0; ichstr < nchstr; ich++, ichstr++) {  /* read ch in str */
      do {
	if (fgets(s,RDMC_MAXLINE,fp) == NULL)
	  return EOF;
	form = sscanf(s,"%f",&z);              /* 5. line: z of 1st channel */
      } while (form < 1);
      do {
	if (fgets(s,RDMC_MAXLINE,fp) == NULL)
	  return EOF;
	form = sscanf(s,"%s",omfile);         /* 6. line: om file of 1st ch */
      } while (form < 1);
      ar->x[ich] = x;
      ar->y[ich] = y;
      ar->z[ich] = z;
      ar->str[ich] = istr+1;
      ar->type[ich] = STD_AMANDA_B_OM;
      ar->clust[ich] = ichstr;
      if (strcmp(omfile,"up.pmt") == 0){
	ar->costh[ich] = 1.0;
        ar->sensit[ich] = 1.0;
      }else if (strcmp(omfile,"down.pmt") == 0){
	ar->costh[ich] = -1.0;
	ar->sensit[ich] = 1.0;
      }else if (strcmp(omfile,"dead.pmt") == 0){
	ar->costh[ich] = 0.0;
	ar->sensit[ich] = 0.0;
      } else{
	ar->costh[ich] = 0.0;
	ar->sensit[ich] = 0.0;
      }
  } /* for ich */
  } /* for istr */
  
  ar->id = AMANDA_B_4;                             /* default: amanda array */
  ar->longitude = -90.0;                  /* geographic longitude of AMANDA */
  ar->lattitude = 0.0;                              /* geographic lattitude */

  {
    array_hdef_t def_trig;

    rdmc_init_array_hdef(&def_trig);
    def_trig.id = 1;
    strcpy(def_trig.tag, "amanda-a");
    def_trig.npars = 1; 
    strcpy(def_trig.pars[0], "type=external");
    def_trig.nwords = 0;
    rdmc_add_trigger_def(ar,&def_trig,ar->n_trigger);

    rdmc_init_array_hdef(&def_trig);
    def_trig.id = 2;
    strcpy(def_trig.tag, "amanda-b");
    def_trig.npars = 1; 
    strcpy(def_trig.pars[0], "type=majority");
    def_trig.nwords = 0;
    rdmc_add_trigger_def(ar,&def_trig,ar->n_trigger);

    rdmc_init_array_hdef(&def_trig);
    def_trig.id  = 3;
    strcpy(def_trig.tag, "spase-1");
    def_trig.npars = 1; 
    strcpy(def_trig.pars[0], "type=external");
    def_trig.nwords = 0;
    rdmc_add_trigger_def(ar,&def_trig,ar->n_trigger);

    rdmc_init_array_hdef(&def_trig);
    def_trig.id  = 4;
    strcpy(def_trig.tag, "spase-2");
    def_trig.npars = 1; 
    strcpy(def_trig.pars[0], "type=external");
    def_trig.nwords = 0;
    rdmc_add_trigger_def(ar,&def_trig,ar->n_trigger);
    
    rdmc_init_array_hdef(&def_trig);
    def_trig.id = 5;
    strcpy(def_trig.tag, "gasp");
    def_trig.npars = 1; 
    strcpy(def_trig.pars[0], "type=external");
    def_trig.nwords = 0;
    rdmc_add_trigger_def(ar,&def_trig,ar->n_trigger);
  }

  /* now set the geo cal flag */
  ar->is_calib.geo=1;

  /* now add one FRESULT line */
  { 
    array_hdef_t fd;
    rdmc_init_array_hdef(&fd);
    rdmc_add_fit_def(ar,&fd,0);
    rdmc_jk_to_fitdef(ar->def_fit, 1);
  }
  return 0;

} /* rarr_uwi() */

/****************************************************************************/
/* revt_uwi() reads the next UWI format event                               */
/****************************************************************************/

int rdmc_revt_uwi(mcfile *fp, mevt *ev, const array *ar)
{

  char s[RDMC_MAXLINE];                                            /* input line */
  int c;                                                      /* first char */
  int header_read = 0;               /* detect if the event header was read */

  if (fp->fmajor < 2)  return ERRLINE; /* only newer formats are supported */

  rdmc_clear_mevt(ev);                           /* clear old event info */

  if (feof(fp->fp))                          /* if file end already reached */
    return EOF;

  ev->nrun = ar->nrun;                       /* copy the run number from fp */

  for (*s = c = getc(fp->fp); c != EOF; *s=c=getc(fp->fp)) {/*for all lines */
    fp->fpos += 1;
    if (fgets(s+1,RDMC_MAXLINE-1,fp->fp) == NULL) {        /* try to read a line */
      return EOF;
    }

    if (strncmp(s,"DH",2) == 0) {                      /* data event header */
      if (header_read++) return ERRLINE;
      if (uwi_rd_DH(s+2, ev)) return ERRLINE;
    } else if (strncmp(s,"EH",2) == 0) {                /* MC event header */
      if (header_read++) return ERRLINE;
      if (uwi_rd_EH(s+2, ev)) return ERRLINE;
    } else if (strncmp(s,"MU",2) == 0) {          /* Muon track information */
      if (uwi_rd_MU(s+2, &(ev->gen), &(ev->ntrack))) return ERRLINE;
    } else if (strncmp(s,"SH",2) == 0) {             /* shower information */
      if (uwi_rd_SH(s+2, &(ev->gen), &(ev->ntrack))) return ERRLINE;
    } else if (strncmp(s,"FT",2) == 0) {                    /* Fit results */
      if (uwi_rd_FT(s+2, &(ev->rec), &(ev->fresult), &(ev->nfit))) 
	return ERRLINE;
    } else if (strncmp(s,"SP",2) == 0) {               /* Spase Fit results */
      if (uwi_rd_SP(s+2, &(ev->rec), &(ev->fresult), &(ev->nfit))) 
	return ERRLINE;
    } else if (strncmp(s,"HT",2) == 0) {                /* Hit information */
      if (uwi_rd_HT(s+2, &(ev->h), &(ev->nhits))) return ERRLINE;
    } else if (strncmp(s,"EN",2) == 0) {                      /* Event end */

   /* now fill necessary data structures */
      rdmc_repair_gen_id(ev); 
      ev->nch = rdmc_count_nch(ev);     /* calc the number of hit channels */
      rdmc_fill_mhit_str(ev,ar);
      ev->nstr = rdmc_count_nstr(ev);
      


      return 0;
    } else                                                   /* unknown id */
      continue;  /* do nothing */

  } /* for all lines */

  if (c == EOF)
    return EOF;

  ev->nch=rdmc_count_nch(ev);                         /* calc the number of hit channels */

  return 0;

} /* revt_uwi() */


/****************************************************************************/
/* skipevt_uwi() skips the next event from a UWI-like file                  */
/****************************************************************************/

int rdmc_skipevt_uwi(mcfile *fp)
{
  char s[RDMC_MAXLINE];                                            /* input line */
  int c;                                                      /* first char */

  fp->fpos += 1;
  if (fgets(s+1,RDMC_MAXLINE-1,fp->fp) == NULL)        /* try to read first line */
    return EOF;

  for (*s = c = getc(fp->fp); c != EOF; *s=c=getc(fp->fp)) {/*for all lines */
    fp->fpos += 1;
    if (fgets(s+1,RDMC_MAXLINE-1,fp->fp) == NULL) {        /* try to read a line */
      return EOF;
    }

    if (strncmp(s,"EN",2) == 0)                                /* Event end */
      return 0;
    
  } /* for all lines */


  return 0;

} /* function skipevt_ascii() */


/****************************************************************************/
/* warr_uwi() may write the geometry to a UWI like geometry file            */
/*  This function is somehow funny: It first looks if there is already a    */
/*  geometry file. If yes, it compares them, and returns a zero if the      */
/*  geometries are the same, and and error if not.                          */
/*  only if there is no geometry file, we will write it.                    */
/****************************************************************************/

int rdmc_warr_uwi(char *geofile,const array *geo)
{
  FILE *fp;
  array stored_array;
  int r;
  int istr, ich, nchstr;
  float xstr, ystr;

  if (geofile == NULL)      /* use the default file (environement variable) */
    geofile = getenv(DEFAULT_UWI_GEO_FILE_ENV);
  if (geofile == NULL)                              /* use the default file */
    geofile = DEFAULT_UWI_GEO_FILE;                            /* (builtin) */
  
  fp = fopen(geofile,"r");                  /* try to open the file to read */
  if (fp != NULL) {                    /* there was already a geometry file */
    fclose(fp);                                      /* just close it again */
    if ((r = rdmc_rarr_uwi(geofile,&stored_array)))    /* read the geometry */
      return r;                                             /* if errornous */
    if (rdmc_comp_array(geo, &stored_array) == 1)        /* compare arrays */
      return ERRLINE;
    else 
      return 0;
  } else {                                 /* if there was no geometry file */
    fp = fopen(geofile,"w");          /* open the geometry file for writing */
    if (fp == NULL) return ERRLINE;
    fprintf(fp, "%.0f\t\t! Depth of detector center\n", geo->depth);
    fprintf(fp, "%i\t\t! Number of strings\n", geo->nstr);
    for (istr = 1; istr <= geo->nstr; istr++) {          /* for all strings */
      nchstr = 0;
      xstr = 0.0; ystr = 0.0;
      for (ich = 0; ich < geo->nch; ich++)                   /* calc nchstr */
	if (geo->str[ich] == istr) {
	  nchstr++;
	  xstr += geo->x[ich];
	  ystr += geo->y[ich];
	} /* if str[ich] == istr */
      if (nchstr > 0) {
	xstr /= nchstr;                             /* mean x, y coordinate */
	ystr /= nchstr;
	fprintf(fp, "\n%.3f %.3f\t! x, y of string %i\n",xstr, ystr, istr);
	fprintf(fp, "%i\t\t! Nr of OMs\n\n",nchstr);
	nchstr = 0;
	for (ich = 0; ich < geo->nch; ich++)
	  if (geo->str[ich] == istr) {
	    nchstr++;
	    fprintf(fp, "%.3f\t\t! z of OM %i\n",geo->z[ich],nchstr);
	    if (geo->costh[ich] < -0.5) fprintf(fp, "down.pmt\n");
	    else if (geo->costh[ich] > 0.5) fprintf(fp, "up.pmt\n");
	    else fprintf(fp, "unknown.pmt\n");
	} /* if str[ich] == istr */
      } /* if nchstr */
     
    } /* for istr */
    fclose(fp);
    return 0;
  } /* if fp == NULL */

} /* warr_uwi() */

/****************************************************************************/
/* whd_uwi() writes the header to UWI file                                  */
/****************************************************************************/

int rdmc_whd_uwi(const mcfile *fp)
{
#if 0
  int day;
  day = gmtime(&(fp->time))->tm_yday;
#endif

  if (fp->info.uwi.mc_id == 'D')
    fprintf(fp->fp, "FH %i DATA unknown\n", UWI_ASCII_VERSION);
  else if (fp->info.uwi.mc_id == 'M')
    fprintf(fp->fp, "FH %i MC %s\n",UWI_ASCII_VERSION, fp->info.uwi.mc_vers);
  else
    fprintf(fp->fp, "FH %i MC unknown\n", UWI_ASCII_VERSION);

  return 0;

} /* whd_uwi() */

/****************************************************************************/
/* function wevt_uwi() writes an event to an UWI file                       */
/****************************************************************************/

int rdmc_wevt_uwi(const mcfile *fp,const mevt *event, const array *ar)
{
  int nmuon = 0;
  int nshower = 0;
  int imuon = -1;
  int shwr_type;
  int itrack;
  int ihit;
  int gpsyear, gpsday;

  /* look how many muons we have */
  for (itrack = 0; itrack < event->ntrack; itrack++)
    if ((event->gen[itrack].id == MUON_PLUS) 
	|| (event->gen[itrack].id == MUON_MINUS)) {
      nmuon++;
      if (imuon < 0) imuon = itrack; /* pointer to first muon */
    } else
      nshower++;

  /* write event header  */
  if (fp->info.uwi.mc_id == 'D') {
    rdmc_mjd_to_gps(event->mjd, &gpsyear, &gpsday);
    fprintf(fp->fp, "DH %i %i %i.%06li %i.%06li %i %li %i\n",
	    gpsyear, gpsday, event->secs,
	    rdmc_nint(event->nsecs*1e-3),
	    event->secs + 86400 * gpsday, 
	    rdmc_nint(event->nsecs*1e-3),
	    event->enr, event->trigger, event->nhits);
  } else {
    if (nmuon > 0)
      fprintf(fp->fp, "EH %i %i %i %.0f %.2f %.2f %.0f\n",
	      event->enr, nmuon, event->nch, 
	      -1.0, /* primary energy */
	      acos(event->gen[imuon].costh)/PID180,
	      event->gen[imuon].phi/PID180,
	      -1.0 /* AMU mass of primary */
	      );
    else
      fprintf(fp->fp, "EH %i %i %i\n", event->enr, 0, event->nch);
  }

  /* write (gen) muons */
  nmuon = 0;
  for (itrack = 0; itrack < event->ntrack; itrack++) {
    float len;
    shwr_type = -1;
    switch (event->gen[itrack].id) {
    case MUON_PLUS:
    case MUON_MINUS:
      nmuon++;
      for (nshower = 1; itrack+nshower < event->ntrack; nshower++)
	if ((event->gen[itrack+nshower].id == MUON_PLUS)
	    || (event->gen[itrack+nshower].id == MUON_MINUS))
	  break;

      fprintf(fp->fp, "MU %i %.0f %i %.2f %.2f ", 
	      nmuon, event->gen[itrack].e * 1e-3, 
	      nshower,
	      acos(event->gen[itrack].costh)/PID180,
	      event->gen[itrack].phi/PID180);
      len = ( event->gen[itrack].length < 0. ) 
	? RDMC_BIG : event->gen[itrack].length ;
      fprintf(fp->fp, "%.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
	      event->gen[itrack].x * 1e2,
	      event->gen[itrack].y * 1e2,
	      event->gen[itrack].z * 1e2,
	      (event->gen[itrack].x + 
	       len  * event->gen[itrack].px) * 1e2,
	      (event->gen[itrack].y + 
	       len * event->gen[itrack].py) * 1e2,
	      (event->gen[itrack].z + 
	       len * event->gen[itrack].pz) * 1e2,
	      0.0 /* closest distance to origin */
	      );
      break;
    } /* switch id */
  } /* for itrack */  

  /* write (gen) showers */
  nmuon = 0;
  nshower = 0;
  for (itrack = 0; itrack < event->ntrack; itrack++) {
    shwr_type = -1;
    switch (event->gen[itrack].id) {
    case MUON_PLUS:
    case MUON_MINUS:
      nmuon++;
      nshower = 0;
      break;
    case BREMS: if (shwr_type < 0) shwr_type = 1;
    case PAIRPROD: if (shwr_type < 0) shwr_type = 2;
    case NUCL_INT: if (shwr_type < 0) shwr_type = 3;
    case DELTAE: if (shwr_type < 0) shwr_type = 4;
    case MU_PAIR: if (shwr_type < 0) shwr_type = 5;
    default:                                 /* store unknown particles */
      nshower++;
      fprintf(fp->fp, "SH %i %i %i %.0f %.0f %.1f %.1f %.1f\n",
	      nmuon, /* muon responsible for this shower */
	      nshower+1, shwr_type,
	      event->gen[itrack].e * 1e-3, 
	      0.0, /* light output from shower */
	      event->gen[itrack].x * 1e2,
	      event->gen[itrack].y * 1e2,
	      event->gen[itrack].z * 1e2);	      
    } /* switch id */
  } /* for itrack */  

  /* write hits */
  for (ihit = 0; ihit < event->nhits; ihit++)
    fprintf(fp->fp, "HT %i %i %i %.1f %.1f %.1f\n",
	    event->h[ihit].ch+1,
	    0, /* nr of p.e. generated */
	    1, /* nr. of output pulses generated */
	    event->h[ihit].amp * UWI_CHANNELS_PER_PE,
	    event->h[ihit].t,
	    event->h[ihit].tot);
	    
  /* write fits */
  for (itrack = 0; itrack < event->nfit; itrack++){
    fprintf(fp->fp, "FT %.1f %.1f %.1f %.1f %.1f",
	    acos(event->rec[itrack].costh)/PID180,
	    event->rec[itrack].phi/PID180,
	    event->rec[itrack].x * 1e2,
	    event->rec[itrack].y * 1e2,
	    event->rec[itrack].z * 1e2);
 
    if ((event->fresult != NULL)              /* if there is an jk record */
	&& (0 <= event->fresult[itrack].id )
	&& (ar->n_fit > event->fresult[itrack].id )  /* there is a fit defined  */
	&& (rdmc_is_this_jk(&(ar->def_fit[event->fresult[itrack].id])
			       ,&(event->fresult[itrack])) ) ) {
      fprintf(fp->fp, " %f", event->fresult[itrack].val[JK_CHI2]);
    }else{
      fprintf(fp->fp, " %f",(float) RDMC_NA);
    }
    fprintf(fp->fp, " %f %i\n" , event->rec[itrack].t,  0);
    
  }
  /* write end of event */
  fprintf(fp->fp,"EN\n");

  return 0;
} /* wevt_uwi() */

/****************************************************************************/
/* rhd_uwi() reads the format relevant informations for UWI like            */
/*          formats                                                         */
/****************************************************************************/

int rdmc_rhd_uwi(mcfile *fp) 
{
  return 0;
} 

int rdmc_uwi_FH(const char *s, mcfile *fp)
{
  int form;
  int format, mcversion;
  char filename[RDMC_MAXLINE];
  char histline[RDMC_MAXLINE];

  form = sscanf(s,"FH %i MC %i", &format, &mcversion);
  if (form >= 2) {                                         /* MC file */
    fp->fmajor = format;
    fp->info.uwi.mc_id = 'M';
    sprintf(fp->info.uwi.mc_vers,"%i",mcversion);
    return 0;
  }
  form = sscanf(s, "FH %i DATA %s", &format, filename);
  if (form >= 1) {                                        /* data file */
    fp->fmajor = format;
    fp->fminor = 0;
    fp->info.uwi.mc_id = 'D';
    strcpy(fp->info.uwi.mc_vers,"0.0");

    sprintf(histline,"uwicalib (%i) %s\n",format, (form > 1)?filename:"");
    if (fp->creator == NULL) {
      fp->creator = malloc(strlen(histline)+2);
      strcpy(fp->creator, "!");
      strcat(fp->creator, histline);
    } else {
      fp->creator = realloc(fp->creator, 
			    strlen(fp->creator) + strlen(histline)+1);
      strcat(fp->creator, "!");
      strcat(fp->creator, histline);
    }
    return 0;
  } 

  return ERRLINE;

} /* rhd_uwi() */

/****************************************************************************/
/* Functions for reading UWI format files                                   */
/****************************************************************************/

/* Read the "EH " line (just ignore the primary) */

int uwi_rd_EH(char *s, mevt *ev)
{
  int form;
  float prim_energy;
  float prim_theta, prim_phi;
  float prim_mass;
  int nparts;
  int nhits;

  form = sscanf(s,"%i %i %i %f %f %f %f",
                &(ev->enr), &(nparts), &nhits,
		&prim_energy, &prim_theta, &prim_phi, &prim_mass);


  if (form == 0) 
    ev->enr = 0;

  return 0;

} /* uwi_rd_EH() */

/* read a Data event header */

int uwi_rd_DH(char *s, mevt *ev)
{
  int form;
  int reg;                  /* Data register */
  int nhits;
  int year, day, sec, nsec; /* time */
  float gmt, fnsecs;
  int enr;
  char nsecstr[RDMC_MAXLINE];
  char nsecfloat[RDMC_MAXLINE];

  form = sscanf(s, "%i %i %i.%s %f %i %i %i",
		&year, &day, &sec, nsecstr, &gmt, &enr, &reg, &nhits);

  switch(form) {
  case 0: year = 70;
  case 1: day = 0;
  case 2: sec = 0;
  case 3: strcpy(nsecstr, "0");
  case 4: gmt = 0.0;
  case 5: enr = 0;
  case 6: reg = 0;
  case 7: nhits = 0;
  }

  sprintf(nsecfloat, "0.%s",nsecstr);
  sscanf(nsecfloat, "%f", &fnsecs);
  nsec = 1e6 * fnsecs;
  nsec *= 1000;

  ev->mjd = rdmc_gps_to_mjd(year, day);  /* convert gps date into mjd */
  ev->secs = sec;
  ev->nsecs = nsec;
  ev->trigger = reg;   /* I hope! that reg contains the trigger bitfield */
                       /* (but where can I get the trigger info?)        */
  ev->enr = enr;

  return 0;

} /* uwi_rd_DH() */

/* read a muon track and reserve the memory for all shower tracks */

int uwi_rd_MU(char *s, mtrack **tr, int *ntrack)
{
  int form;
  float l;
  int muon_number, num_showers;
  float muon_energy;
  float theta_mu, phi_mu;
  float xstart, ystart, zstart;
  float xend, yend, zend;
  mtrack muon;
  int i;

  form = sscanf(s, "%i %f %i %f %f %f %f %f %f %f %f",
		&muon_number, &muon_energy, &num_showers,
		&theta_mu, &phi_mu, 
		&xstart, &ystart, &zstart,
		&xend, &yend, &zend);

  switch (form) {
  case 8: 
    xend = xstart;
    yend = ystart;
    zend = zstart;
  case 11:
  case 12: break;
  default: return ERRLINE;
  }

  if (num_showers < 1) return ERRLINE;

  *tr = (mtrack *)realloc(*tr, sizeof(mtrack) * (*ntrack+num_showers));

  muon.id = MUON_PLUS;
  muon.e = muon_energy * 1e3; /* file: GeV, rdmc: MeV */
  muon.x = xstart*1e-2;
  muon.y = ystart*1e-2;
  muon.z = zstart*1e-2;
  muon.phi = phi_mu * PID180;
  muon.costh = - cos(theta_mu*PID180);
  l = (xstart - xend) * (xstart - xend)
    + (ystart - yend) * (ystart - yend)
    + (zstart - zend) * (zstart - zend);
  if (l < 0.) 
    l = RDMC_NA;
  else
    l = sqrt(l);
  muon.length = l * 1e-2; /* file: cm, rdmc: m */

  if (l < 1e-12) /* if the track is too short for a good calculation */
    rdmc_tau_tr
(*tr);     /* calc direction cosinuus from phi and costh */
  else {
    muon.px = (xend - xstart) / l;/* else calc them from end-start */
    muon.py = (yend - ystart) / l;
    muon.pz = (zend - zstart) / l;
  }

  muon.t = 0.;
  muon.nmuon = 1;
  
  for (i = 0; i < num_showers; i++)            /* copy the muon to the event */
    memcpy((*tr)+*ntrack+i, &muon, sizeof(mtrack));   /* (dummy for showers) */

  for (i = 1; i < num_showers; i++)
    (*tr)[*ntrack+i].id = -1;                         /* dummies for showers */

  (*ntrack) += num_showers;
  
  return 0;

} /* uwi_rd_MU() */

/* read a shower and fill it to the (in uwi_rd_MU allocated) track */

int uwi_rd_SH(char *s, mtrack **tr, int *ntrack)
{
  int form;
  int which_muon, which_shower, showr_type;
  float shwr_energy, shwr_photons;
  float x, y, z;
  int i;
  mtrack shower;

  form = sscanf(s, "%i %i %i %f %f %f %f %f",
		&which_muon, &which_shower,
		&showr_type, 
		&shwr_energy, &shwr_photons,
		&x, &y, &z);

  if (form < 8) return ERRLINE;

  switch (showr_type) {
  case 1: shower.id = BREMS; break;
  case 2: shower.id = PAIRPROD; break;
  case 3: shower.id = NUCL_INT; break;
  case 4: shower.id = DELTAE; break;
  case 5: shower.id = MU_PAIR; break;
  default: shower.id = -1;
  }

  shower.e = shwr_energy * 1e3;
  shower.x = x * 1e-2;
  shower.y = y * 1e-2;
  shower.z = z * 1e-2;

  /* lets look for the muon number to fill into the track direction */

  for (i = 0; (which_muon > 0) && (i < *ntrack); i++)
    if (((*tr)[i].id == MUON_PLUS) || ((*tr)[i].id == MUON_MINUS))
      which_muon--;

  if (i >= *ntrack) return ERRLINE;

  shower.px = (*tr)[i].px;
  shower.py = (*tr)[i].py;
  shower.pz = (*tr)[i].pz;
  shower.phi = (*tr)[i].phi;
  shower.costh = (*tr)[i].costh;
  
  shower.length = 0; /* a shower has zero length (lets assume)  */
  shower.nmuon = 0;

  shower.t = 0; /* normally, we should calculate the time which 
			 * corresponds to the  track I do not do this.
			 */

  /* insert the shower at the right place (after the corresponding muon) */

  memcpy((*tr) + i + which_shower - 2, &shower, sizeof(mtrack));

  return 0;

} /* uwi_rd_SH() */

int uwi_rd_HT(char *s, mhit **h, int *nhits)
{
  int form;
  int om_num, n_pulses;
  float adc, n_pe;
  float tdc0, tdc1, tdc2, tdc3, tdc4;
  float tot0, tot1, tot2, tot3, tot4;

  form = sscanf(s, "%i %f %i %f %f %f %f %f %f %f %f %f %f %f",
		&om_num, &n_pe, &n_pulses,
		&adc,
		&tdc0, &tot0, 
		&tdc1, &tot1, 
		&tdc2, &tot2, 
		&tdc3, &tot3, 
		&tdc4, &tot4);

  /*
   * We will skip all hits without hits!!! (That means, without TDC)
   */

  if (form < 4)
    return ERRLINE;
  if (n_pulses < 0) return ERRLINE;

  if ((n_pulses == 0) && (form == 6))
    n_pulses = 1;

  if (n_pulses > 5) n_pulses = 5;
  if (form < (4 + 2*n_pulses) )
    n_pulses = (form - 4 )/2;
    
  *h = (mhit *)realloc(*h, sizeof(mhit) * (*nhits + n_pulses));
    
  switch(n_pulses) {
  case 5:
    (*h)[*nhits+4].ch = om_num-1;
    (*h)[*nhits+4].str = 0;  /* default: no string (will be set at the end) */
    (*h)[*nhits+4].mt = RDMC_NA;
    (*h)[*nhits+4].ma = RDMC_NA;
    (*h)[*nhits+4].t = tdc4; /* units are nanonseconds */
    (*h)[*nhits+4].tot = tot4;
    (*h)[*nhits+4].amp = 0.;
  case 4:
    (*h)[(*nhits)+3].ch = om_num-1;
    (*h)[(*nhits)+3].str = 0; 
    (*h)[(*nhits)+3].mt = RDMC_NA;
    (*h)[(*nhits)+3].ma = RDMC_NA;
    (*h)[(*nhits)+3].t = tdc3;
    (*h)[(*nhits)+3].tot = tot3;
    (*h)[(*nhits)+3].amp = 0.;
  case 3:
    (*h)[(*nhits)+2].ch = om_num-1;
    (*h)[(*nhits)+2].str = 0; 
    (*h)[(*nhits)+2].mt = RDMC_NA;
    (*h)[(*nhits)+2].ma = RDMC_NA;
    (*h)[(*nhits)+2].t = tdc2;
    (*h)[(*nhits)+2].tot = tot2;
    (*h)[(*nhits)+2].amp = 0.;
  case 2:
    (*h)[(*nhits)+1].ch = om_num-1;
    (*h)[(*nhits)+1].str = 0; 
    (*h)[(*nhits)+1].mt = RDMC_NA;
    (*h)[(*nhits)+1].ma = RDMC_NA;
    (*h)[(*nhits)+1].t = tdc1;
    (*h)[(*nhits)+1].tot = tot1;
    (*h)[(*nhits)+1].amp = 0.;
  case 1:
    (*h)[*nhits].ch = om_num-1;
    (*h)[*nhits].str = 0; 
    (*h)[*nhits].mt = RDMC_NA;
    (*h)[*nhits].ma = RDMC_NA;
    (*h)[*nhits].t = tdc0;
    (*h)[*nhits].tot = tot0;
    (*h)[*nhits].amp = adc / UWI_CHANNELS_PER_PE;
    break;
  case 0: 
    break;
  default: return ERRLINE;
  }

  (*nhits) += n_pulses;

  return 0;
} /* uwi_rd_HT() */

int uwi_rd_SP(char *s, mtrack **fit, mevt_special_t **fitres, int *nfit)
{ /* Tmiller: SP gps_time theta phi x y Nparticles dummy(currently 0) */
  /* serap:   SP gmt, theta, phi, x, y, s30, xxx */

  int form;
  float theta_fit, phi_fit;
  float x_fit, y_fit;
  float s30,xxx;
  double gmt;
  
  form = sscanf(s, "%lf %f %f %f %f %f %f",
		&gmt,&theta_fit, &phi_fit, 
		&x_fit, &y_fit, &s30,
		&xxx);

  switch (form) {
  case 5: s30 = 0;
  case 6: xxx = 0.;
  case 7: break;
  default: return ERRLINE;
  }

  *fit = (mtrack *)realloc(*fit, sizeof(mtrack) * (*nfit+1));
  *fitres = (mevt_special_t *) 
    realloc( *fitres , sizeof(mevt_special_t)  *  (*nfit+1));
  rdmc_init_mtrack(&((*fit)[*nfit]));
  rdmc_init_fit_jk(&((*fitres)[*nfit]),RDMC_NA);
  
  
  (*fit)[*nfit].id = MUON_PLUS;
  (*fit)[*nfit].e = 0.0;
  (*fit)[*nfit].x = x_fit;
  (*fit)[*nfit].y = y_fit;
  (*fit)[*nfit].z = 0.0;
  (*fit)[*nfit].phi = phi_fit * PID180;
  (*fit)[*nfit].costh = cos(theta_fit*PID180);
  (*fit)[*nfit].length = RDMC_BIG; 
  (*fit)[*nfit].tag = *nfit; 

  rdmc_tau_tr(*fit + *nfit); /* calc direction cosinus from phi and costh */

  (*fit)[*nfit].nmuon = 1;  /* we dont know this, lets assume one muon */

  (*fitres)[*nfit].id = 0;
  (*fitres)[*nfit].val[JK_FITID] = 9000;
  (*fitres)[*nfit].val[JK_CHI2] = 0.;
  (*fitres)[*nfit].val[JK_RCHI2] = 0.;
  (*fitres)[*nfit].val[JK_PROB] = 0.;
  (*fitres)[*nfit].val[JK_SIGTH] = 0.;
  (*fitres)[*nfit].val[JK_COVMIN] = s30;
  (*fitres)[*nfit].val[JK_COVMAX] = xxx;
  (*fitres)[*nfit].val[JK_CUTFLAG] = -1;
  
  (*nfit)++;
  
  return 0;
  
} /* uwi_rd_FT() */

int uwi_rd_FT(char *s, mtrack **fit, mevt_special_t  **fitres, int *nfit)
{
  int form;
  float theta_fit, phi_fit;
  float x_fit, y_fit, z_fit;
  float chi2;
  int n_in_time;
  float t_fit;

  form = sscanf(s, "%f %f %f %f %f %f %f %i",
		&theta_fit, &phi_fit, 
		&x_fit, &y_fit, &z_fit,
		&chi2, &t_fit, &n_in_time);

  switch (form) {
  case 5: chi2 = 0.0;
  case 6: t_fit = 0.0;
  case 7: n_in_time = 0;
  case 8: break;
  default: return ERRLINE;
  }

  *fit = (mtrack *)realloc(*fit, sizeof(mtrack) * (*nfit+1));
  *fitres = (mevt_special_t *) 
    realloc( *fitres , sizeof(mevt_special_t)  *  (*nfit+1));
  rdmc_init_mtrack(&((*fit)[*nfit]));
  rdmc_init_fit_jk(&((*fitres)[*nfit]),RDMC_NA);



  (*fit)[*nfit].id = MUON_PLUS;
  (*fit)[*nfit].e = 0.0;
  (*fit)[*nfit].x = x_fit*1e-2;
  (*fit)[*nfit].y = y_fit*1e-2;
  (*fit)[*nfit].z = z_fit*1e-2;
  (*fit)[*nfit].t = t_fit;
  (*fit)[*nfit].phi = phi_fit * PID180;
  (*fit)[*nfit].costh = cos(theta_fit*PID180);
  (*fit)[*nfit].length = RDMC_BIG; 
  (*fit)[*nfit].tag = *nfit; 

  rdmc_tau_tr(*fit + *nfit); /* calc direction cosinuus from phi and costh */

  (*fit)[*nfit].nmuon = 1;     /* we dont know this, lets assume one muon */

  (*fitres)[*nfit].id = 0;
  (*fitres)[*nfit].val[JK_FITID] = 0;
  (*fitres)[*nfit].val[JK_CHI2] = chi2;
  (*fitres)[*nfit].val[JK_RCHI2] = -1.0;
  (*fitres)[*nfit].val[JK_PROB] = 0.0;
  (*fitres)[*nfit].val[JK_SIGTH] = 0.0;
  (*fitres)[*nfit].val[JK_COVMIN] = 0.0;
  (*fitres)[*nfit].val[JK_COVMAX] = 0.0;
  (*fitres)[*nfit].val[JK_CUTFLAG] = -1;
  
  (*nfit)++;
  
  return 0;
  
} /* uwi_rd_FT() */

#endif /* UWI_ASCII_F */

/****************************************************************************/
/********************************** E O F ***********************************/
/****************************************************************************/
/* 
   This is just for EMACS:
   Local Variables:
   compile-command: "cd .. ; make -k rdmc" 
   End: 
*/
