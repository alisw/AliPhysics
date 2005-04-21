#ifndef _RDMC_F2K_H
#define _RDMC_F2K_H

#include "rdmc.h"

/******************************************************/
/* define the posible lines which contribute to f2k */
/******************************************************/
#define F2K_MAX_EVENTTYPES 10
#define F2K_MAX_LINETYPES 100

enum F2K_LINETYPES_T {
  V_LINE=1  
  , HI_LINE 
  , COMMENT_LINE
  , ARRAY_LINE
  , FIT_DEF_LINE
  , STAT_DEF_LINE
  , USER_DEF_LINE
  , TRIG_DEF_LINE
  , TRIG_PAR_LINE
  , KH_LINE
  , OM_LINE
  , KADC_LINE
  , KTOT_LINE
  , KTDC_LINE
  , KUTC_LINE
  , TBEGIN_LINE
  , EM_LINE
  , US_LINE
  , HT_LINE
  , WF_LINE
  , CH_LINE
  , TR_LINE
  , FIT_LINE
  , FRESULT_LINE
  , TRIG_LINE
  , STATUS_LINE
  , USES_LINE
  , EE_LINE
  , TEND_LINE
  , END_LINE
  , DUMMY_LINE
};

enum F2K_TAG_SEARCH_T {COMP_STRINGWISE=1, /* stringwise comp */
		       COMP_STRINGWISE_NOTRAIL,
		       COMP_CHARWISE, /* one char should fit */
		       COMP_CHARPUNCT /* same, char is isspunct -> speed */ 
#if 0
		       ,COMP_NOTHING /* do not compare at all ->always true */ 
#endif
};


/************** functions to actually parse the lines ****/


typedef int (f2k_parser_fun_t)(mcfile *fp, array *a, 
			       mevt *e, void *misc);
/* the void pointer misc is a trick to pass anny target to fill */

 
typedef struct  { 
  char tag[RDMC_MAXTOKENLENGTH]; 
  enum F2K_LINETYPES_T line_id;
  enum F2K_TAG_SEARCH_T searchtype;
  f2k_parser_fun_t *parser;
} f2000_line_t;

/******************************************************/
/*** Definitions of f2k blocks                     ***/
/******************************************************/

typedef int (f2k_event_reader_fun_t)(mcfile *fp, array *a, mevt *e);
typedef int (f2k_event_writer_fun_t)(const mcfile *fp, const array *a, const mevt *e);

typedef struct{
  enum RDMC_EVENT_T id;
   const  f2000_line_t *opener[F2K_MAX_LINETYPES];/* lines that are required to indicate such an event */
  const  f2000_line_t *inner[F2K_MAX_LINETYPES]; /* lines that re alowd within such n event */
  const  f2000_line_t *closer[F2K_MAX_LINETYPES];/* lines that define the end of such an event */
  f2k_event_reader_fun_t *reader;
  f2k_event_writer_fun_t *writer;
} f2000_event_t;


/******************************************************/
/**** buffer for decoding f2k ************************/
/******************************************************/
#define F2K_BUFFERSIZE 32768 /* initial size of the raw text buffer */
#define F2K_LINE_BUFFERSIZE 2048 /* initial size of the raw text buffer */



typedef struct rdmc_f2k_buffer_s { /* a structure to keep track of raw lines */
  unsigned long int used; /*length of the buffer (excluding last \0) */
  unsigned long int ntot;   /* number of allocated chars in buff */ 
  char *buff;             /* this buffers the raw text  */
  /* now the decoding of the lines */
  unsigned long int lines;   /* number of gathered lines in buff */ 
  unsigned long int lines_tot;   /* number of allocated lines in buff */ 
  unsigned long int iline;    /* the current line under investigation */
  const f2000_line_t **type_def;  /* pointer to the global array definition */
  char **line_pt;               /* pointer to each line in buffer */
} rdmc_f2k_buffer_t; 
/* the last element of line_pt is always set to NULL */


/* **************************************************** */
/* now this are all known event functions  */
/* *****************************************************/

// the different formats differe only in the events
extern const f2000_event_t f2000_preamble_1x1; /* works for 1.1, 1.2, 2004.1*/
extern const f2000_event_t f2000_mhead_1x1; /* works for 1.1, 1.2, 2004.1*/
extern const f2000_event_t f2000_mfoot_1x1;

extern const f2000_event_t *f2000_events_1x1[];
extern const f2000_event_t *f2000_events_1x2[];
extern const f2000_event_t *f2000_events_2004x1[];


/* ******************************************************************* */
/* Functions in f2k_utl.c */
/* ******************************************************************* */



/* functions that convert special strings to rdmc values and vice versa */
f2k_parser_fun_t rdmc_f2k_dummy_parser;
f2k_event_writer_fun_t rdmc_f2k_dummy_event_writer;

void rdmc_push_f2k_buffer(rdmc_f2k_buffer_t *b, char *s, const f2000_line_t * type_def);
void rdmc_init_f2k_buffer(rdmc_f2k_buffer_t *b);
void rdmc_reset_f2k_buffer(rdmc_f2k_buffer_t *b);
void rdmc_unlink_f2k_buffer(rdmc_f2k_buffer_t *b);
void rdmc_clear_f2k_buffer(rdmc_f2k_buffer_t *b);

char ** rdmc_f2k_tokenize(char *line, int *nargs);

/* stores specific error information int f2k file pointer for print */
void rdmc_f2k_errorlog(mcfile *fp);

/* parser routines to catch ?, *, N, NaN, -inf, +inf, na */
int rdmc_amanda_strtoi(const char *str, int default_nan);
double rdmc_amanda_strtof(char *str, double default_nan);
char * rdmc_amanda_itostr(int val, int default_nan);
char * rdmc_amanda_ftostr(double val, double default_nan);

/* parser to catch some non float values in special blocks */
double rdmc_amanda_sptof(char *str);

/* parses the status flags at the end of the HT line */
int rdmc_amanda_mhit_stat(mhit_stat_t *hstat,char *stat_s);

/* return the rdmc value for "up" and "dn" .. and vice versa */
float rdmc_amanda_fori(char *ori);
char * rdmc_amanda_sori(float ori);

/* en/decode the time in a sec.ns fomat */
int rdmc_amanda_strtimetoi(const char *stime, int *sec, int *nsec); /* returns 0 on succes */
char * rdmc_amanda_itimetostr(int sec, int nsec);

/* create a string the USES lines for itoken */
char * rdmc_amanda_uses_to_str(int n_uses, mevt_uses_t *uses, int id);

/* tries to patches any year to avoid y2k problems returns corrected year */
int rdmc_f2k_y2k(int year);
/* ******************************************************************* */
/* ******************************************************************* */
/* Functions in f2k_YxZ.c */


f2k_event_reader_fun_t rdmc_rhd_f2k_1x1;
f2k_event_reader_fun_t rdmc_mhead_f2k_1x1;
f2k_event_reader_fun_t rdmc_mevt_f2k_1x1;
f2k_event_reader_fun_t rdmc_mfoot_f2k_1x1;


f2k_parser_fun_t  rdmc_amanda_HI_1x1;
f2k_parser_fun_t  rdmc_amanda_ARRAY_1x1;
f2k_parser_fun_t  rdmc_amanda_FIT_DEF_1x1;
f2k_parser_fun_t  rdmc_amanda_STAT_DEF_1x1;
f2k_parser_fun_t  rdmc_amanda_USER_DEF_1x1;
f2k_parser_fun_t  rdmc_amanda_TRIG_DEF_1x1;
f2k_parser_fun_t  rdmc_amanda_TRIG_PAR_1x1;
f2k_parser_fun_t  rdmc_amanda_OM_1x1;
f2k_parser_fun_t  rdmc_amanda_KH_1x1;
f2k_parser_fun_t  rdmc_amanda_KADC_1x1;
f2k_parser_fun_t  rdmc_amanda_KTDC_1x1;
f2k_parser_fun_t  rdmc_amanda_KTOT_1x1;
f2k_parser_fun_t  rdmc_amanda_KUTC_1x1;
f2k_parser_fun_t  rdmc_amanda_TBEGIN_1x1;
f2k_parser_fun_t  rdmc_amanda_EM_1x1;
f2k_parser_fun_t  rdmc_amanda_TR_1x1;
f2k_parser_fun_t  rdmc_amanda_CH_1x1;
f2k_parser_fun_t  rdmc_amanda_HT_1x1;
f2k_parser_fun_t  rdmc_amanda_US_1x1;
f2k_parser_fun_t  rdmc_amanda_STATUS_1x1;
f2k_parser_fun_t  rdmc_amanda_trigblock_1x1;
f2k_parser_fun_t  rdmc_amanda_TRIG_1x1;
f2k_parser_fun_t  rdmc_amanda_fitblock_1x1;
f2k_parser_fun_t  rdmc_amanda_FIT_1x1;
f2k_parser_fun_t  rdmc_amanda_FRESULT_1x1;
f2k_parser_fun_t  rdmc_amanda_USES_1x1;

f2k_parser_fun_t  rdmc_amanda_HT_1x2;
f2k_parser_fun_t  rdmc_amanda_WF_1x2;

f2k_parser_fun_t  rdmc_amanda_WF_2004x1;

/* ******************************************************************* */
/* write functions */
/* ******************************************************************* */
/* Functions in f2k_YxZ.c */


/* functions handle all writing of 1.1 and 1.2 and 2004.1 */
f2k_event_writer_fun_t rdmc_whead_f2k_1x1_2;
f2k_event_writer_fun_t rdmc_wfoot_f2k_1x1_2;
f2k_event_writer_fun_t rdmc_wevt_f2k_1x1_2;

int rdmc_wrhist_f2k_1x1(const mcfile *fp, const char *s, const char *pre);
int rdmc_wrcomment_f2k_1x1(const mcfile *fp, const char *s);


/* ******************************************************************* */
/* Line defs  */
/* ******************************************************************* */

extern const f2000_line_t V_line_1x1; 
extern const f2000_line_t COMMENT_line_1x1; 
extern const f2000_line_t HI_line_1x1; 
extern const f2000_line_t ARRAY_line_1x1; 
extern const f2000_line_t TRIG_DEF_line_1x1; 
extern const f2000_line_t TRIG_PAR_line_1x1; 
extern const f2000_line_t FIT_DEF_line_1x1; 
extern const f2000_line_t STAT_DEF_line_1x1; 
extern const f2000_line_t USER_DEF_line_1x1; 
extern const f2000_line_t KH_line_1x1; 
extern const f2000_line_t OM_line_1x1; 
extern const f2000_line_t KADC_line_1x1; 
extern const f2000_line_t KTDC_line_1x1; 
extern const f2000_line_t KTOT_line_1x1; 
extern const f2000_line_t KUTC_line_1x1; 
extern const f2000_line_t TBEGIN_line_1x1; 
extern const f2000_line_t EM_line_1x1; 
extern const f2000_line_t TR_line_1x1; 
extern const f2000_line_t CH_line_1x1; 
extern const f2000_line_t HT_line_1x1; 
extern const f2000_line_t US_line_1x1; 
extern const f2000_line_t STATUS_line_1x1; 
extern const f2000_line_t FIT_block_1x1; 
extern const f2000_line_t TRIG_block_1x1; 
extern const f2000_line_t TRIG_line_1x1; 
extern const f2000_line_t FIT_line_1x1; 
extern const f2000_line_t FRESULT_line_1x1; 
extern const f2000_line_t USES_line_1x1; 
extern const f2000_line_t EE_line_1x1; 
extern const f2000_line_t TEND_line_1x1; 
extern const f2000_line_t END_line_1x1; 

extern const f2000_line_t HT_line_1x2; 
extern const f2000_line_t WF_line_1x2; 

extern const f2000_line_t WF_line_2004x1; 

#endif



