/****************************************************************************/
/* This is rdmc.h - some routines for reading mc and preprocessed data      */
/* Here are the following units used: meter, nanosecond, radian, MeV        */
/****************************************************************************/
/* Version: $Header: /net/local/cvsroot/siegmund/rdmc/rdmc.h,v 1.188 2004/04/21 14:37:46 wiebusch Exp $ */

#ifndef _RDMC_H
#define _RDMC_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <float.h>
////#include <malloc.h>

#define NEW_USES 0
#define NEW_MHEAD 1

/*****************************/
/* general definitions - *******/
/*****************************/

/* limits */

#define RDMC_MAXCHANNELS 16384
#define MAXCHANNELS RDMC_MAXCHANNELS /* OBSOLETE */

#define RDMC_MAXLINE 10000  /* maximal input line for the asccii format */
#define F2000_MAXLINE 10000

#define RDMC_MAXDIGI 4096              /* max. number of digitizations */

#define RDMC_MAXTRIGID  8*(sizeof(long)) /* size of a bitmask */

#define RDMC_MAXTOKENLENGTH 64    /* maximal token length in f2000 */
  /* #define RDMC_MAXTOKEN_PER_LINE 250*/  /* maximum user words per line */
#define RDMC_MAXTOKEN_PER_LINE 1000  /* maximum user words per line */
#define RDMC_MAXUSER_PER_LINE RDMC_MAXTOKEN_PER_LINE  /* maximum user words per line */


/* special; values */
#define RDMC_BIG 1e6                              /* a big value (in meter) */
#define RDMC_SMALL 1e-9                                    /* a small value */
#define BIG RDMC_BIG  /*obsolete */
#define SMALL RDMC_SMALL /* obsolete */

/* some special values converting the f2000 */
#define RDMC_NA -1          /* "not available" number (a quick hack) */
#define RDMC_NOT_AV RDMC_NA                  /* an alias, dont use it */
#define RDMC_PARENT_NA -1 /* this is the parent id for noise hits */
#define RDMC_PARENT_NOISE -2 /* this is the parent id for noise hits */
#define RDMC_PARENT_AFPULS -3 /* this is the parent id for afterpulse hits */	
#define RDMC_REPEAT_CH -2  /* no adc value because it was already given away */
#define RDMC_TDC_NA  -RDMC_BIG 
#define RDMC_SPECIAL_NA -RDMC_BIG
#define RDMC_WF_BASELINE_NA  -RDMC_BIG
#define RDMC_RSORTAMP_ALL 0
#define RDMC_RSORTAMP_FIRST 1
#define RDMC_RSORTAMP_TOT 2

/* special nan values: */
#define RDMC_LONGI_NA RDMC_BIG   /* ARRAY lattitude */
#define RDMC_LATTI_NA RDMC_BIG   /* ARRAY longitude */
#define RDMC_DEPTH_NA RDMC_BIG   /* ARRAY depth longitude */

#define RDMC_DINFTY DBL_MAX         /* double infinity */
#define RDMC_FINFTY FLT_MAX         /* singkle infinity */
#define RDMC_IINFTY INT_MAX         /* integer infinity */
#define RDMC_LINFTY LONG_MAX         /* long integer infinity */

#define RDMC_SMALL_THRESH 0.01
#define RDMC_DEFAULT_DEPTH -1730.00
#define RDMC_MIN_ION_ENERGY 100e3

/****** error messages */
enum RDMC_ERROR_T {
  RDMC_EOF=EOF,             /* end of file from <stdio.h>, usually -1 */
  RDMC_ILF=-10 ,                      /* error value: invalid format! */
  RDMC_UNKNOWN_FORMAT=-11,            /* format could not be detected */
  RDMC_UNKNOWN_FORMAT_VERSION=-12,  /* format version is not suported */
  RDMC_INCONSISTENT_GEOMETRY=-13, /*something fishy in the file header*/
  RDMC_LINE_NOT_PARSED=-14,              /* input line is not  parsed */
  RDMC_TOO_MANY_CHANNELS=-15,     /* OM number larger RDMC_MAXCHANNEL */
  RDMC_HEADERTAG_NOTFOUND=-16,    /* the current mevt special was not */
                                           /*   defined in the header */
  RDMC_EVENT_NOT_RECOGNIZED=-17,  /* the f2k block was not identified */
  RDMC_EE_NOT_RECOGNIZED=-18,   /* the end of an f2k block  not found */
  RDMC_LINE_EMPTY=-19, /* the line is only white (not really an error */
  RDMC_LIBRARY_ERROR=-20,  /* a programming error in the rdmc library */
  RDMC_IO_OK=0                                    /* everythings fine */
};

/***** format ids */
#define BAIKAL_BIN_F 0          /* format id for the baikal-like bin format */
#define DUMAND_ASCII_F 1      /* format id for the dumand-like ASCII format */
#define UWI_ASCII_F 2      /* format id for the ASCII used by Jacobsen etal */
#define STOCKHOLM_BIN_F 3     /* format id for swedish amanda binary format */
#define AMANDA_ASCII_F 4        /* the proposed AMANDA format (aka "F2000") */
#define F2000_F AMANDA_ASCII_F
#define RDMC_DEFAULT_ASCII_F  AMANDA_ASCII_F

/****** Current implemented format versions ********/

#define DUMAND_ASCII_VERSION  -7                /* ASCII version format id */
#define UWI_ASCII_VERSION 3
#define AMANDA_ASCII_VERSION 2004      /* minor number; major expected to be 1 */
#define AMANDA_ASCII_MINOR 1        /* minor number; major expected to be 1 */

/* baikal format dialects */
#define BAIKAL_MOS -1                       /* id for the moscow bin format */
#define BAIKAL_OLD 0                   /* id for the old zeuthen bin format */
#define BAIKAL_NEW 1                   /* id for the new zeuthen bin format */

#define DEFAULT_UWI_GEO_FILE "amanda.geo"
#define DEFAULT_UWI_GEO_FILE_ENV "RDMC_UWI_GEO_FILE"
#ifdef RDMC_DIRTY_AMANDA_ADC /* very ugly */
#define UWI_CHANNELS_PER_PE 120.         /* number of adc channels per p.e. */
                                                   /* (uwi format) */
#else
#define UWI_CHANNELS_PER_PE 1.         /* number of adc channels per p.e. */
#endif
/* some conventions */

#define RDMC_READ_MODE 0
#define RDMC_WRITE_MODE 1

/****** event types */
enum RDMC_EVENT_T {
  RDMC_EVENT_MUON=1, 
  RDMC_EVENT_HEADER_PREF,
  RDMC_EVENT_HEADER,
  RDMC_EVENT_FOOT
};

  /* detector ids */
enum RDMC_GEO_ID_T {
 AMANDA = 2000,
 AMANDA_A = 2004,
 AMANDA_B_4 = 2006,
 AMANDA_B_10 = 2007,
 AMANDA_B_11 = 2008,
 AMANDA_B_13 = 2009,
 AMANDA_II = 2010,
 AMANDA_KM3 = 2099,
 JULIA = 6091,
 BAIKAL = 1000,
 NT_36 = 1003 ,                    /* array id for the baikal = telescope */
 NT_36s = 1004,
 NT_72 = 1005,
 NT_96 = 1006,
 NT_144 = 1007,
 NT_192 = 1008,
 NT_200 = 1099
};

/* om ids */
#define QUASAR 1000000                                     /* Quasar pmt id */
#define STD_AMANDA_A_OM 2030030
#define STD_AMANDA_B_OM 2020021
#define STD_AMANDA_B_10_OM 2020121
#define AMANDA_HYBRID  2020221
#define RDMC_DEFAULT_OM STD_AMANDA_B_OM

/* calibration status */
  /* the entry array_calib_t is a bitmask of these */
#define RDMC_CALIB_TDC 0x01
#define RDMC_CALIB_ADC 0x02
#define RDMC_CALIB_TOT 0x04
#define RDMC_CALIB_GEO 0x08


/* particle ids */
#define GAMMA 1
#define E_PLUS 2
#define E_MINUS 3
#define NU 4
#define MUON_PLUS 5                                /* particle id for a mu+ */
#define MUON_MINUS 6                               /* particle id for a mu- */
#define MUON MUON_PLUS /* default: mu+ */
#define PI_0 7
#define PI_PLUS 8
#define PI_MINUS 9
#define P_PLUS 14
#define P_MINUS 15
#define TAU_PLUS 33 
#define TAU_MINUS 34
#define TAU TAU_PLUS 

#define NU_EL 201
#define NU_MU 202
#define NU_TAU 203
#define NU_EL_BAR 204
#define NU_MU_BAR 205
#define NU_TAU_BAR 206

#define BREMS 1001                        /* particle id for bremsstrahlung */
#define DELTAE 1002                       /* delta electron */
#define PAIRPROD 1003                     /* pair production */
#define NUCL_INT 1004                     /* nuclear interaction */
#define MU_PAIR 1005                      /* muon pair */
#define HADRONS 1006
#define FIBERLASER 2100       /* Fiber laser (Baikal or Amanda calib laser) */
#define N2LASER 2101                      /* N2 laser firing in water or ice */
#define YAGLASER 2201                    /* YAG laser firing in water or ice */
#define Z_PRIMARY 3000                /* CR primary with Z xxx start here */
#define A_PRIMARY 3500                /* CR primary with A xxx start here */

/* which command to use for retrieving of remote files */
#define ALLOW_REMOTE /* allow writing/reading of remote files */
#define RSH_CMD_ENV "RDMC_RSH"   /* name of the environment variable */

/* sorting tags which are allowed by rdmc_sort_hits */
enum RDMC_SORT_T{
  RDMC_NO_SORT=0,                                     /* no hit sorting */
  RDMC_CH_SORT=1,                              /* sort hits by channels */
  RDMC_TM_SORT=2,                              /* sort hits by hit time */
  RDMC_CZ_SORT=3,                        /* sort hits by hit channel depth */
  RDMC_ID_SORT=4,                        /* sort hits by mhit-id */
  RDMC_MT_SORT=5                        /* sort hits by mt-id */
};

  /* types of possible MCinfo inquireies */
  /* used by rdmc_which_gen in rdmc_genutl.c */

enum RDMC_GENINFO_T {
  MOST_HITS_MUON=1 /* muon responsible for  most hits in the event */
  , HIGHEST_ENERGY_MUON /* the most energetic muon */
  , CLOSEST_MUON  /* the closest track to the detector centre */
  , FIRST_HIT_MUON  /* the muon which makes the first hit */
  , MOST_HITS_SHOWER /* shower responsible for  most hits in the event */
  , HIGHEST_ENERGY_SHOWER /* the most energetic shower  */
  , ANY_PRIMARY    /* the first Primary in the file */
  , CR_PRIMARY        /* the  CR PRIMARY */
  , NEUTRINO_PRIMARY  /* the first(index) neutrino primary in event */
};

 /********************************** data structures *************************/



/************** structures for the maintainance of different formats ********/

typedef struct {
  char parse_line[RDMC_MAXLINE];  /* line buffer during I/O f2k */
  int unread;                       /* this flag is 1 if the line is unread */
  int errline;                  /* file line where a parse error occurs */
  int nolex; /* if on, the  event is not parsed and can  not be used further */
  struct rdmc_f2k_buffer_s *f2k_buffer; /* this is defined in f2k.h */
} f2000_fp_t;

typedef struct {
  int mc_id;                         /* mc id character */
              /* 'U' = unknown, 'D' =data, 'M' = MC, 'G' = event generator */
  char mc_vers[RDMC_MAXLINE];                        /* mc program version */
  int igtrack;                                        /* presently unused -1*/
  int igen;                                 /* id for the generator program */
  int igeo;                           /* identifier for expoeriment type    */
  int daswrun;                           /* Random generator seed           */
  time_t time;                               /* mc creation date (obsolete) */
  int nrun;                                 /*run number */
} dumand_fp_t;

typedef struct {
  int mc_id;                         /* mc id character */
              /* 'U' = unknown, 'D' =data, 'M' = MC, 'G' = event generator */
  char mc_vers[RDMC_MAXLINE];                        /* mc program version */
} uwi_fp_t;

typedef struct {
  int mc;                                          /* 0 - data file, 1 - mc */
  int swap;                       /* 1 - Bytes are swapped, 0 - not swapped */
  int fortran_recs;                    /* was it a formatted fortran file? */
  int nw_rec_arr;    /* number of words in reco sub header (baikal_bin only)*/
  int nw_rec_evt;   /* number of words in reco sub record  (baikal_bin only)*/
  int rec_vers;                    /* reco program version (baikal_bin only)*/
  int min_vers;          /* minimization procedure applied (baikal_bin only)*/
  int enr;                              /* current event number OBSOLETE */
  int mc_id;                         /* mc id character */
              /* 'U' = unknown, 'D' =data, 'M' = MC, 'G' = event generator */
  char mc_vers[RDMC_MAXLINE];                        /* mc program version */
  int igtrack;                                        /* presently unused -1*/
  int igen;                                 /* id for the generator program */
  time_t time;                               /* mc creation date (obsolete) */
  int nrun;                                 /*run number */
} baikal_fp_t;

/* a union type for the mcfile structure */
typedef union {
    f2000_fp_t f2000; 
    dumand_fp_t dum; 
    uwi_fp_t uwi; 
    baikal_fp_t bai; 
} rdmc_special_format_info_t;


typedef struct {
  FILE *fp;                                     /* this is the file pointer */
  int mode;                             /* read or write mode */
  int format;                         /* Format: baikal/siegmund/amanda/... */
  int fpos;                           /* file position (line nr or byte nr) */
  int fmajor;                           /* format version */
  int fminor;                                       /* minor version number */
  int zipped;                                /* == 1 if the file was zipped */
  int errline;                  /* source line where an format error occurs */
  char last_line[RDMC_MAXLINE];     /* buffer for the last input line */
                                      /* set to [0]=0 if nothing inside */
  int sloppy;          /* sloppy access to data -> force events after error */
  char *creator;            /* string that points to the History lines */
  char *comment;           /* string that points to the top comment lines */
  rdmc_special_format_info_t info;  /* format specific stuff */
} mcfile;

/* now the event header structures  for CALIBRATION */

typedef struct {
  float beta_t;                             /* TDC resolution, nsec/channel */
  float t_0;                                            /* time shift, nsec */
  float alpha_t;    /* amplitude dependend time calibration, nsec*sqrt(ADC) */
  float ped;                                /* amplitude pedestal, channels */
  float beta_a;                             /* ADC resolution, p.e./channel */
  float kappa;                              /* amplitude nonlinearity, p.e. */
  float ped_tot;                                 /* tot  pedestal, channels */
  float beta_tot;                           /* TOT resolution, p.e./channel */
  float kappa_tot;                               /* tot  nonlinearity, p.e. */
  int flag;        /* != 0 if this is a valid calib entry -> RDMC_CALIB_... */
} array_calib_t;

typedef struct {
  char utc_src[RDMC_MAXTOKENLENGTH];  /* source of this calibration */
  int secs;                                  /* time correction to.be.  */
  int nsecs; 				  /* applied to the event time */
} array_calib_utc_t;  

typedef struct {
  int geo;                                    /* 1 if calibration was done */
  int adc;                                    /* 1 if calibration was done */
  int tdc;                                    /* 1 if calibration was done */
  int tot;                                    /* 1 if calibration was done */
  int utc;                                    /* 1 if calibration was done */
} array_calib_stat_t;                         /* 1 if calibration was done */

/* now various header definitions */

typedef struct {
  int id;
  char tag[RDMC_MAXTOKENLENGTH]; 
  int nwords;
  int npars;
#if !NEW_USES
  char words[RDMC_MAXTOKEN_PER_LINE][RDMC_MAXTOKENLENGTH];
  char pars[RDMC_MAXTOKEN_PER_LINE][RDMC_MAXTOKENLENGTH];
#else
  char *words[RDMC_MAXTOKENLENGTH];
  char *pars[RDMC_MAXTOKENLENGTH];
#endif
} array_hdef_t;


typedef struct {
  int id;                                                    /* geometry id */
  int nch;                                            /* number of channels */
  int nstr;                                            /* number of strings */
  float longitude;                /* detector location longitude, [degrees] */
  float lattitude;                /* detector location lattitude, [degrees] */
  float depth;                                       /* detector depth, [m] */
  int nrun;                                                   /* run number */
  /* geometry follows */
  int str[RDMC_MAXCHANNELS];                               /* string number */
  int clust[RDMC_MAXCHANNELS];                            /* channel number */
  float x[RDMC_MAXCHANNELS];                      /* x coordinate in meters */
  float y[RDMC_MAXCHANNELS];                      /* y coordinate in meters */
  float z[RDMC_MAXCHANNELS];                      /* z coordinate in meters */
  float costh[RDMC_MAXCHANNELS];  /* cos theta of the channel (+1=up/-1=dn) */
  int type[RDMC_MAXCHANNELS];                                   /* OM type */
  int serial[RDMC_MAXCHANNELS];                 /* serial number of channel */
  float thresh[RDMC_MAXCHANNELS];                  /* thresholds of all pmt */
  float sensit[RDMC_MAXCHANNELS];               /* sensitivities of all pmt */
  /* values of the dynamic allocation */
  int n_trigger;                                   /* number of triggers */
  int n_user;                                /* number of userdef lines */
  int n_fit;                                /* number of fit def lines */
  int n_stat;                 /* number of status defs */
  int n_mcinfo;                                   /* number of mcinfo's */
  /* time of the  file */
  time_t tbegin;                                       /* mc creation date */
  time_t tend;                                          /* mc  end date */
  /* calibration */
  array_calib_stat_t is_calib;          /* which calibrations have been done*/
  array_calib_t cal[RDMC_MAXCHANNELS];        /* time/amplitude calibration */
  array_calib_utc_t cal_utc;            /* utc calibration */
  /* trigger */
  array_hdef_t *def_trig;
  /* user */
  array_hdef_t  *def_user;     /* NULL or an array of user line defs */
  /* fits */
  array_hdef_t  *def_fit;        /* NULL or an array of fit line defs */
  /* status ==  housekeeping */
  array_hdef_t *def_stat;           /* see above */
  /* mcinfo, still a dummy */
  array_hdef_t *def_mcinfo;   /* mcinfos types */
 /* uniqe array id attached by init and clear */
  unsigned long int array_id; 
  /* comments */
  char *comment;          /* comment line: a '\0' terminated string or NULL */
  void *tmp;             /* a stack to hold anything */
} array;               /* geometric and technic array info of the telescope */

#if NEW_MHEAD
  /* now the new arry structure */
typedef struct {
  int id;                                         /* unique channel id */
  int type;                                    /* type of this channel */
  float thresh;                             /* diskriminator threshold */
  array_calib_t cal;            /* time/amplitude calibration constants*/
} mch_def_t;


typedef struct {
  int id;                                                  /* channel id */
  int type;                                                   /* OM type */
  int serial;                                     /* serial number of OM */
  int str;                                              /* string number */ 
  int clust;                                      /* OM number in string */
  int x;                                       /* x coordinate in meters */
  int y;                                       /* y coordinate in meters */
  int z;                                       /* z coordinate in meters */
  int costh;                    /* cos theta  orientation  (+1=up/-1=dn) */
  int sensit;                             /* rel. sensitivity of the pmt */
  int n_ch;                                 /* number of readout channels */
  mch_def_t *ch;                         /* field of nch readoutchannels */
} mom_def_t;

typedef struct {
  unsigned long int array_id;   /* uniqe array id created by init and clear */
  int id;                                                    /* geometry id */
  float longitude;                /* detector location longitude, [degrees] */
  float lattitude;                /* detector location lattitude, [degrees] */
  float depth;                                       /* detector depth, [m] */
  int nstr;                                            /* number of strings */
  int n_om;                                                  /* number of OM */
  int n_trigger;                                      /* number of triggers */
  int n_user;                                         /* number of userdef  */
  int n_fit;                                     /* number of fit def lines */
  int n_stat;                                      /* number of status defs */
  int n_mcinfo;                                       /* number of mcinfo's */
  mom_def_t *om;                                       /* field of nom OM's */
  array_calib_stat_t is_calib;          /* which calibrations have been done*/
  array_calib_utc_t cal_utc;                             /* utc calibration */
  array_hdef_t *def_trig;                           /* n_trigger triggerdefs*/
  array_hdef_t  *def_user;                             /*  n_user line defs */
  array_hdef_t  *def_fit;                                /*  n_fit fit defs */
  array_hdef_t *def_stat;                               /*  n_stat fit defs */
  array_hdef_t *def_mcinfo;                           /*  n_mcinfo fit defs */
  char *comment;          /* comment line: a '\0' terminated string or NULL */
  void *tmp;             /* a stack to hold anything */
} mhead_t;           /* structure for the file header/definitions  */

#endif

/* now the event related structures */

typedef struct {
#if NEW_USES
  int nuses; /* number of hits */
  int *id;   /* the hit id */
#else
  int hitid;  /* id of the hit */
  int useid; /* id of the fit or the trigger */
#endif
} mevt_uses_t;



typedef struct {
  int id;  /* identifier/index of e.g. the trigger */
  int nval;
  float val[RDMC_MAXTOKEN_PER_LINE];  /*the values defined e.g. in the TRIG_DEF line */
#if NEW_USES
  mevt_uses_t uses;
#endif
} mevt_special_t;

typedef struct {
  int id;                                    /* particle id (muon+: id = 5) */
  float e;                                                /* energy  in MeV */
  float x,y,z;                           /* parameter of a track point in m */
  float t;                      /* time in nsec corresponding the point xyz */
  float costh;                           /* theta cosinus of the muon track */
  float phi;                                        /* phi value in radians */
  float px, py, pz;                                   /* direction cosinuus */
  float length;                  /* length of the track, or BIG if infinity */
  int nmuon;                        /* obsolete (baikal only) */
  int parent;                     /*track which is mother of this one */
  int tag;                       /* unique tag for this track */
  float weight;                  /* weight of this track */
  int nuser;                  /* DUMMY presently ! */
  mevt_special_t *user;                                 /* the user words */
  char *comment;          /* comment line: a '\0' terminated string or NULL */
  void *tmp;             /* a stack to hold anything */
} mtrack;                         /* track info from mc, reconstructed, etc */

typedef struct {
  int n_tdc_edges;
  int tdc_flag;
} mhit_stat_t;

typedef struct {
  int ch;                                           /* channel number (0,..)*/
  int str;                                       /* string number (1...str) */
  float t;                                   /* time in nsec (Leading edge) */
  float amp;                                   /* amplitude of each channel */
  float tot;                                  /* time over threshold (nsec) */
  int mt;                                    /* muon id for this hit (time) */
  int ma;                               /* muon id for this hit (amplitude) */
  int id;                                 /* a unique hit id */
  float weight;                  /* weight of this hit */
  mhit_stat_t hstat;       /* some status flags for the channel of this hit */
  int nuser;                  /* DUMMY presently ! */
  mevt_special_t *user;                                 /* the user words */
  char *comment;          /* comment line: a '\0' terminated string or NULL */
  void *tmp;             /* a stack to hold anything */
} mhit;                                                    /* hit structure */



typedef struct {
  int id;                                                 /* a unique WF id */
  int om;                                                      /* OM number */
  int ch;                                               /* digitizer number */
  int pa;                                                      /* parent id */
  int ndigi;                                       /* numer of digitizations*/
  float t_start;                                 /* time of 1. digitization */
  float t_bin;                                        /* width of time bins */
  float baseline;                            /* running mean of baseline */
  float *digi;                                   /* amplitude digitizations */
} waveform;                                           /* waveform structure */


typedef struct {
  int nrun;                                                   /* run number */
  int enr;                                                  /* event number */
  int mjd;                                        /* modified julian day  */
  int secs;                                        /* time in sec after mjd */
  int nsecs;                                        /* Nanosec part of time */
  float t_offset;       /* time offset; all hits should be centered by time */
  float weight;                                     /* weight of this event */
  int ee;                                /* 1: end of event corectly tagged */
  int nch;                                     /* number of channels hitted */
  int nstr;                                     /* number of hitted strings */
  int nhits;                                              /* number of hits */
  int nwf;                              /* number of wafeform digitizations */
  int ntrack;                                        /* number of mc tracks */
  int nfit;                       /* number of fitted tracks AND fitresults */
  int nuser;                                      /* number of user blocks  */
  int nstat;                                     /* number of status blocks */
  int ntrig;                                          /* number of triggers */
  int nmcinfo;                           /* number of mcinfo blocks (dummy) */
#if !NEW_USES
  int ntrig_uses;                       /* number of uses words for triggers*/
  int nfit_uses;                           /* number of uses words for fits */
#endif
  enum RDMC_SORT_T sort_status;   /* sort status of hits, def: RDMC_NO_SORT */
  unsigned long int event_id;  /* uniqe event id attached by init and clear */
  unsigned long trigger;                           /* bitmap of the trigger */
  unsigned long trigid;                       /* bitmap of the trigger id's */
  mhit *h;                     /* pointer to all hits (allocated by revt()) */
  waveform *wf;                    /* pointer to all waveform digitizations */
  mtrack *gen;                  /* pointer to the generated tracks, or NULL */
  mtrack *rec;              /* pointer to the reconstructed tracks, or NULL */
  mevt_special_t *fresult;                       /* the  reco results words */
  mevt_special_t *user;                                   /* the user words */
  mevt_special_t *status;                               /* the status words */
  mevt_special_t *ptrig;                               /* the trigger words */
  mevt_special_t *mcinfo;                       /* the mcinfo words (dummy) */
#if !NEW_USES
  mevt_uses_t *trig_uses;                         /* the trigger uses words */
  mevt_uses_t *fit_uses;                              /* the fit uses words */
#endif
  char *comment;          /* comment line: a '\0' terminated string or NULL */
  void *tmp;             /* a stack to hold anything */
} mevt;                 


/* structure that holds converions between strings and numbers */
typedef struct {
  int id;
  char *name;
} rdmc_idtable_t;

/* some tables defined in rdmc_local.c */
extern const rdmc_idtable_t rdmc_pmt_idtable[];
extern const rdmc_idtable_t rdmc_sphere_idtable[];
extern const rdmc_idtable_t rdmc_datatrans_idtable[];
extern const rdmc_idtable_t rdmc_detector_idtable[];
extern const rdmc_idtable_t  rdmc_particle_idtable[];


/********************************** functions *******************************/

/*###########################################*/
/***** in rdmc.c ********/
/*###########################################*/
void rdmc_write_parameters(mcfile *fp, int argc, char **arg, const char *version);
        /* add the current program name  to the creator's list (HISTORY) */

void rdmc_append_comment(char **comment, const char *src);
/* reallocates comment (it must be NULL, or already malloced) */
/* then appends the string src via strcat */

void rdmc_concat_comment(char **comment, const char *src, const int form);
/* reallocates comment (it must be NULL, or already malloced) */
/* then appends the string src via strcat */
/* prvide with form=mcfile.format, e.g. F2000_F */



/*###########################################*/
/***** in rdmc_mcopen.c ********/
/*###########################################*/
mcfile *rdmc_mcopen(const char *name, const char *mode, int format); /* opens a mc/data file */
int rdmc_mcclose(mcfile *fp);              /* closes a mc/data file */


/*###########################################*/
/**** in rdmc_mcfile.c *****/
/*###########################################*/
void rdmc_init_mcfile(mcfile *fp, int format, int mode, FILE *ft); /* init */
void rdmc_free_mcfile(mcfile *fp); /*free */
int rdmc_mccpfp(const mcfile *src, mcfile *dest);   /* copy run parameters */


#if NEW_MHEAD
/*###########################################*/
/**** in rdmc_mhead.c *****/
/*###########################################*/
void rdmc_init_mhead(mhead_t *a);
void rdmc_clear_mhead(mhead_t *a);
void rdmc_free_mhead(mhead_t *a);

/***** now function for the header_def and header_def_par */
/* functions return 0 on success,  1 else */
int rdmc_mhead_add_user_def(mhead_t *mh, array_hdef_t *user, int iuse);
int rdmc_mhead_del_user_def(mhead_t *mh , int iuse);
int rdmc_mhead_add_stat_def(mhead_t *mh, array_hdef_t *sta, int ista);
int rdmc_mhead_del_stat_def(mhead_t *mh, int ista);
int rdmc_mhead_add_fit_def(mhead_t *mh, array_hdef_t *fd, int ifd);
int rdmc_mhead_del_fit_def(mhead_t *mh, int ifd);
int rdmc_mhead_add_mcinfo_def(mhead_t *mh, array_hdef_t *mc, int imc);
int rdmc_mhead_del_mcinfo_def(mhead_t *mh, int imc);
int rdmc_mhead_add_trigger_def(mhead_t *mh, array_hdef_t *t, int it);
int rdmc_mhead_del_trigger_def(mhead_t *mh, int it);

int rdmc_mhead_add_mom(mhead_t *mh, mom_def_t *om, int ipos);



/*###########################################*/
/**** in rdmc_mom.c *****/
/*###########################################*/
void rdmc_init_mom(mom_def_t *ch);
void rdmc_clear_mom(mom_def_t *ch);
void rdmc_free_mom(mom_def_t *ch);


/*###########################################*/
/**** in rdmc_mch.c *****/
/*###########################################*/

void rdmc_init_mch(mch_def_t *ch);
void rdmc_clear_mch(mch_def_t *ch);
void rdmc_free_mch(mch_def_t *ch);

#endif

/*###########################################*/
/**** in rdmc_array_calib.c *****/
/*###########################################*/

void rdmc_init_array_calib(array_calib_t *c);
void rdmc_clear_array_calib(array_calib_t *c);
void rdmc_free_array_calib(array_calib_t *c);

void rdmc_init_array_calib_stat(array_calib_stat_t *cs);
void rdmc_clear_array_calib_stat(array_calib_stat_t *cs);
void rdmc_free_array_calib_stat(array_calib_stat_t *cs);
void rdmc_init_array_calib_utc(array_calib_utc_t *cu);
void rdmc_clear_array_calib_utc(array_calib_utc_t *cu);
void rdmc_free_array_calib_utc(array_calib_utc_t *cu);


/*###########################################*/
/**** in array2mhead.c *****/
/*###########################################*/
/* target has to be allocated and initilized/cleared (freed memory) */ 
int rdmc_array2mhead(const array *ar, mhead_t *mh);
int rdmc_mhead2array(const mhead_t *mh, array *ar);


/*###########################################*/
/**** in rdmc_array.c *****/
/*###########################################*/
void rdmc_init_array(array *ar);                               /* structures */
void rdmc_free_array(array *ar);
void rdmc_clear_array(array *ar);               /* clear the array structure */
int rdmc_rarr(mcfile *fp, array *ar);  /* read the array info from the begin */
int rdmc_warr(mcfile *fp, const array *ar);     /* writes array info to file */
int rdmc_comp_array(const array *a1, const array *a2);/* compares two geometries */
/* returns 0 if they are equal */
/* else it returns the ored result of: */
#define RDMC_ARRAY_COMP_OK 0x0000            /* all OK*/
#define RDMC_ARRAY_COMP_HEADER 0x01         
#define RDMC_ARRAY_COMP_GEO 0x02
#define RDMC_ARRAY_COMP_OMS 0x04
#define RDMC_ARRAY_COMP_CALIB 0x08
#define RDMC_ARRAY_COMP_TRIGGER 0x10
#define RDMC_ARRAY_COMP_USER 0x20
#define RDMC_ARRAY_COMP_FIT 0x40
#define RDMC_ARRAY_COMP_STATUS 0x80
#define RDMC_ARRAY_COMP_MC 0x100

void rdmc_cp_array(array *out, const array *in);/* copies two geometries */

/***** now function for the header_def and header_def_par */
/* functions return 0 on success,  1 else */
int rdmc_array_add_user_def(array *ar, array_hdef_t *user, int iuse);
int rdmc_array_del_user_def(array *ar , int iuse);
int rdmc_array_add_stat_def(array *ar, array_hdef_t *sta, int ista);
int rdmc_array_del_stat_def(array *ar, int ista);
int rdmc_array_add_fit_def(array *ar, array_hdef_t *fd, int ifd);
int rdmc_array_del_fit_def(array *ar, int ifd);
int rdmc_array_add_mcinfo_def(array *ar, array_hdef_t *mc, int imc);
int rdmc_array_del_mcinfo_def(array *ar, int imc);
int rdmc_array_add_trigger_def(array *ar, array_hdef_t *t, int it);
int rdmc_array_del_trigger_def(array *ar, int it);


  /* temporary until array is really nuked */
#define rdmc_add_user_def rdmc_array_add_user_def
#define rdmc_del_user_def rdmc_array_del_user_def
#define rdmc_add_stat_def rdmc_array_add_stat_def
#define rdmc_del_stat_def rdmc_array_del_stat_def
#define rdmc_add_fit_def rdmc_array_add_fit_def
#define rdmc_del_fit_def rdmc_array_del_fit_def
#define rdmc_add_mcinfo_def rdmc_array_add_mcinfo_def
#define rdmc_del_mcinfo_def rdmc_array_del_mcinfo_def
#define rdmc_add_trigger_def rdmc_array_add_trigger_def
#define rdmc_del_trigger_def rdmc_array_del_trigger_def


/* copies all information for channel in_i from geometry *in  */
/* to channel out_i into geometry *out  */
void rdmc_channel_cp(array  *out,int out_i,const array *in ,int in_i); 

/* in order to create objects one needs unique id's*/
/* these functions return them */


/* this function returns a uniqe id AND builds a uniqe tag string */
/* on basis of tag_root. The resuklting tag is filled into tag */
/* The function should return RDMC_NA if this operation fails, */
/* e.g. if there are too many triggers defined */
/* the tag usually contains the tag_root appended by a number */
/* here the trigger number */
int rdmc_unique_trigger_id(char *tag, const array *ar, const char *tag_root);
int rdmc_unique_user_id(char *tag, const array *ar, const char *tag_root);
int rdmc_unique_stat_id(char *tag, const array *ar, const char *tag_root);


/*###########################################*/
/**** in rdmc_hdef.c *****/
/*###########################################*/

/* add a new initilized header def */
int  rdmc_add_array_hdef(array_hdef_t **list,
			       int *count,
			       array_hdef_t *new_ahdt, int ipos); 
/* remove  header def */
int  rdmc_del_array_hdef(array_hdef_t **list, 
			       int *count, int ipos); 


/* copy two structures */
void  rdmc_cp_hdef(array_hdef_t *out, const array_hdef_t *in); 
/* init the structures rdmc_init_array_hdef*/
void  rdmc_init_array_hdef(array_hdef_t *def); 
void  rdmc_clear_array_hdef(array_hdef_t *def); 
void  rdmc_free_array_hdef(array_hdef_t *def); 

/* compares two structures, retyurns 0 if they are the same or the 
   number of differences */
int rdmc_comp_array_hdef(const array_hdef_t *d1, 
			  const array_hdef_t *d2);
/* search an array of ndef array_hdef_t elements */
/* find the number of a header definition */
/* returns the number  0..(ndef-1) or RDMC_NA */
int rdmc_get_hdef_tag(const array_hdef_t *hd, int ndef, const char *tag); /*according to the tag */
int rdmc_get_hdef_id(const array_hdef_t *hd, int ndef,int id); /*according to an id */

/* get a unique id for a new def object */
int rdmc_unique_hdef_id(const array_hdef_t *p, int ndef);

/* get a unique tag for a new def object. It is the tag_root itself
 * as long it is unique, or the tag_root with a number at the end  */
int rdmc_unique_hdef_tag(char *tag,const array_hdef_t *p, int ndef,const char *tag_root);

/* returns the index of the string token in the list of hdef parameters */
/* if the token is not found RDMC_NA is returned */
int rdmc_token_in_hdef(const array_hdef_t *p, const char *token);


/*###########################################*/
/*********** in rdmc_mevt.c **************/
/*###########################################*/
int rdmc_revt(mcfile *fp, array *ar, mevt *ev);          /* reads an event */
int rdmc_skipevt(mcfile *fp);                /* skips over the next event record */
int rdmc_wevt(mcfile *fp, const mevt *ev, const array *ar);   /* writes an event */

/**** init =  set to default values ****/
/***  clear = free also memory    *****/
void rdmc_init_mevt(mevt *ev);                    /* init the mevt structure */
void rdmc_clear_mevt(mevt *ev);                  /* clear the mevt structure */
void rdmc_free_mevt(mevt *ev);           /* free the internal mevt structure */

/* counts number of hit channels in an event */
/* returns number of hit channels */
int rdmc_count_nch(const mevt *e);
/* counts number of hit strings in an event */
/* returns number of hit strings */
int rdmc_count_nstr(const mevt *e);
/* fill the string number into mhit */
int rdmc_fill_mhit_str(mevt *ev, const array *ar);
/* ( rdmc_fill_mhit_str rdmc_count_nstr, replace the old function ns_evt) */


/* Routine to copy datastructures from one event (*in) to another (*out) */
/* *out has to be allocated by the user before via: */
/*            out = (mevt *) malloc(sizeof(mevt);  */
/* the space required for internal substructures e.g. mtrack,mhit */
/* is mallocated during execution of this routine */
/* return value is 0 */
/* previously allocated memory e.g. in tracks is NOT freed ! */
/* in case you want to free first call: rdmc_clear_mevt(out) first */
int rdmc_cp_mevt(mevt *out, mevt *in);

/* now functions which add specific values */
/* return value is 0 if OK */
/* and 1 on failure */
int rdmc_add_gen(mevt *ev, mtrack *tr, int itrack);
int rdmc_del_gen(mevt *ev , int itrack);
int rdmc_add_fit(mevt *ev, mtrack *ft, mevt_special_t *fr, int ifit);
int rdmc_del_fit(mevt *ev , int ifit);
int rdmc_add_user(mevt *ev, mevt_special_t *us, int iu);
int rdmc_del_user(mevt *ev , int iu);
int rdmc_add_status(mevt *ev, mevt_special_t *fr, int ifr);
int rdmc_del_status(mevt *ev , int ifr);
int rdmc_add_mcinfo(mevt *ev, mevt_special_t *mci, int imci);
int rdmc_del_mcinfo(mevt *ev , int imci);
int rdmc_add_trigger(mevt *ev, mevt_special_t *tr, int it, int id);
int rdmc_del_trigger(mevt *ev , int it, int id);
int rdmc_add_trig_uses(mevt *ev, mevt_uses_t *tu, int ituse);
int rdmc_del_trig_uses(mevt *ev , int ituse);
int rdmc_add_fit_uses(mevt *ev, mevt_uses_t *fu, int ifuse); /*ifuse is iuse number */
int rdmc_del_fit_uses(mevt *ev , int ifuse);
/* delete all uses for a ifit, itrig or ihit */
int rdmc_del_ifit_uses(mevt *ev, int fitid);
int rdmc_del_itrig_uses(mevt *ev, int trigid);
int rdmc_del_ihit_uses(mevt *ev, int hitid);

int rdmc_add_mhit(mevt *ev, mhit *h, int ihit);
int rdmc_del_mhit(mevt *ev , int ihit);

int rdmc_add_WF(mevt *ev, waveform *wf, int iwf);
int rdmc_del_WF(mevt *ev , int iwf);


/* now the generic function to add  initilized header def */
int  rdmc_add_mevt_mtrack(mtrack **list, int *count, mtrack *new_mtrack,
			  int ipos); 
int  rdmc_add_mevt_special(mevt_special_t **list, int *count
			   ,mevt_special_t *new_mst, int ipos); 
int  rdmc_add_mevt_mhit(mhit **list, int *count, mhit *new_mhit, int ipos); 
int  rdmc_add_mevt_uses(mevt_uses_t **list, int *count
			,mevt_uses_t *new_mut, int ipos); 
/* remove   */
int  rdmc_del_mevt_mtrack(mtrack **list, int *count, int ipos); 
int  rdmc_del_mevt_special(mevt_special_t **list, int *count, int ipos); 
int  rdmc_del_mevt_mhit(mhit **list, int *count, int ipos); 
int  rdmc_del_mevt_uses(mevt_uses_t **list, int *count, int ipos); 


/* in order to create objects one needs unique id's*/
/* these functions return them */
int rdmc_unique_mhit_id(const mevt *ev);
int rdmc_unique_gen_id(const mevt *ev);

/* In case of corrupted hit/track ids (<0 or double counts) */
/* this functions assign new ids to all hits (mhit.id) */
/*     or all MC tracks (mtrack.tag)  */
/* if inconsitiencies are detected corresponding references */
/* e.g. uses, mhit.mt , mtrack.parent are deleted or set to RDMC_NA */
/* returns the number of id's which wer inconsistent */
int  rdmc_repair_mhit_id(mevt *ev); 
int  rdmc_repair_gen_id(mevt *ev); 

  /* return the index of the hit hitid */
int rdmc_mevt_ihit(mevt *ev, int hitid);

/*###########################################*/
/*********** in rdmc_mevt_special.c **************/
/*###########################################*/
/* init the special structure to RDMC_NA */
/* maxtoken is for speedup, if avaliable  */
/* - if not, set  maxtoken < 0 and RDMC_MAXTOKEN_PER_LINE is taken*/
void rdmc_init_mevt_special(mevt_special_t *s, int maxtoken); 
void rdmc_clear_mevt_special(mevt_special_t *s, int maxtoken); 
void rdmc_free_mevt_special(mevt_special_t *s); 
/* copies the special structure in to the previusly allocated out */
void rdmc_cp_mevt_special(mevt_special_t *out ,mevt_special_t  *in);
/*###########################################*/
/*********** in rdmc_mevt_uses.c **************/
/*###########################################*/
void rdmc_init_mevt_uses(mevt_uses_t *u);
void rdmc_clear_mevt_uses(mevt_uses_t *u);
void rdmc_free_mevt_uses(mevt_uses_t *u);
/* copies nuses uses_t words from in to the previously allocated out */
void rdmc_cp_mevt_nuses(mevt_uses_t *out, mevt_uses_t *in, int nuses);

/* sorts the uses array for a) uses id and the hit_id */
void rdmc_sort_uses(mevt_uses_t *use, int nuse);

/* store all currently used hits as uses block as the given fit/trig num */
int rdmc_create_fit_uses(mevt *e, int ifit);
int rdmc_create_trig_uses(mevt *e, int itrig);

  /* count the number of channels which are used by fit or trigger i*/
int rdmc_count_fit_nch_uses(mevt *e, array *a, int i);


/*###########################################*/
/*********** in rdmc_mtrack.c **************/
/*###########################################*/
/* init the structure */
void rdmc_init_mtrack(mtrack *tr);                      
void rdmc_clear_mtrack(mtrack *tr);             /* free also */
void rdmc_free_mtrack(mtrack *tr);             /* free only */

/* copies the track *in to *out by overwriting.
 * If in == out, nothing is done.  */
void rdmc_cp_mtrack(mtrack *out, const mtrack *in);
/* calcs the direction cosine from the entries tr.costh, tr.phi */
/*  The entries  tr.px,tr.py,tr.pz are updated          */
void rdmc_tau_tr(mtrack *tr);
/* calcs the track angles costh, phi angle from tr.px,tr.py,tr.pz */
/* the entries tr.costh, tr.phi are updated */
void rdmc_tauinv_tr(mtrack *tr);


/*###########################################*/
/*********** in rdmc_mhit.c **************/
/*###########################################*/
void rdmc_init_mhit(mhit *h);                                   /* several */
void rdmc_init_mhit_stat(mhit_stat_t *s);
void rdmc_clear_mhit(mhit *h);             /* free also */
void rdmc_free_mhit(mhit *h);             /* free also */

void rdmc_merge_hits(mhit *h1, mhit *h2, int ch);
/* merges two hits  *hit_1 *hit_2 of the same channel */
/*     filling them into *hit_1 */
/*  - amplitudes are added */
/*  - the earliest time of hit is taken*/
/*  - mt is used from the earliest hit */
/*  - ma is taken from the hit with the largest amplitude */
/*  - ch and str  remain unchanged */
/*  - since rdmc-1.9 mhit->tot is added substracting the overlap*/
/*               tot = tot_1 + tot_2 - (overlap) */
/* the  amplitude of the merged hit is added only if ch=-1
 * else the amplitude of the channel ch is used.
 * If there is no such channel, a zero amplitude is assumed
 */

void rdmc_cp_mhit(mhit *out, mhit *in);
/* copies the hit *in to *out by overwriting.
 * If in == out, nothing is done.
 */



/*###########################################*/
/*********** in rdmc_WF.c **************/
/*###########################################*/
void rdmc_init_WF(waveform *wf);                 /* several */
void rdmc_clear_WF(waveform *wf);             /* free also */
void rdmc_free_WF(waveform *wf);               /* free also */


void rdmc_cp_WF(waveform *out, waveform *in);
/* copies the waveform *in to *out by overwriting.
 * If in == out, nothing is done.
 */



/*###########################################*/
/**** in rdmc_jk *******/
/*###########################################*/
/* this implements functions that provide compatibility to the old 
   jk structure of rdmc-1.x */
/* the indices of the special structure values in "rdmc-jk" */
#define JK_N_TOKEN  8                 /* number of FRESULT tokens */
#define JK_FITID 0                    /* an id of the fit algorithm */
#define JK_RCHI2 1                    /* chi2/NDF or Likelihood/nch */
#define JK_PROB 2                     /* some reco result parameter */
#define JK_SIGTH 3                    /* some reco result parameter */
#define JK_COVMIN 4                   /* some reco result parameter */
#define JK_COVMAX 5                   /* some reco result parameter */
#define JK_CUTFLAG 6                  /* some reco result parameter */
#define JK_CHI2 7                             /* chi2 or likelihood */

/* initilizes a header definition as  "rdmc-jk"  */
void rdmc_jk_to_fitdef(array_hdef_t *jk_def, int id);

/* init values of fitresult as expected for "rdmc-jk" (similar to init_jk) */
void rdmc_init_fit_jk(mevt_special_t *result,int id);

/* checks if fresult ifit is a rdmc-jk fresult returns 1 if yes else 0*/
int rdmc_is_fresult_jk(const array *ar, const mevt *ev, int ifit);

/* checks if special fresult is defined as rdmc-jk  returns 1 if yes 
   else 0*/
int rdmc_is_this_jk(const array_hdef_t *jk_def, const mevt_special_t *result);

/* returns 1 if fitdef idef (0..ar->n_fit-1) is of jk type */
int rdmc_is_fitdef_jk(const array *ar,int idef);

#if 1
extern const array_hdef_t rdmc_jk_fit_def;
#endif

/*###########################################*/
/********* in rdmc_time.c **********/
/*###########################################*/

time_t rdmc_o_dateconv(int date);  /* converts a YYMMDD date into unix time */
int rdmc_o_rdateconv(time_t time); /* converts unix time into YYMMDD format */

int rdmc_gps_to_mjd(int gpsyear, int gpsday);    /* convert GPS date to mjd */
void rdmc_mjd_to_gps(int mjd, int *gpsyear, int *gpsday);     /* mjd to GPS */
void rdmc_tjd_to_mjd(int tjd, double sec, int *mjd, int *secs, int *ns);
	/* convert truncated julian days to mjd */
void rdmc_mjd_to_tjd(int mjd, int secs, int ns, int *tjd, double *sec); 
	/* convert mjd to tjd */
void rdmc_jd_to_mjd(double jd, int *mjd, int *secs, int *ns);   
        /* convert jd to mjd */
void rdmc_mjd_to_jd(int mjd, int secs, int ns, double *jd);    
        /* convert mjd to jd */
time_t rdmc_mjd_to_unixtime(int mjd, int secs); /* convert mjd in unix secs */
int rdmc_unixtime_to_mjd(time_t unix_time);     /* convert unix secs in mjd */
int rdmc_unixtime_to_mjdsecs(time_t unix_time);    /* unix secs in mjd secs */
void rdmc_jd_to_gmst(double jd, double *gmst); 
       /* Julian Days to Greenwich mean sidereal time */
void rdmc_gmst_to_jd(double gmst, int intjd, double *jd);     
       /* convert gmst to jd */

/*###########################################*/
/********* in rdmc_local.c **********/
/*###########################################*/

char *rdmc_which_format(mcfile *fp);  /* returns string with the format name*/


/*###########################################*/
/********* in amanda.c **********/
/*###########################################*/


/*###########################################*/
/********* in f2k_utl.c **********/
/*###########################################*/


int rdmc_is_f2000_supported(int major, int minor);   /* returns 1 if f2000 */
                                             /*format version is supportet */
/* return the rdmc value for the OM/PMT string and vive versa */
int rdmc_amanda_iomid(const char *str, int ar_id);
int rdmc_amanda_ipmtid(const char *str);
const char *rdmc_amanda_spmtid(int type);
const char *rdmc_amanda_somid(int id);
/* convert rdmc particle ids to f2000 names and vice versa */
int rdmc_amanda_ipartid(const char *name);
const char *rdmc_amanda_spartid(int id);
/* convert rdmc detector ids to f2000 names and vice versa */
const  char * rdmc_amanda_sdet(int geo_id); /* returns a pointer to the name */
int rdmc_amanda_idet(const char *detnam);     /* return the id */
/* WARNING the returned pointers are static local strings, which values may change so do not work on these but imediately copy the contnt to yozur private variable */


/*###########################################*/
/* routines in rdmc_phys.c           */
/*###########################################*/

/*****************************************************************************/
/* global const, variables                                                   */
/*****************************************************************************/

extern const double rdmc_c_vacuum; /* m/nsec; light speed in vacuum */
extern const double rdmc_rez_c_vacuum;        /* 1/c  [ ns/m ] */

#ifndef RDMC_OLD_REFIDX
#define RDMC_OLD_REFIDX 0 /* use the old value of the refraction index */
#endif

#if RDMC_OLD_REFIDX
extern const double c_water;                 /* m/nsec; light speed in water */
extern const double rez_c_water;                        /* 1/c_water [ns/m ] */
extern const double n_water;                    /* refraction index of water */
extern const double cer;                                   /* cerenkov angle */
extern const double tg_cer;                     /* tangens of cerenkov angle */
extern const double sn_cer;                                  /* sin cerencov */
extern const double rez_sn_cer;                            /* 1/sin cerencov */
extern const double cs_cer;                                  /* cos cerencov */
#else
extern const double rdmc_n_water; /* refraction index of water */
extern const double rdmc_n_ice_p;/* phase refraction index of water */
extern const double rdmc_n_ice_g;/* group refraction index of water */
extern const double rdmc_c_water; /* m/nsec;  speed in water */
extern const double rdmc_rez_c_water;/* nsec/m;  speed in water */
extern const double rdmc_c_ice_p; /* m/nsec; speed in water */
extern const double rdmc_rez_c_ice_p;/* nsec/m;  speed in water */
extern const double rdmc_c_ice_g; /* m/nsec;  speed in water */
extern const double rdmc_rez_c_ice_g; /* nsec/m;  speed in water */
extern const double rdmc_cs_cer_wat;       /* cos cerencov */
extern const double rdmc_cs_cer_ice_p;     /* cos cerencov */
extern const double rdmc_cs_cer_ice_g;     /* cos cerencov */
extern const double rdmc_cer_wat;       /* cerenkov angle in rad */
extern const double rdmc_cer_ice_p;   /* cerenkov angle in rad */
extern const double rdmc_cer_ice_g;   /* cerenkov angle in rad */
extern const double rdmc_sn_cer_wat;   /* sin cerencov */
extern const double rdmc_sn_cer_ice_p; /* sin cerencov */
extern const double rdmc_sn_cer_ice_g; /* sin cerencov */
extern const double rdmc_rez_sn_cer_wat;    /* 1/sin cerencov */
extern const double rdmc_rez_sn_cer_ice_p; /* 1/sin cerencov */
extern const double rdmc_rez_sn_cer_ice_g; /* 1/sin cerencov */
extern const double rdmc_tg_cer_wat; /*tan  cer */
extern const double rdmc_tg_cer_ice_p; /*tan cer */
extern const double rdmc_tg_cer_ice_g; /*tan cer */
extern const double rdmc_tg_cer_tkoeff; /* koefficient in rdmc_art()  */
#endif

int rdmc_track_closest(double *dist, double poca_xyz[3], 
		    const mtrack *it, const double x, 
		    const double y, const double z);
/* calculates the distance of closest approach to point of coordinates x,y,z*/
/* it returns the values: perpendicular distance *dist */
/* vector poca_xyz[3] with coordinates of the track point of closest aproach */
/* both dist or xyz may be NULL pointers if you do not want to calculate this*/


int rdmc_track_delta(float *delta, float *dist,
		     const mtrack *tr1, const mtrack *tr2);
/* calculates the angular difference (solid angle) between two tracks */
/* and the closest distance between them */
int rdmc_deltaang(double theta1, double phi1, double theta2, 
			double phi2, double *delta);
/* calculate the angular difference between two directions given by */
/* four angles */



/* calculate the cherenkov arrival times */
double rdmc_art_xyz(const mtrack *tr, double x, double y, double z); 
                                          /* Cerenkov light arrival to point */

double rdmc_art(const array *ar, const mtrack *tr, int i); 
                                       /* Cerenkov light arrival to i-th pmt */

double rdmc_pnt_art(const array *ar, const mtrack *tr, int i);
                                         /* arrival time from a point source */

double rdmc_om_dist(const array *a, int ich1, int ich2); /* distance between */
	                                                 /* 2 channels */

double rdmc_dist(double x1, double y1, double z1,        /* distance between */
		 double x2, double y2, double z2);               /* 2 points */

/* rdmc_angtr() calculates the cosinus of the direction angles between two  */
/*   tracks: returns between 1. (parrallel) and -1. (antiparallel)          */
double rdmc_angtr(const mtrack *t1, const mtrack *t2);

double rdmc_tdist_xyz(const mtrack *it, double x, double y, double z);
                           /* calcs the min distance from a track to a point */
double rdmc_tdist(const array *ar, const mtrack *it, int iom);
                         /* calcs the min distance between a track and an om */
double rdmc_tdist_orig(const mtrack *it); 
            /* calcs the min distance between a track and the origin (0,0,0) */

double rdmc_vdist(const array *ar, const mtrack *it, int iom);
/* calcs the distance between a the vertex of track and a pmt       */

void rdmc_tang(const array *ar, const mtrack *it, int iom, 
	    double *dist, double *csang);


void rdmc_vang(const array *ar, const mtrack *it, int iom, 
	    double *dist, double *csaxis, double *csang);

/*calculates the perpendicular distance, direction angle and  cos(angle) */
/* of incident light if dist, csaxis or csang are NULL, */
/* they are not calculated */


/*###########################################*/
/* routines in rdmc_utl.c           */
/*###########################################*/

void rdmc_add_timeoffset(mevt *e, float offset);
/* adds a time offset to all time related values of an event */
/* this are presently ONLY the hit times, and the fit/gen track times */

/* sort hits in a event by a call to quicksort either by time */
/*  or by channel number */
int rdmc_sort_hits(mevt *ev, const array *ar, enum RDMC_SORT_T method);

/* remove fit no ifit (counting is: 0..nfit-1) from event */
/* including all associated information, like fit_uses */
void rdmc_remove_fit(mevt *e, int ifit);

/* remove_groups() removes all hits which are not from the leading muon */
/* (invert=0) or which are from the leading muon (invert=1) */
/* returns the number of hits after the procedure */
/* Noise is treated according to the last flag */
/* ev->nch and ev->nstr is NOT updated !!! */
int rdmc_remove_groups(mevt *e, int invert,int keep_noise);

/* remove a hit from the event */
/* including all associated information, like fit_uses */
/* ev->nch and ev->nstr is NOT updated !!! */
void rdmc_remove_hit(mevt *e, int ihit);

/* Return mtrack->tag of track-like parent of a track */
int rdmc_tracklike_parent(const mevt *e,const mtrack *t);

/* returns 1 if this is a secondary shower else 0 */
int rdmc_is_secondary(mtrack *t);
int rdmc_is_neutrino(mtrack *t);
int rdmc_is_cr_primary(mtrack *t);

/* returns 1 if this mtrack is a point-like or track-like light source */
int rdmc_is_point_like(const mtrack *t);
int rdmc_is_track_like(const mtrack *t);

/****************************************************************************/
/* checks if this hit is a  ifold  [scip] coincident hit                   */
/* returns 1 if yes, 0if no                                                 */
/****************************************************************************/
int rdmc_is_coinc(mevt *e, const array *ar, int ihit, 
			 int fold, float twin);
int rdmc_is_scoinc(mevt *e, const array *ar, int ihit, 
			 int fold, float twin);

/* calculates the number of neighbouring hits */
/* which are coincident in the time window twin */
/* and are within the radius */
/* returns -1 on error */
int rdmc_neighbours(array *a, mevt *e, int ihit, double radius, double twin);


int rdmc_count_geo(const array *a);
/* counts number of strings in an geometry-array *a and returns that number */

void rdmc_make_str(array *a);   
/* cleans up the geometry data structure *a */
/* the number of strings (a->nstr) is recalculated */
/* and the string number for each chanel (a->str[i]) are  recalculated */

void rdmc_make_clust(array *a); 
/* calculates the iclust number for each channel in the data structure *a */
/* the old value a->clust[i] is deleted  (Baikal specific) */

int rdmc_give_clust(const array *a,int ich);
/* this function calculates the value iclust for channel ich in array *a */
/* the function returns the result iclust or -1 on error */
/*  (Baikal specific) */

double rdmc_count_pe(const mevt *e);
/* calculate total pe in this event */

int rdmc_omhit(float *pesum, int iom, const mevt *event);
/* calculates the sum of PE in OM iom ([0..nom-1]) */
/* return 0 if the om was not hit else 1 */

double rdmc_count_ompe(const mevt *e,const array *a
		       ,float pe_oms[],int hit_oms[]);
/* calculates for each om if it was hit, */
/* and how many PE (ADCs summed over all hits) this are */
/* the  arrays pe_oms[] and hit_oms[] are initlized to 0 */
/* but have to be provided by the user (size: ar->nch) */
/* it fills 1 into hit_om[i] if the channel i was hit */
/*     and fills total pe of all hits in this OM into pe_oms[i] */
/* returns the total number of pe of all OM                     */
/* the function can be called with pe_oms or hit_oms set to NULL */
/* then the corresponding values are not calculated (efficiency) */

int rdmc_centre(const mevt *ev, const array *geom, double rc[3]
		, double amp_pwr);
/* calculate the vector rc to the centre of gravity of the event */
/* amp_pwr is the power by which the amplitude is to be weighted */


/*###########################################*/
/* routines in rdmc_genutl.c           */
/*###########################################*/
  /* returns the index of the track in mevt.gen or RDMC_NA
     if not found */

int rdmc_which_gen(mevt *e, enum RDMC_GENINFO_T selection);

  /* get the leading muon from the tracks */
  /*  OLD routine, kept for backward compatibility */
  /* eqivalent to rdmc_which_gen(e, MOST_HITS_MUON); */
int rdmc_get_leading_muon(mevt *e);  /* OBSOLETE */

/* ************************************************************************ */
/* sort_generated_tracks() sorts all generated tracks with respect to the   */
/*  number of hits caused by the track,  */
/* and removes tracks without hits   if  keep_tracks==0 && keep_all==0 */
/* and keeps all tracks    if   keep_all==1 */
/* and removes secondaries (no muon!) without hits   if  keep_tracks=1 */
/****************************************************************************/
void rdmc_sort_generated_tracks(mevt *e, int keep_tracks, int keep_all);



/*###########################################*/
/* routines in rdmc_sky.c           */
/*###########################################*/
void rdmc_spdetec2radec(int mjd, int sec, int nsec, double th, double ph,
		      double *ra, double *dec); 
void rdmc_radec2spdetec(int mjd, int sec, int nsec, double ra, double dec,
			double *theta, double *phi);
void rdmc_old_spdetec2radec(int mjd, int sec, int nsec, double th, double ph,
		      double *ra2, double *dec2); 
  /* transform south pole detector to equatorial coordinates */
void rdmc_radec2llbb(double ra, double dec, double *ll, double *bb);
  /* transform equatorial to galactic coordinates */
void rdmc_llbb2radec(double ll, double bb, double *ra, double *dec);
  /* transform galactic to equatorial coordinates */
void rdmc_llbb2sllsbb(double ll, double bb, double *sll, double *sbb);
  /* transform galactic to supergalactic coordinates */
void rdmc_sllsbb2llbb(double sll, double sbb, double *ll, double *bb);
  /* transform supergalactic to galactic coordinates */
void rdmc_radec2ecliptic(int mjd, int sec, int nsec, double ra, double dec, 
			 double *lamda, double *beta);
  /* transform equatorial to ecliptic coordinates */
void rdmc_ecliptic2radec(int mjd, int sec, int nsec, double lamda, 
			 double beta, double *ra, double *dec);
  /* transform ecliptic to equatorial coordinates */
void rdmc_lunarradec(int mjd, int sec, int nsec, double *ra, double *dec);
  /* calculate the position of the moon in equatorial coordinates */
void rdmc_solarradec(int mjd, int sec, int nsec, double *ra, double *dec);
  /* calculate the position of the sun in equatorial coordinates */
void rdmc_solar_anomaly(int mjd, int sec, int nsec, 
			double *lamda_s, double *M_s);
  /* calculate coordinates needed for rdmc_solarradec and rdmc_lunarradec */
void rdmc_get_obliquity(int mjd, int sec, int nsec, double *obl);
  /* calculate the obliquity needed for ecliptic coordinates */


/*###########################################*/
/* routines in astrocoord.c OBSOLETE          */
/*###########################################*/
void rdmc_coord(float gpssec, int gpsday, int gpsyear, float th, float ph,
	   float *ra, float *dec, float *ll, float *bb);
  /* OBSOLETE use functuions in rdmc_sky.c */



/*###########################################*/
/* routines in rdmc_clean_hits.c           */
/*###########################################*/

typedef struct {
  int modify_mhit;    /* flag to indicate modification of mhit  during fit*/
} rdmc_clean_hit_t;

/*****************   PS hit  selection structure   ************************/
//                              +--------------------------------+
//                              | new structure for hit cleaning |
//                              +--------------------------------+
typedef struct {
  char cut;
  int ia1,ia2;
  float fa1,fa2,fa3;
} hit_sel_options;
typedef struct {
  int id;
  int modify_mhit;    /* flag to indicate modification of mhit  during fit*/
  hit_sel_options* opt;
} rdmc_clean_hit_ps;
/*****************  * * * * * * * * * * * * * * *  ************************/

extern  rdmc_clean_hit_ps rdmc_clean_hit;

//                              +--------------------------------+
//                              | init. selection options        |
//                              +--------------------------------+
void hit_sel_options_init(rdmc_clean_hit_ps** phit_options);
//                              +--------------------------------+
//                              |  init. new  selection param.   |
//                              +--------------------------------+
void ps_new_sel_opt(rdmc_clean_hit_ps* phit_options);
//                              +--------------------------------+
//                              |  clear selection options       |
//                              +--------------------------------+
void ps_clear_sel_opt(rdmc_clean_hit_ps* phit_options);
//                              +--------------------------------+
//                              |  parse  selection options      |
//                              +--------------------------------+
int ps_parse_sel_opt(const char *optarg,const char *progname
                     ,rdmc_clean_hit_ps* phit_options);
//                              +--------------------------------+
//                              |  print selection options       |
//                              +--------------------------------+
void ps_usage_cleaning();
void ps_print_sel_opt(const char *progname
                      ,rdmc_clean_hit_ps* phit_options);
//                              +--------------------------------+
//                              | do hit cleaning                |
//                              +--------------------------------+
int ps_exec_cleaning(mevt *e, array *a, rdmc_clean_hit_ps* phit_options);


void rdmc_usage_cleaning(void); /* print all cleaning options == help */
int rdmc_parse_cleaning(const char *optarg,const char *progname);  
/* parse the cleaning command line */
void rdmc_print_cleaning(const char *progname);    /* print the cleaning options */

/* make a copy of the hit data structure return 1 on error */
int rdmc_clean_backup_hits(mevt *e);
/* restore the backuped copy  return 1 on error */
int rdmc_clean_restore_hits(mevt *e);

/* function does recalculate mevt.nch and mevt.nstr */
int rdmc_exec_cleaning(mevt *e, array *a);    /* delete hits to be cleaned */

/* clear the internal list -> experts use only */
void rdmc_clear_cleaning_list(void); /* reset/clean cleaning list */

/* rountines do not recalculate mevt.nch and mevt.nstr */
int rdmc_rm_h(array *a, mevt *e, int hit_nr,int nhits);   /* remove one hit */
                                                    /* hit_nr: 0..(nhit-1) */
int rdmc_rm_rnd_h(array *a, mevt *e);             /* remove one random hit */
                                                    /* returns ihit*/
int rdmc_rm_rnd_nh(array *a, mevt *e, int n);        /* remove n random hits */
                                                    /* returns number of hits*/
int rdmc_rm_earliest_h(array *a, mevt *e,int num); 
                                           /* remove num earliest hit*/
                                           /* return number of removed hits */ 
int rdmc_rm_latest_h(array *a, mevt *e,int num);
                                           /* remove num latest hit */
                                           /* return number of removed hits */ 
int rdmc_rm_interval_h(array *a, mevt *e,float before, float after);
                                           /* remove outside time window */
                                           /* return number of removed hits */
int rdmc_rm_amp_h(array *a, mevt *e,float amp_low, float amp_high, 
		  int low_channel, int high_channel);
                                     /*remove hits amp_low > amp > amp_high */
                                           /* return number of removed hits */
int rdmc_rm_tot_h(array *a, mevt *e,float tot_low, float tot_high,
		  int low_channel, int high_channel);
                                     /*remove hits tot_low > tot > tot_high */
                                           /* return number of removed hits */
int rdmc_rm_nocoinc_h(array *a, mevt *e
		      , int coinc, float twin); 
                                         /*remove  coinc-fold  hits*/
                                        /* the time window twin is aplied */
                                        /* return number of removed hits */ 
int rdmc_rm_snocoinc_h(array *a, mevt *e
		       , int coinc, float twin); 
                                         /*remove SKIP coinc-fold  hits*/
                                        /* the time window twin is aplied */
                                        /* return number of removed hits */ 
int rdmc_rm_isolate_h(array *a, mevt *e,float window);
                                             /*time separated hits by window*/
                                           /* return number of removed hits */ 
int rdmc_rm_coinc_h(array *a, mevt *e,float window);
                                             /*time coincident hits in window*/
                                           /* return number of removed hits */ 
int rdmc_rm_xtalk_h(array *a, mevt *e, float window, float adc_min, 
		    float ratio);
  /* remove LE-coincident cross-talk */
int rdmc_rm_xtalk_map_h(array *a, mevt *e, float adc_min, float ratio);
  /* remove x-talk using Klug & Hanson channel maps, Klug's time windows */
                                           /* return number of removed hits */ 
int rdmc_rm_local_h(array *a, mevt *e, double radius, double twin, int nhits);
                          /* local space and time isolated hits are removed */
                              /* at least nhits within radius and time twin */
int rdmc_rm_additional_h(array *a, mevt *e,int from);
                                           /* remove additional hits in 
					     each channel*/
                                           /* return number of removed hits */ 
int rdmc_rm_channel_h(array *a, mevt *e,int ich);
                            /* OBSOLETE remove all hits in channel ich*/
int rdmc_rm_channels_h(array *a, mevt *e,int ich1, int ich2);
                          /*remove all hits from channel ich1 to ich2*/
int rdmc_rm_string_h(array *a, mevt *e,int istr);
                                           /*remove all hits in string istr*/
int rdmc_rm_uses_h(array *a, mevt *e,int itrig);
                                           /* remove hits not used by trigger 
					     itrig*/
                                           /* return number of removed hits */ 
int rdmc_rm_fuses_h(array *a, mevt *e,int ifit);
                                           /* remove hits not used by fit 
					     ifit*/
                                           /* return number of removed hits */ 
int rdmc_rm_imhoff_h (array *a, mevt *e, int n);
					   /* remove the hits with highest/lowest */
					   /* (t-t_bar)*(z-z_bar) */
int rdmc_rm_imhoff_s_h (array *a, mevt *e, double r, int s);
                                           /* remove the hits with a */
                                           /* (t-t_bar)*(z-z_bar)  */
                                           /*  above/above and  below/below*/
                                           /* r*sigma of the (t)(z) */
                                           /* if s>/=/<0 */
int rdmc_rm_early_amp_h(array *a, mevt *e,int num, float amp_high);
                                     /* within first num hits, hits with */
                                   /* amplitude larger than amp are removed */ 
int rdmc_rm_amp_early_h(array *a, mevt *e,int num, float amp_high);
                                     /* within first num hits, hits with */
                                   /* amplitude larger than amp are removed */ 
int rdmc_rm_dt_rho_h(array *a, mevt *e, int ifit, float rho, float tmin, float tmax);                              /* removes  (truncates) hits with */
              /*distance rho and outside time-window tmin,tmax from fit */

int rdmc_rm_inrho_h(array *a, mevt *e, int ifit, float rho);                              /* removes  (truncates) hits with radius smaller than rho from fit */
int  rdmc_rm_tunnel_h(array *a, mevt *e, int mode);
int rdmc_rm_split_h (array *a, mevt *e, int splitorder, int splitmode);

/*###########################################*/
/* routines in pandel_track.c and pandel_point.c      */
/*###########################################*/
/* naming convention: */
/*   _pt_... : track 
     _pp_...:  point

     _td_ : time delay
     _ph_ : phit _pnh_: nohit 

     lg : logarithm

     norm : normalisation = integral  to infinity of unormalize function.
     int : cumulative integral to time t
     diff : d/dt differential

     mpe : multi pe
     psa : poisson saturated

     patched: a patched pandel approach
*/

/* init the physics constants */
/*if one does not which to overwrite a certain value 
 then call with a large negative value (<= RDMC_SPECIAL_NA)*/
typedef struct {
  double td_tau;     /* tau for tracks */
  double td_lam;     /* lambda for tracks */
  double td_att;     /* Absorption lenght */
  double ps_tau;     /* tau for tracks */
  double ps_lam;     /* lambda for tracks */
  double td_sigma;   /* pmt jitter -> patched functions */
  double td_dist_p1; /* scale for the distance */
  double td_dist_p0_cs0; /* const for distance ped (P0) */
  double td_dist_p0_cs1; /* const for distance ped (P1*cs_ori) */
  double td_dist_p0_cs2; /* const for distance ped (P2*cs_ori^2) */
  double td_ph_eps_pe0;
  double td_ph_eps_pe1;
  double td_ph_eps_ori_n0;
  double td_ph_eps_ori_pow;
  double td_ph_dist_a;
  double td_ph_dist_b;
  double td_ph_dist_l;
  double td_ph_dist_e0;
} rdmc_pandel_par_t;


extern const  rdmc_pandel_par_t 
  rdmc_pandel_par_h0, /* No irregulaarities in the hole */
  rdmc_pandel_par_h1, /* hole irregularities scat=100cm */
  rdmc_pandel_par_h2, /* hole irregularities scat=50cm */
  rdmc_pandel_par_h3, /* hole irregularities scat=30cm */
  rdmc_pandel_par_h4; /* hole irregularities scat=10cm */


void rdmc_td_init(rdmc_pandel_par_t *pandel_par);

/* now the functions */
double rdmc_pt_td(double delay, double perp_dist, double cs_ori);
double rdmc_pt_lgtd(double delay, double perp_dist, double cs_ori);

double rdmc_pt_td_norm(double perp_dist, double cs_ori);
double rdmc_pt_lgtd_norm(double perp_dist, double cs_ori);

double rdmc_pt_td_int(double delay, double perp_dist, double cs_ori);
double rdmc_pt_lgtd_int(double delay, double perp_dist, double cs_ori);
double rdmc_pt_td_diff(double delay, double perp_dist, double cs_ori);


double rdmc_pt_td_mpe(double delay,double perp_dist,double cs_ori,double pe);
double rdmc_pt_lgtd_mpe(double delay,double perp_dist,double cs_ori,double pe);

double rdmc_pt_td_mpe_int(double delay,double perp_dist,double cs_ori,double pe);
double rdmc_pt_lgtd_mpe_int(double delay,double perp_dist,double cs_ori,double pe);
double rdmc_pt_td_mpe_diff(double delay,double perp_dist,double cs_ori,double pe);

double rdmc_pt_td_psa(double delay, double perp_dist, double cs_ori, double mean_pe);
double rdmc_pt_lgtd_psa(double delay, double perp_dist, double cs_ori, double mean_pe);

double rdmc_pt_td_psa_int(double delay, double perp_dist, double cs_ori, double mean_pe);
double rdmc_pt_lgtd_psa_int(double delay, double perp_dist, double cs_ori, double mean_pe);
double rdmc_pt_td_psa_diff(double delay, double perp_dist, double cs_ori, double mean_pe);


double rdmc_pt_td_patched(double delay, double perp_dist, double cs_ori);
double rdmc_pt_lgtd_patched(double delay, double perp_dist, double cs_ori);

double rdmc_pt_td_mpe_patched(double delay, double perp_dist, double cs_ori, double pe);
double rdmc_pt_lgtd_mpe_patched(double delay, double perp_dist, double cs_ori, double pe);

double rdmc_pt_td_psa_patched(double delay, double perp_dist, double cs_ori, double mean_pe);
double rdmc_pt_lgtd_psa_patched(double delay, double perp_dist, double cs_ori, double mean_pe);

double rdmc_pt_ph(double perp_dist, double cs_ori, double energy, double sensit);
double rdmc_pt_pnh(double perp_dist, double cs_ori, double energy, double sensit);
double rdmc_pt_ph_dt(double perp_dist, double cs_ori, double energy, double sensit, double tmin, double tmax);
double rdmc_pt_pnh_dt(double perp_dist, double cs_ori, double energy, double sensit, double tmin, double tmax);

double rdmc_pp_td(double delay, double perp_dist, double cs_ori, double cs_axis);
double rdmc_pp_lgtd(double delay, double perp_dist, double cs_ori, double cs_axis);
double rdmc_pp_td_norm(double perp_dist, double cs_ori, double cs_axis);
double rdmc_pp_lgtd_norm(double perp_dist, double cs_ori, double cs_axis);


double rdmc_pp_td_int(double delay, double perp_dist, double cs_ori, double cs_axis);
double rdmc_pp_lgtd_int(double delay, double perp_dist, double cs_ori, double cs_axis);
double rdmc_pp_td_diff(double delay, double perp_dist, double cs_ori, double cs_axis);


double rdmc_pp_td_mpe(double delay,double perp_dist,double cs_ori, double cs_axis,double pe);
double rdmc_pp_lgtd_mpe(double delay,double perp_dist,double cs_ori, double cs_axis,double pe);

double rdmc_pp_td_mpe_int(double delay,double perp_dist,double cs_ori, double cs_axis,double pe);
double rdmc_pp_lgtd_mpe_int(double delay,double perp_dist,double cs_ori, double cs_axis,double pe);
double rdmc_pp_td_mpe_diff(double delay,double perp_dist,double cs_ori, double cs_axis,double pe);

double rdmc_pp_td_psa(double delay, double perp_dist, double cs_ori, double cs_axis, double mean_pe);
double rdmc_pp_lgtd_psa(double delay, double perp_dist, double cs_ori, double cs_axis, double mean_pe);

double rdmc_pp_td_psa_int(double delay, double perp_dist, double cs_ori, double cs_axis, double mean_pe);
double rdmc_pp_lgtd_psa_int(double delay, double perp_dist, double cs_ori, double cs_axis, double mean_pe);
double rdmc_pp_td_psa_diff(double delay, double perp_dist, double cs_ori, double cs_axis, double mean_pe);

double rdmc_pp_td_psd_patched(double delay, double dist, double cs_ori, double cs_axis, double mean_pe);
double rdmc_pp_td_patched_int(double delay, double dist, double cs_ori, double cs_axis);


double rdmc_pp_td_patched(double delay, double perp_dist, double cs_ori, double cs_axis);
double rdmc_pp_lgtd_patched(double delay, double perp_dist, double cs_ori, double cs_axis);

double rdmc_pp_td_mpe_patched(double delay, double perp_dist, double cs_ori, double cs_axis, double pe);
double rdmc_pp_lgtd_mpe_patched(double delay, double perp_dist, double cs_ori, double cs_axis, double pe);

double rdmc_pp_td_psa_patched(double delay, double perp_dist, double cs_ori, double cs_axis, double mean_pe);
double rdmc_pp_lgtd_psa_patched(double delay, double perp_dist, double cs_ori, double cs_axis, double mean_pe);

#if 0 /* Not all functions for point sources are implemented yet */ 
double rdmc_pp_ph(double perp_dist, double cs_axis, double energy, double sensit);
double rdmc_pp_pnh(double perp_dist, double cs_axis, double energy, double sensit);
#endif

/*###########################################*/
/* routines in rdmc_ndirect.c           */
/*###########################################*/
/****************************************************************************/
/* rdmc_get_direct_hits calcultes the number of direct early and delayed hits*/
/*   The number of arguments is variable                                    */
/*   A format string is passed  */
/*      for each additional argument it has the form: t-:t+;f   */
/*       t- and t+ define the time window                        */
/*       f is a format id:                                        */
/*         'c'   asks for a number of channels              */
/*         's'   asks for a number of strings              */
/*         'd'   asks for the smallest perpendicular distance between the */
/*               track and a hit */
/*         'D'   asks for the largest perp distance between track and a hit */
/*         'R'   asks for the average perp distance between track and a hit */
/*         'L'   asks for the length of hits projected on the track */
/*                = Zdist(baikal)  */
/*         'M'   asks for the projected length on the track of all hits */
/*                but one: leaving out the most outlying direct hit  */
/*         'O'   Get the largest polar angle in the projection of hits */
/*               onto the plane perp to the track direction */
/*         'x','y','z' = centre of gravity of direct hits */
/*   example : */
/*     get_direct_hits("-10000.:-5.;c -5.:35.;c -5.:20.;s 75.:10000.;c" */
/*	    ,&nearly,&ndir,&nsdir,&nlate);*/
/* the function returns gdir */
/****************************************************************************/
double rdmc_get_direct_hits(const array *a, const mtrack *tr, mevt *e
			    , const char *fmt, ...);

/*###########################################*/
/* routines in rdmc_math.c           */
/*###########################################*/
double rdmc_lgamma(double x);  /* logarithm of the gamma-function */
double rdmc_gamac(double a, double x); /* incomplete gamma function */
long rdmc_nint(double a);                      /* double to nearest integer */
void rdmc_vecprod(const double a[3], const double b[3], double c[3]);
double rdmc_scalprod(const double a[3], const double b[3]);
/*###########################################*/
/* routines in random.c           */
/*###########################################*/
double rdmc_rand(void);             /* random number s from numerical recipes*/
double rdmc_poidev(double xm);        /* poisson value with mean xm */
double rdmc_gasdev(void);                /* gaussian random number */


/*###########################################*/
/* routines in  messages.c                   */
/*###########################################*/

void rdmc_msgprintf(const char *fmt, ...);
/* msgprintf() prints out a message                                         */
void rdmc_errorprintf(const char *fmt, ...);
/* errorprintf() prints out an error message                                */
void rdmc_warnprintf(const char *fmt, ...);
/* warnprintf() prints out a warning message                                */

/*###########################################*/
/* routines in  rdmc_error.c                   */
/*###########################################*/

void rdmc_err_print(mcfile *fp, enum RDMC_ERROR_T ierr);
/* gives more verbose info after an rdmc error occured */

/*###########################################*/
/* routines in  rdmc_poem.c                   */
/*###########################################*/

const char *rdmc_poem(void);


/*###########################################*/
/* routines for Fortran in f_rdmc.c          */
/*###########################################*/


void frdmc_mcopen_(char *name, int* mode, int *format, int *res);
void frdmc_mcclose_(int *mode, int *res);
void frdmc_rarr_   (int *res);
void frdmc_revt_   (int *res);
void frdmc_skipevt_(int *res);
void frdmc_warr_   (int *res);
void frdmc_wevt_   (int *res);

/* filling adding and ddeletion functions */
void frdmc_fill_gen_(int *itrack, int *result); 
void frdmc_fill_rec_(int *itrack, int *result);
void frdmc_fill_user_(int *iuser, int *result);
void frdmc_fill_usdef_(int *iuser, int *result);
void frdmc_fill_wf_(int *iwf, int *result);
void frdmc_new_gen_(int *itrack, int *result);
void frdmc_new_rec_(int *itrack, int *idef, int *result);
void frdmc_new_user_(int *iuser, int *idef, int *result);
void frdmc_new_trig_(int *itrig, int *idef, int *result);
void frdmc_new_wf_(int *iwf, int *result);
void frdmc_new_usdef_(int *iuser, int *result);
void frdmc_new_trdef_(int *itrig, int *result);
void frdmc_new_ftdef_(int *ifit, int *result);
void frdmc_del_gen_(int *itrack, int *result);
void frdmc_del_user_(int *iuser, int *result);
void frdmc_del_wf_(int *iwf, int *result);

void frdmc_add_history_(char *history, int *result);
void frdmc_add_hdcomment_(char *comment, int *result);
void frdmc_add_evcomment_(char *comment, int *result);

void frdmc_init_dummy_(void);
void frdmc_new_farray_(void);
void frdmc_new_fmevt_(void);
void frdmc_new_hits_(void);

int frdmc_return_eventpointer_(void);

#if 1 /* OBSOLETE */
void frdmc_get_trigs_(int *res, int *MxTrig, int *NTrigpar, int *Trig); 
#endif


/*** some interfaces to rdmc utility functions ****/
int frdmc_track_closest_(double *dist, double xyz[3],
                      double *x, double *y, double *z);
int frdmc_get_leading_muon_(void);
int frdmc_is_secondary_(void);
int frdmc_is_point_like_(void);
int frdmc_is_track_like_(void);

void frdmc_amanda_spartid_(int *id, char *string);

/*###########################################*/
/*****OBSOLETE***************/
/*###########################################*/

#ifdef __cplusplus
}
#endif

#endif /* _RDMC_H */
