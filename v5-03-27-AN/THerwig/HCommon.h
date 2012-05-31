#ifndef HerwigCommon
#define HerwigCommon

#ifndef __CFORTRAN_LOADED
#include <cfortran.h>
#endif

extern "C" {
// Translation of Fortran commons from the Herwig6
// f77 program into c++ structures to be used in ROOT
// and declaration of Fortran functions as extern
// C functions to be called from the class Herwig6
// author: j. g. contreras jgcn@moni.mda.cinvestav.mx
// date: december 22, 2000

struct dcpx {double dr,di;};

int const NMXHEP = 4000;

typedef struct {
  int NEVHEP;
  int NHEP;
  int ISTHEP[NMXHEP];
  int IDHEP[NMXHEP];
  int JMOHEP[NMXHEP][2];
  int JDAHEP[NMXHEP][2];
  double PHEP[NMXHEP][5];
  double VHEP[NMXHEP][4];
} HepevtCommon;

#define HEPEVT COMMON_BLOCK(HEPEVT, hepevt)
  COMMON_BLOCK_DEF(HepevtCommon, HEPEVT);


typedef struct {
  int IPART1;
  int IPART2;
} HwbeamCommon;

#define HWBEAM COMMON_BLOCK(HWBEAM, hwbeam)
  COMMON_BLOCK_DEF(HwbeamCommon, HWBEAM);

typedef struct {
  char PART1[8];
  char PART2[8];
} HwbmchCommon;

#define HWBMCH COMMON_BLOCK(HWBMCH, hwbmch)
  COMMON_BLOCK_DEF(HwbmchCommon, HWBMCH);


typedef struct {
  double EBEAM1;
  double EBEAM2;
  double PBEAM1;
  double PBEAM2;
  int    IPROC;
  int    MAXEV;
} HwprocCommon;

#define HWPROC COMMON_BLOCK(HWPROC, hwproc)
  COMMON_BLOCK_DEF(HwprocCommon, HWPROC);


typedef struct {
  double AFCH[2][16];
  double ALPHEM;
  double B1LIM;
  double BETAF;
  double BTCLM;
  double CAFAC;
  double CFFAC;
  double CLMAX;
  double CLPOW;
  double CLSMR[2];
  double CSPEED;
  double ENSOF;
  double ETAMIX;
  double F0MIX;
  double F1MIX;
  double F2MIX;
  double GAMH;
  double GAMW;
  double GAMZ;
  double GAMZP;
  double GEV2NB;
  double H1MIX;
  double PDIQK;
  double PGSMX;
  double PGSPL[4];
  double PHIMIX;
  double PIFAC;
  double PRSOF;
  double PSPLT[2];
  double PTRMS;
  double PXRMS;
  double QCDL3;
  double QCDL5;
  double QCDLAM;
  double QDIQK;
  double QFCH[16];
  double QG;
  double QSPAC;
  double QV;
  double SCABI;
  double SWEIN;
  double TMTOP;
  double VFCH[2][16];
  double VCKM[3][3];
  double VGCUT;
  double VQCUT;   
  double VPCUT;
  double ZBINM;
  double EFFMIN;
  double OMHMIX;
  double ET2MIX;
  double PH3MIX;
  double GCUTME;
  int    IOPREM;
  int    IPRINT;
  int    ISPAC;
  int    LRSUD;
  int    LWSUD;
  int    MODPDF[2];
  int    NBTRY;
  int    NCOLO;
  int    NCTRY;
  int    NDTRY;
  int    NETRY;
  int    NFLAV;
  int    NGSPL;
  int    NSTRU;
  int    NSTRY;
  int    NZBIN;
  int    IOP4JT[2];
  int    NPRFMT;
  int AZSOFT;
  int AZSPIN;
  int CLDIR[2];
  int HARDME;
  int NOSPAC;
  int PRNDEC;
  int PRVTX;
  int SOFTME;
  int ZPRIME;
  int PRNDEF;
  int PRNTEX;
  int PRNWEB;
} HwpramCommon;

#define HWPRAM COMMON_BLOCK(HWPRAM, hwpram)
  COMMON_BLOCK_DEF(HwpramCommon, HWPRAM);

typedef struct {
  char AUTPDF[2][20];
  char BDECAY[4];
} HwprchCommon;

#define HWPRCH COMMON_BLOCK(HWPRCH, hwprch)
  COMMON_BLOCK_DEF(HwprchCommon, HWPRCH);

int const NMXPAR = 500;

typedef struct {
  int  NEVPAR;
  int  NPAR;
  int  ISTPAR[NMXPAR];
  int  IDPAR[NMXPAR];
  int  JMOPAR[NMXPAR][2];
  int  JDAPAR[NMXPAR][2];
  double  PPAR[NMXPAR][5];
  double  VPAR[NMXPAR][4];
} HwpartCommon;

#define HWPART COMMON_BLOCK(HWPART, hwpart)
  COMMON_BLOCK_DEF(HwpartCommon, HWPART);

typedef struct {
  double DECPAR[NMXPAR][2];
  double PHIPAR[NMXPAR][2];
  double RHOPAR[NMXPAR][2];
  int TMPAR[NMXPAR];
} HwparpCommon;

#define HWPARP COMMON_BLOCK(HWPARP, hwparp)
  COMMON_BLOCK_DEF(HwparpCommon, HWPARP);

int const MODMAX = 5;

typedef struct {
  double  ALPFAC;
  double  BRHIG[12];
  double  ENHANC[12];
  double  GAMMAX;
  double  RHOHEP[NMXHEP][3];
  int     IOPHIG;
  int     MODBOS[MODMAX];
} HwboscCommon;

#define HWBOSC COMMON_BLOCK(HWBOSC, hwbosc)
  COMMON_BLOCK_DEF(HwboscCommon, HWBOSC);

typedef struct {
  int     JCOPAR[NMXPAR][4];
} HwparcCommon;

#define HWPARC COMMON_BLOCK(HWPARC, hwparc)
  COMMON_BLOCK_DEF(HwparcCommon, HWPARC);

typedef struct {
  double ANOMSC[2][2];
  double HARDST;
  double PTINT[2][3];
  double XFACT;
  int    INHAD;
  int    JNHAD;
  int    NSPAC[7];
  int    ISLENT;
  int    BREIT;
  int    FROST;
  int    USECMF;
} HwbrchCommon;

#define HWBRCH COMMON_BLOCK(HWBRCH, hwbrch)
  COMMON_BLOCK_DEF(HwbrchCommon, HWBRCH);

typedef struct {
  double AVWGT;
  double EVWGT;
  double GAMWT;
  double TLOUT;
  double WBIGST;
  double WGTMAX;
  double WGTSUM;
  double WSQSUM;
  int    IDHW[NMXHEP];
  int    IERROR;
  int    ISTAT;
  int    LWEVT;
  int    MAXER;
  int    MAXPR;
  int    NOWGT;
  int    NRN[2];
  int    NUMER;
  int    NUMERU;
  int    NWGTS;
  int    GENSOF;
} HwevntCommon;

#define HWEVNT COMMON_BLOCK(HWEVNT, hwevnt)
  COMMON_BLOCK_DEF(HwevntCommon, HWEVNT);

typedef struct {
  double ASFIXD;
  double CLQ[6][7];
  double COSS;
  double COSTH;
  double CTMAX;
  double DISF[2][13];
  double EMLST;
  double EMMAX;
  double EMMIN;
  double EMPOW;
  double EMSCA;
  double EPOLN[3];
  double GCOEF[7];
  double GPOLN;
  double OMEGA0;
  double PHOMAS;
  double PPOLN[3];
  double PTMAX;
  double PTMIN;
  double PTPOW;
  double Q2MAX;
  double Q2MIN;
  double Q2POW;
  double Q2WWMN;
  double Q2WWMX;
  double QLIM;
  double SINS;
  double THMAX;
  double Y4JT;
  double TMNISR;
  double TQWT;
  double XX[2];
  double XLMIN;
  double XXMIN;
  double YBMAX;
  double YBMIN;
  double YJMAX;
  double YJMIN;
  double YWWMAX;
  double YWWMIN;
  double WHMIN;
  double ZJMAX;
  double ZMXISR;
  int    IAPHIG;
  int    IBRN[2];
  int    IBSH;
  int    ICO[10];
  int    IDCMF;
  int    IDN[10];
  int    IFLMAX;
  int    IFLMIN;
  int    IHPRO;
  int    IPRO;
  int    MAPQ[10];
  int    MAXFL;
  int    BGSHAT;
  int    COLISR;
  int    FSTEVT;
  int    FSTWGT;
  int    GENEV;
  int    HVFCEN;
  int    TPOL;
  int     DURHAM;
} HwhardCommon;

#define HWHARD COMMON_BLOCK(HWHARD, hwhard)
  COMMON_BLOCK_DEF(HwhardCommon, HWHARD);

int const NMXRES = 500;

typedef struct {
  double RLTIM[NMXRES+1];
  double RMASS[NMXRES+1];
  double RSPIN[NMXRES+1];
  int    ICHRG[NMXRES+1];
  int    IDPDG[NMXRES+1];
  int    IFLAV[NMXRES+1];
  int    NRES;
  int    VTOCDK[NMXRES+1];
  int    VTORDK[NMXRES+1];
  int    QORQQB[NMXRES+1];
  int    QBORQQ[NMXRES+1];
} HwpropCommon;

#define HWPROP COMMON_BLOCK(HWPROP, hwprop)
  COMMON_BLOCK_DEF(HwpropCommon, HWPROP);

typedef struct {
  char RNAME[NMXRES+1][8];
  char TXNAME[NMXRES+1][2][37];
} HwunamCommon;

#define HWUNAM COMMON_BLOCK(HWUNAM, hwunam)
  COMMON_BLOCK_DEF(HwunamCommon, HWUNAM);

int const NMXDKS = 4000;
int const NMXMOD = 200;

typedef struct {
  double BRFRAC[NMXDKS];
  double CMMOM[NMXDKS];
  double DKLTM[NMXRES];
  int    IDK[NMXDKS];
  int    IDKPRD[NMXDKS][5];
  int    LNEXT[NMXDKS];
  int    LSTRT[NMXRES];
  int    NDKYS;
  int    NME[NMXDKS];
  int    NMODES[NMXRES];
  int    NPRODS[NMXDKS];
  int    DKPSET;
  int    RSTAB[NMXRES+1];
} HwupdtCommon;

#define HWUPDT COMMON_BLOCK(HWUPDT, hwupdt)
  COMMON_BLOCK_DEF(HwupdtCommon, HWUPDT);


typedef struct {
  double REPWT[5][4][4];
  double SNGWT;
  double DECWT;
  double QWT[3];
  double PWT[12];
  double SWTEF[NMXRES];
} HwuwtsCommon;

#define HWUWTS COMMON_BLOCK(HWUWTS, hwuwts)
  COMMON_BLOCK_DEF(HwuwtsCommon, HWUWTS);


int const NMXCDK = 4000;

typedef struct {
  double CLDKWT[NMXCDK];
  double CTHRPW[12][12];
  double PRECO;
  double RESN[12][12];
  double RMIN[12][12];
  int    LOCN[12][12];
  int    NCLDK[NMXCDK];
  int    NRECO;
  int    CLRECO;
} HwucluCommon;

#define HWUCLU COMMON_BLOCK(HWUCLU, hwuclu)
  COMMON_BLOCK_DEF(HwucluCommon, HWUCLU);


typedef struct {
  double EXAG;
  double GEV2MM;
  double HBAR;
  double PLTCUT;
  double VMIN2;
  double VTXPIP[4];
  double XMIX[2];
  double XMRCT[2];
  double YMIX[2];
  double YMRCT[2];
  int    IOPDKL;
  int    MAXDKL;
  int    MIXING;
  int    PIPSMR;
} HwdistCommon;

#define HWDIST COMMON_BLOCK(HWDIST, hwdist)
  COMMON_BLOCK_DEF(HwdistCommon, HWDIST);

int const NMXQDK=20;

typedef struct {
  double VTXQDK[NMXQDK][4];
  int    IMQDK[NMXQDK];
  int    LOCQ[NMXQDK];
  int    NQDK;
} HwqdksCommon;

#define HWQDKS COMMON_BLOCK(HWQDKS, hwqdks)
  COMMON_BLOCK_DEF(HwqdksCommon, HWQDKS);

int const NMXSUD = 1024;

typedef struct {
  double ACCUR;
  double QEV[6][NMXSUD];
  double SUD[6][NMXSUD];
  int    INTER;
  int    NQEV;
  int    NSUD;
  int    SUDORD;
} HwusudCommon;

#define HWUSUD COMMON_BLOCK(HWUSUD, hwusud)
  COMMON_BLOCK_DEF(HwusudCommon, HWUSUD);

typedef struct {
  double TANB;
  double ALPHAH;
  double COSBPA;
  double SINBPA;
  double COSBMA;
  double SINBMA;
  double COSA;
  double SINA;
  double COSB;
  double SINB;
  double COTB;
  double ZMIXSS[4][4];
  double ZMXNSS[4][4];
  double ZSGNSS[4]; 
  double LFCH[16];
  double RFCH[16];
  double SLFCH[4][16];
  double SRFCH[4][16]; 
  double WMXUSS[2][2];
  double WMXVSS[2][2]; 
  double WSGNSS[2];
  double QMIXSS[2][2][6];
  double LMIXSS[2][2][6];
  double THETAT;
  double THETAB;
  double THETAL;
  double ATSS;
  double ABSS;
  double ALSS;
  double MUSS;
  double FACTSS;
  double GHWWSS[3];
  double GHZZSS[3];
  double GHDDSS[4];
  double GHUUSS[4];
  double GHWHSS[3];
  double GHSQSS[2][2][6][4];
  double XLMNSS;
  double RMMNSS;
  double IMSSM;
  double SENHNC[24];
  double SSPARITY;
  int    SUSYIN;
} HwsusyCommon;

#define HWSUSY COMMON_BLOCK(HWSUSY, hwsusy)
  COMMON_BLOCK_DEF(HwsusyCommon, HWSUSY);

typedef struct {
  double LAMDA1[3][3][3];
  double LAMDA2[3][3][3];
  double LAMDA3[3][3][3];
  int    HRDCOL[5][2];
  int    RPARTY;
  int    COLUPD;
} HwrparCommon;

#define HWRPAR COMMON_BLOCK(HWRPAR, hwrpar)
  COMMON_BLOCK_DEF(HwrparCommon, HWRPAR);

typedef struct {
  double PMBN1;
  double PMBN2;
  double PMBN3;
  double PMBK1;
  double PMBK2;
  double PMBM1;
  double PMBM2;
  double PMBP1;
  double PMBP2;
  double PMBP3;
} HwminbCommon;

#define HWMINB COMMON_BLOCK(HWMINB, hwminb)
  COMMON_BLOCK_DEF(HwminbCommon, HWMINB);

int const NMXCL = 500;

typedef struct {
  double PPCL[NMXCL][5];
  int    IDCL[NMXCL];
  int    NCL;
} HwclusCommon;

#define HWCLUS COMMON_BLOCK(HWCLUS, hwclus)
  COMMON_BLOCK_DEF(HwclusCommon, HWCLUS);

// herwig 6.507

typedef struct Hwgrav_t {
  double GRVLAM;
  double EMGRV;
  double GAMGRV;
} HwgravCommon;

#define HWGRAV COMMON_BLOCK(HWGRAV, hwgrav)
  COMMON_BLOCK_DEF(HwgravCommon, HWGRAV);

typedef struct {
  double VIPWID[3];
  double DXRCYL;
  double DXZMAX;
  double DXRSPH;
  int    WZRFR;
  int    FIX4JT;
  int    IMSSM;
  int    IHIGGS;
  int    PARITY;
  int    LRSUSY;
} Hw6202Common;

#define HW6202 COMMON_BLOCK(HW6202, hw6202)
  COMMON_BLOCK_DEF(Hw6202Common, HW6202);

typedef struct {
  double ABWGT;
  double ABWSUM;
  double AVABW;
  int    NNEGWT;
  int    NNEGEV;
  int    NEGWTS;
} Hw6203Common;

#define HW6203 COMMON_BLOCK(HW6203, hw6203)
  COMMON_BLOCK_DEF(Hw6203Common, HW6203);

int const IMAXCH = 20;

typedef struct {
  double MJJMIN;
  double CHNPRB[IMAXCH];
  int    IOPSTP;
  int    IOPSH;
  int    OPTM;
  int    CHON[IMAXCH];
} Hw6300Common;

#define HW6300 COMMON_BLOCK(HW6300, hw6300)
  COMMON_BLOCK_DEF(Hw6300Common, HW6300);

int const NXMRS = 49;
int const NQMRS = 37;
int const NPMRS = 8;

typedef struct {
  double FMRS[NQMRS+1][NXMRS][NPMRS][3];
} HwpmrsCommon;

#define HWPMRS COMMON_BLOCK(HWPMRS, hwpmrs)
  COMMON_BLOCK_DEF(HwpmrsCommon, HWPMRS);

typedef struct {
  int    CIRCOP;
  int    CIRCAC;
  int    CIRCVR;
  int    CIRCRV;
  int    CIRCCH;
} HwcircCommon;

#define HWCIRC COMMON_BLOCK(HWCIRC, hwcirc)
  COMMON_BLOCK_DEF(HwcircCommon, HWCIRC);

int const NCFMAX = 3;
int const NMODE2 = 500;
int const NMODE3 = 500;
int const NDIAGR = 8;
int const NMODEB = 50;
int const NMODE4 = 4;

typedef struct {
  double ABMODE[NMODEB][2];
  double BBMODE[NMODEB][12][2];
  double PBMODE[NMODEB][12];
  double WTBMAX[NMODEB][12];
  int    IDBPRT[NMODEB];
  int    IBDRTP[NMODEB];
  int    IBMODE[NMODEB];
  int    NBMODE;
} HwdspbCommon;

#define HWDSPB COMMON_BLOCK(HWDSPB, hwdspb)
  COMMON_BLOCK_DEF(HwdspbCommon, HWDSPB);

typedef struct {
  double A2MODE[NMODE2][2];
  double P2MODE[NMODE2];
  double WT2MAX[NMODE2];
  int    ID2PRT[NMODE2];
  int    I2DRTP[NMODE2];
  int    N2MODE;
} Hwdsp2Common;

#define HWDSP2 COMMON_BLOCK(HWDSP2, hwdsp2)
  COMMON_BLOCK_DEF(Hwdsp2Common, HWDSP2);

typedef struct {
  double A3MODE[NMODE3][NDIAGR][2];
  double B3MODE[NMODE3][NDIAGR][2];
  double P3MODE[NMODE3];
  double WT3MAX[NMODE3];
  double SPN3CF[NMODE3][NCFMAX][NCFMAX];
  int    ID3PRT[NMODE3];
  int    I3MODE[NMODE3][NDIAGR];
  int    I3DRTP[NMODE3][NDIAGR];
  int    N3MODE;
  int    NDI3BY[NMODE3];
  int    N3NCFL[NMODE3];
  int    I3DRCF[NMODE3][NDIAGR];
} Hwdsp3Common;

#define HWDSP3 COMMON_BLOCK(HWDSP3, hwdsp3)
  COMMON_BLOCK_DEF(Hwdsp3Common, HWDSP3);

typedef struct {
  double A4MODE[NMODE4][12][2];
  double B4MODE[NMODE4][12][2];
  double P4MODE[NMODE4][12][12];
  double WT4MAX[NMODE4][12][12];
  int    ID4PRT[NMODE4];
  int    I4MODE[NMODE4][2];
  int    N4MODE;
} Hwdsp4Common;

#define HWDSP4 COMMON_BLOCK(HWDSP4, hwdsp4)
  COMMON_BLOCK_DEF(Hwdsp4Common, HWDSP4);

typedef struct {
  int    NDECSY;
  int    NSEARCH;
  int    LRDEC;
  int    LWDEC;
  int    SYSPIN;
  int    THREEB;
  int    FOURB;
  char   TAUDEC[6];
} HwdspnCommon;

#define HWDSPN COMMON_BLOCK(HWDSPN, hwdspn)
  COMMON_BLOCK_DEF(HwdspnCommon, HWDSPN);

int const NMXSPN = 50;

typedef struct {
  dcpx   MESPN[NMXSPN][NCFMAX][2][2][2][2];
  dcpx   RHOSPN[NMXSPN][2][2];
  double SPNCFC[NMXSPN][NCFMAX][NCFMAX];
  int    IDSPN[NMXSPN];
  int    JMOSPN[NMXSPN];
  int    JDASPN[NMXSPN][2];
  int    NSPN;
  int    ISNHEP[NMXHEP];
  int    NSNTRY;
  int    DECSPN[NMXSPN];
  int    NCFL[NMXSPN];
  int    SPCOPT;
} HwspinCommon;

#define HWSPIN COMMON_BLOCK(HWSPIN, hwspin)
  COMMON_BLOCK_DEF(HwspinCommon, HWSPIN);

typedef struct {
  int    JAK1;
  int    JAK2;
  int    ITDKRC;
  int    IFPHOT;
} HwstauCommon;

#define HWSTAU COMMON_BLOCK(HWSTAU, hwstau)
  COMMON_BLOCK_DEF(HwstauCommon, HWSTAU);

int const MAXHRP = 100;

typedef struct {
  double LHWGT[MAXHRP];
  double LHWGTS[MAXHRP];
  double LHXSCT[MAXHRP];
  double LHXERR[MAXHRP];
  double LHXMAX[MAXHRP];
  double LHMXSM;
  int    LHIWGT[MAXHRP];
  int    LHNEVT[MAXHRP];
  int    ITYPLH;
  int    LHSOFT;
  int    LHGLSF;    
} HwguprCommon;

#define HWGUPR COMMON_BLOCK(HWGUPR, hwgupr)
  COMMON_BLOCK_DEF(HwguprCommon, HWGUPR);

typedef struct {
  int    PRESPL;
} Hw6500Common;

#define HW6500 COMMON_BLOCK(HW6500, hw6500)
  COMMON_BLOCK_DEF(Hw6500Common, HW6500);

typedef struct {
  int    ITOPRD;
} Hw6504Common;

#define HW6504 COMMON_BLOCK(HW6504, hw6504)
  COMMON_BLOCK_DEF(Hw6504Common, HW6504);

typedef struct {
  double PDFX0;
  double PDFPOW;
} Hw6506Common;

#define HW6506 COMMON_BLOCK(HW6506, hw6506)
  COMMON_BLOCK_DEF(Hw6506Common, HW6506);

}

extern "C" {
    void  hwuepr_();
    void  hwigin_();
    void  hwuinc_();
    void  hweini_();
    void  hwuine_();
    void  hwepro_();
    void  hwbgen_();
    void  hwdhob_();
    void  hwcfor_();
    void  hwcdec_();
    void  hwdhad_();
    void  hwdhvy_();
    void  hwmevt_();
    void  hwufne_();
    void  hwefin_();
    void  hwusta_(const char * name, int);
    void  hwiodk_(int);
}

// subroutines to be call by JIMMY

extern "C" {
  void  jminit_();
  void  jimmin_();
  void  jmefin_();
}

#endif




