// Translation of Fortran commons from the Herwig6
// f77 program into c++ structures to be used in ROOT
// and declaration of Fortran functions as extern
// C functions to be called from the class Herwig6
// author: j. g. contreras jgcn@moni.mda.cinvestav.mx
// date: december 22, 2000

struct dcpx {double dr,di;};

int const NMXHEP = 2000;

struct Hepevt_t {
  int NEVHEP;
  int NHEP;
  int ISTHEP[NMXHEP];
  int IDHEP[NMXHEP];
  int JMOHEP[NMXHEP][2];
  int JDAHEP[NMXHEP][2];
  double PHEP[NMXHEP][5];
  double VHEP[NMXHEP][4];
};

struct Hwbeam_t {
  int IPART1;
  int IPART2;
};

struct Hwbmch_t {
  char PART1[8];
  char PART2[8];
};

struct Hwproc_t {
  double EBEAM1;
  double EBEAM2;
  double PBEAM1;
  double PBEAM2;
  int    IPROC;
  int    MAXEV;
};

struct Hwpram_t {
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
};

struct Hwprch_t {
  char AUTPDF[2][20];
  char BDECAY[4];
};

int const NMXPAR = 500;

struct Hwpart_t {
  int  NEVPAR;
  int  NPAR;
  int  ISTPAR[NMXPAR];
  int  IDPAR[NMXPAR];
  int  JMOPAR[NMXPAR][2];
  int  JDAPAR[NMXPAR][2];
  double  PPAR[NMXPAR][5];
  double  VPAR[NMXPAR][4];
};

struct Hwparp_t {
  double DECPAR[NMXPAR][2];
  double PHIPAR[NMXPAR][2];
  double RHOPAR[NMXPAR][2];
  int TMPAR[NMXPAR];
};

int const MODMAX = 5;

struct Hwbosc_t {
  double  ALPFAC;
  double  BRHIG[12];
  double  ENHANC[12];
  double  GAMMAX;
  double  RHOHEP[NMXHEP][3];
  int     IOPHIG;
  int     MODBOS[MODMAX];
};

struct Hwparc_t {
  int     JCOPAR[NMXPAR][4];
};

struct Hwbrch_t {
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
};

struct Hwevnt_t {
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
};

struct Hwhard_t {
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
};

int const NMXRES = 500;

struct Hwprop_t {
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
};

struct Hwunam_t {
  char RNAME[NMXRES+1][8];
  char TXNAME[NMXRES+1][2][37];
};

int const NMXDKS = 4000;
int const NMXMOD = 200;

struct Hwupdt_t {
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
};

struct Hwuwts_t {
  double REPWT[5][4][4];
  double SNGWT;
  double DECWT;
  double QWT[3];
  double PWT[12];
  double SWTEF[NMXRES];
};

int const NMXCDK = 4000;

struct Hwuclu_t {
  double CLDKWT[NMXCDK];
  double CTHRPW[12][12];
  double PRECO;
  double RESN[12][12];
  double RMIN[12][12];
  int    LOCN[12][12];
  int    NCLDK[NMXCDK];
  int    NRECO;
  int    CLRECO;
};

struct Hwdist_t {
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
};

int const NMXQDK=20;

struct Hwqdks_t {
  double VTXQDK[NMXQDK][4];
  int    IMQDK[NMXQDK];
  int    LOCQ[NMXQDK];
  int    NQDK;
};

int const NMXSUD = 1024;

struct Hwusud_t {
  double ACCUR;
  double QEV[6][NMXSUD];
  double SUD[6][NMXSUD];
  int    INTER;
  int    NQEV;
  int    NSUD;
  int    SUDORD;
};

struct Hwsusy_t {
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
};

struct Hwrpar_t {
  double LAMDA1[3][3][3];
  double LAMDA2[3][3][3];
  double LAMDA3[3][3][3];
  int    HRDCOL[5][2];
  int    RPARTY;
  int    COLUPD;
};

struct Hwminb_t {
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
};

int const NMXCL = 500;

struct Hwclus_t {
  double PPCL[NMXCL][5];
  int    IDCL[NMXCL];
  int    NCL;
};

// herwig 6.507

struct Hwgrav_t {
  double GRVLAM;
  double EMGRV;
  double GAMGRV;
};

struct Hw6202_t {
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
};

struct Hw6203_t {
  double ABWGT;
  double ABWSUM;
  double AVABW;
  int    NNEGWT;
  int    NNEGEV;
  int    NEGWTS;
};

int const IMAXCH = 20;

struct Hw6300_t {
  double MJJMIN;
  double CHNPRB[IMAXCH];
  int    IOPSTP;
  int    IOPSH;
  int    OPTM;
  int    CHON[IMAXCH];
};


int const NXMRS = 49;
int const NQMRS = 37;
int const NPMRS = 8;

struct Hwpmrs_t {
  double FMRS[NQMRS+1][NXMRS][NPMRS][3];
};

struct Hwcirc_t {
  int    CIRCOP;
  int    CIRCAC;
  int    CIRCVR;
  int    CIRCRV;
  int    CIRCCH;
};

int const NCFMAX = 3;
int const NMODE2 = 500;
int const NMODE3 = 500;
int const NDIAGR = 8;
int const NMODEB = 50;
int const NMODE4 = 4;

struct Hwdspb_t {
  double ABMODE[NMODEB][2];
  double BBMODE[NMODEB][12][2];
  double PBMODE[NMODEB][12];
  double WTBMAX[NMODEB][12];
  int    IDBPRT[NMODEB];
  int    IBDRTP[NMODEB];
  int    IBMODE[NMODEB];
  int    NBMODE;
};

struct Hwdsp2_t {
  double A2MODE[NMODE2][2];
  double P2MODE[NMODE2];
  double WT2MAX[NMODE2];
  int    ID2PRT[NMODE2];
  int    I2DRTP[NMODE2];
  int    N2MODE;
};

struct Hwdsp3_t {
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
};

struct Hwdsp4_t {
  double A4MODE[NMODE4][12][2];
  double B4MODE[NMODE4][12][2];
  double P4MODE[NMODE4][12][12];
  double WT4MAX[NMODE4][12][12];
  int    ID4PRT[NMODE4];
  int    I4MODE[NMODE4][2];
  int    N4MODE;
};

struct Hwdspn_t {
  int    NDECSY;
  int    NSEARCH;
  int    LRDEC;
  int    LWDEC;
  int    SYSPIN;
  int    THREEB;
  int    FOURB;
  char   TAUDEC[6];
};

int const NMXSPN = 50;

struct Hwspin_t {
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
};

struct Hwstau_t {
  int    JAK1;
  int    JAK2;
  int    ITDKRC;
  int    IFPHOT;
};

int const MAXHRP = 100;

struct Hwgupr_t {
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
};

struct Hw6500_t {
  int    PRESPL;
};

struct Hw6504_t {
  int    ITOPRD;
};

struct Hw6506_t {
  double PDFX0;
  double PDFPOW;
};



extern "C" {
  void  hwigin_();
  void  hwuinc_();
  void  hwusta_(char * name, int);
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
}

// subroutines to be call by JIMMY

extern "C" {
  void  jminit_();
  void  jimmin_();
  void  jmefin_();
}






