#ifndef THERWIG6_H
#define THERWIG6_H

// declaration of c++ Class THerwig6 to be used in ROOT
// this is a c++ interface to the F77 Herwig6 program
// author: j. g. contreras jgcn@moni.mda.cinvestav.mx
// date: december 22, 2000

/*

 Class THerwig6 is an interface to the Herwig program

C-----------------------------------------------------------------------
C                           H E R W I G
C
C            a Monte Carlo event generator for simulating
C        +---------------------------------------------------+
C        | Hadron Emission Reactions With Interfering Gluons |
C        +---------------------------------------------------+
C I.G. Knowles(*), G. Marchesini(+), M.H. Seymour($) and B.R. Webber(#)
C-----------------------------------------------------------------------
C with Minimal Supersymmetric Standard Model Matrix Elements by
C                  S. Moretti($) and K. Odagiri($)
C-----------------------------------------------------------------------
C R parity violating Supersymmetric Decays and Matrix Elements by
C                          P. Richardson(&)
C-----------------------------------------------------------------------
C matrix element corrections to top decay and Drell-Yan type processes
C                         by G. Corcella(+)
C-----------------------------------------------------------------------
C Deep Inelastic Scattering and Heavy Flavour Electroproduction by
C                  G. Abbiendi(@) and L. Stanco(%)
C-----------------------------------------------------------------------
C and Jet Photoproduction in Lepton-Hadron Collisions by J. Chyla(~)
C-----------------------------------------------------------------------
C(*)  Department of Physics & Astronomy, University of Edinburgh
C(+)  Dipartimento di Fisica, Universita di Milano
C($)  Rutherford Appleton Laboratory
C(#)  Cavendish Laboratory, Cambridge
C(&)  Department of Physics, University of Oxford
C(@)  Dipartimento di Fisica, Universita di Bologna
C(%)  Dipartimento di Fisica, Universita di Padova
C(~)  Institute of Physics, Prague
C-----------------------------------------------------------------------
C                  Version 6.100 - 16th December 1999
C-----------------------------------------------------------------------
C Main reference:
C    G.Marchesini,  B.R.Webber,  G.Abbiendi,  I.G.Knowles,  M.H.Seymour,
C    and L.Stanco, Computer Physics Communications 67 (1992) 465.
C-----------------------------------------------------------------------
C Please send e-mail about  this program  to one of the  authors at the
C following Internet addresses:
C    I.Knowles@ed.ac.uk        Giuseppe.Marchesini@mi.infn.it
C    M.Seymour@rl.ac.uk        webber@hep.phy.cam.ac.uk
C-----------------------------------------------------------------------
*/

/* declarations from ROOT */
#include "TGenerator.h"

typedef enum
{
   kHwCharm         =  1704,
   kHwBeauty        =  1705,
   kHwCharmMCATNLO  = -1704,
   kHwBeautyMCATNLO = -1705,
   kHwJetsMCATNLO   = -1396
} Process_t;

class TObjArray;

// Translation of Fortran commons from the Herwig6
// f77 program into c++ structures to be used in ROOT
// and declaration of Fortran functions as extern
// C functions to be called from the class Herwig6
// author: j. g. contreras jgcn@moni.mda.cinvestav.mx
// date: december 22, 2000

typedef struct {double dr,di;} dcpx;

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
} Hepevt_t;

typedef struct {
  int IPART1;
  int IPART2;
} Hwbeam_t;

typedef struct {
  char PART1[8];
  char PART2[8];
} Hwbmch_t;

typedef struct {
  double EBEAM1;
  double EBEAM2;
  double PBEAM1;
  double PBEAM2;
  int    IPROC;
  int    MAXEV;
} Hwproc_t;

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
} Hwpram_t;

typedef struct {
  char AUTPDF[2][20];
  char BDECAY[4];
} Hwprch_t;

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
} Hwpart_t;

typedef struct {
  double DECPAR[NMXPAR][2];
  double PHIPAR[NMXPAR][2];
  double RHOPAR[NMXPAR][2];
  int TMPAR[NMXPAR];
} Hwparp_t;

int const MODMAX = 50;

typedef struct {
  double  ALPFAC;
  double  BRHIG[12];
  double  ENHANC[12];
  double  GAMMAX;
  double  RHOHEP[NMXHEP][3];
  int     IOPHIG;
  int     MODBOS[MODMAX];
} Hwbosc_t;

typedef struct {
  int     JCOPAR[NMXPAR][4];
} Hwparc_t;

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
} Hwbrch_t;

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
  int    EV1PR;
  int    EV2PR;
} Hwevnt_t;

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
} Hwhard_t;

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
} Hwprop_t;

typedef struct {
  char RNAME[NMXRES+1][8];
  char TXNAME[NMXRES+1][2][37];
} Hwunam_t;

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
} Hwupdt_t;

typedef struct {
  double REPWT[5][4][4];
  double SNGWT;
  double DECWT;
  double QWT[3];
  double PWT[12];
  double SWTEF[NMXRES];
} Hwuwts_t;

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
} Hwuclu_t;

typedef struct {
  double EXAG;
  double GEV2MM;
  double HBAR;
  double PLTCUT;
  double VMIN2;
  double VTXPIP[5];
  double XMIX[2];
  double XMRCT[2];
  double YMIX[2];
  double YMRCT[2];
  int    IOPDKL;
  int    MAXDKL;
  int    MIXING;
  int    PIPSMR;
} Hwdist_t;

int const NMXQDK=20;

typedef struct {
  double VTXQDK[NMXQDK][4];
  int    IMQDK[NMXQDK];
  int    LOCQ[NMXQDK];
  int    NQDK;
} Hwqdks_t;

int const NMXSUD = 1024;

typedef struct {
  double ACCUR;
  double QEV[6][NMXSUD];
  double SUD[6][NMXSUD];
  int    INTER;
  int    NQEV;
  int    NSUD;
  int    SUDORD;
} Hwusud_t;

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
  double DMSSM;
  double SENHNC[24];
  double SSPARITY;
  int    SUSYIN;
} Hwsusy_t;

typedef struct {
  double LAMDA1[3][3][3];
  double LAMDA2[3][3][3];
  double LAMDA3[3][3][3];
  int    HRDCOL[5][2];
  int    RPARTY;
  int    COLUPD;
} Hwrpar_t;

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
} Hwminb_t;

int const NMXCL = 500;

typedef struct {
  double PPCL[NMXCL][5];
  int    IDCL[NMXCL];
  int    NCL;
} Hwclus_t;

// herwig 6.507

typedef struct {
  double GRVLAM;
  double EMGRV;
  double GAMGRV;
} Hwgrav_t;

typedef struct {
  double VIPWID[3];
  double DXRCYL;
  double DXZMAX;
  double DXRSPH;
  int	   WZRFR;
  int	   FIX4JT;
  int	   IMSSM;
  int	   IHIGGS;
  int	   PARITY;
  int	   LRSUSY;
} Hw6202_t;

typedef struct {
  double ABWGT;
  double ABWSUM;
  double AVABW;
  int	   NNEGWT;
  int	   NNEGEV;
  int	   NEGWTS;
} Hw6203_t;

int const IMAXCH = 20;

typedef struct {
  double MJJMIN;
  double CHNPRB[IMAXCH];
  int	   IOPSTP;
  int	   IOPSH;
  int	   OPTM;
  int	   CHON[IMAXCH];
} Hw6300_t;


int const NXMRS = 49;
int const NQMRS = 37;
int const NPMRS = 8;

typedef struct {
  double FMRS[NQMRS+1][NXMRS][NPMRS][3];
} Hwpmrs_t;

typedef struct {
  int	   CIRCOP;
  int	   CIRCAC;
  int	   CIRCVR;
  int	   CIRCRV;
  int	   CIRCCH;
} Hwcirc_t;

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
  int	   IDBPRT[NMODEB];
  int	   IBDRTP[NMODEB];
  int	   IBMODE[NMODEB];
  int	   NBMODE;
} Hwdspb_t;

typedef struct {
  double A2MODE[NMODE2][2];
  double P2MODE[NMODE2];
  double WT2MAX[NMODE2];
  int	   ID2PRT[NMODE2];
  int	   I2DRTP[NMODE2];
  int	   N2MODE;
} Hwdsp2_t;

typedef struct {
  double A3MODE[NMODE3][NDIAGR][2];
  double B3MODE[NMODE3][NDIAGR][2];
  double P3MODE[NMODE3];
  double WT3MAX[NMODE3];
  double SPN3CF[NMODE3][NCFMAX][NCFMAX];
  int	   ID3PRT[NMODE3];
  int	   I3MODE[NMODE3][NDIAGR];
  int	   I3DRTP[NMODE3][NDIAGR];
  int	   N3MODE;
  int	   NDI3BY[NMODE3];
  int	   N3NCFL[NMODE3];
  int	   I3DRCF[NMODE3][NDIAGR];
} Hwdsp3_t;

typedef struct {
  double A4MODE[NMODE4][12][2];
  double B4MODE[NMODE4][12][2];
  double P4MODE[NMODE4][12][12];
  double WT4MAX[NMODE4][12][12];
  int	   ID4PRT[NMODE4];
  int	   I4MODE[NMODE4][2];
  int	   N4MODE;
} Hwdsp4_t;

typedef struct {
  int	   NDECSY;
  int	   NSEARCH;
  int	   LRDEC;
  int	   LWDEC;
  int	   SYSPIN;
  int	   THREEB;
  int	   FOURB;
  char   TAUDEC[6];
} Hwdspn_t;

int const NMXSPN = 50;

typedef struct {
  dcpx   MESPN[NMXSPN][NCFMAX][2][2][2][2];
  dcpx   RHOSPN[NMXSPN][2][2];
  double SPNCFC[NMXSPN][NCFMAX][NCFMAX];
  int	   IDSPN[NMXSPN];
  int	   JMOSPN[NMXSPN];
  int	   JDASPN[NMXSPN][2];
  int	   NSPN;
  int	   ISNHEP[NMXHEP];
  int	   NSNTRY;
  int	   DECSPN[NMXSPN];
  int	   NCFL[NMXSPN];
  int	   SPCOPT;
} Hwspin_t;

typedef struct {
  int	   JAK1;
  int	   JAK2;
  int	   ITDKRC;
  int	   IFPHOT;
} Hwstau_t;

int const MAXHRP = 100;

typedef struct {
  double LHWGT[MAXHRP];
  double LHWGTS[MAXHRP];
  double LHXSCT[MAXHRP];
  double LHXERR[MAXHRP];
  double LHXMAX[MAXHRP];
  double LHMXSM;
  int	   LHIWGT[MAXHRP];
  int	   LHNEVT[MAXHRP];
  int	   ITYPLH;
  int	   LHSOFT;
  int	   LHGLSF;
} Hwgupr_t;

typedef struct {
  int    PRESPL;
} Hw6500_t;

typedef struct {
  int    ITOPRD;
} Hw6504_t;

typedef struct {
  double PDFX0;
  double PDFPOW;
} Hw6506_t;




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
  void  hwiodk_(int);
}

// JIMMY4.2
extern "C" {
  void  jimmin_();
  void  jminit_();
  void  jmefin_();
}





/* THerwig6 class declaration */
class THerwig6 : public TGenerator {
//----------------------------------------------------------------------------
//  functions:
//----------------------------------------------------------------------------
public:
				// ****** constructors and destructor
  THerwig6();
  THerwig6(const THerwig6 & source);
  THerwig6 & operator=(const THerwig6 & /*source*/) {
    Fatal("THerwig6","Assignment operator not implemented yet");
    return *this;
  }
  virtual ~THerwig6();

  // acces to hep common block
  Hepevt_t*   GetHepevt        () const     { return fHepevt; }
  int         GetNEVHEP        () const     { return fHepevt->NEVHEP; }
  int         GetNhep          () const     { return fHepevt->NHEP; }
  int         GetISTHEP    (int i)const     { return fHepevt->ISTHEP[i-1]; }
  int         GetIDHEP     (int i)const     { return fHepevt->IDHEP[i-1]; }
  int         GetJMOHEP (int i, int j) const
    { return fHepevt->JMOHEP[i-1][j-1]; }
  int         GetJDAHEP (int i, int j) const
    { return fHepevt->JDAHEP[i-1][j-1]; }
  double      GetPHEP   (int i, int j) const
    { return fHepevt->PHEP[i-1][j-1]; }
  double      GetVHEP   (int i, int j) const
    { return fHepevt->VHEP[i-1][j-1]; }

  // access to Herwig6 common-blocks
  // WARNING: Some arrays start in 1, others in 0. Look up the manual!

  // /HWBEAM/

  Hwbeam_t*   GetHwbeam        ()           { return fHwbeam; }
  int         GetIPART1        () const     { return fHwbeam->IPART1; }
  int         GetIPART2        () const     { return fHwbeam->IPART2; }

  // /HWBMCH/
  Hwbmch_t*   GetHwbmch        ()           { return fHwbmch; }
  char*       GetPART1         () const     { return fHwbmch->PART1; }
  char*       GetPART2         () const     { return fHwbmch->PART2; }


  // /HWPROC/
  Hwproc_t*   GetHwproc        ()           { return fHwproc; }
  double      GetEBEAM1        () const     { return fHwproc->EBEAM1; }
  double      GetEBEAM2        () const     { return fHwproc->EBEAM2; }
  double      GetPBEAM1        () const     { return fHwproc->PBEAM1; }
  double      GetPBEAM2        () const     { return fHwproc->PBEAM2; }
  int         GetIPROC         () const     { return fHwproc->IPROC; }
  int         GetMAXEV         () const     { return fHwproc->MAXEV; }

  // /HWPRAM/
  Hwpram_t*   GetHwpram        ()           { return fHwpram; }
  double      GetQCDLAM        () const     { return fHwpram->QCDLAM; }
  void        SetQCDLAM   (double q)        { fHwpram->QCDLAM = q; }
  double      GetVQCUT         () const     { return fHwpram->VQCUT; }
  void        SetVQCUT    (double v)        { fHwpram->VQCUT = v; }
  double      GetVGCUT         () const     { return fHwpram->VGCUT; }
  void        SetVGCUT    (double v)        { fHwpram->VGCUT = v; }
  double      GetVPCUT         () const     { return fHwpram->VPCUT; }
  void        SetVPCUT    (double v)        { fHwpram->VPCUT = v; }
  double      GetCLMAX         () const     { return fHwpram->CLMAX; }
  void        SetCLMAX    (double c)        { fHwpram->CLMAX = c; }
  double      GetCLPOW         () const     { return fHwpram->CLPOW; }
  void        SetCLPOW    (double c)        { fHwpram->CLPOW = c; }
  double      GetPSPLT    (int i) const     { return fHwpram->PSPLT[i-1];}
  void        SetPSPLT    (int i, double p) { fHwpram->PSPLT[i-1] = p;}
  double      GetQDIQK         () const     { return fHwpram->QDIQK; }
  void        SetQDIQK    (double q)        { fHwpram->QDIQK = q; }
  double      GetPDIQK         () const     { return fHwpram->PDIQK; }
  void        SetPDIQK    (double p)        { fHwpram->PDIQK = p; }
  double      GetQSPAC         () const     { return fHwpram->QSPAC; }
  void        SetQSPAC    (double q)        { fHwpram->QSPAC = q; }
  double      GetPTRMS         () const     { return fHwpram->PTRMS; }
  void        SetPTRMS    (double p)        { fHwpram->PTRMS = p; }
  double      GetENSOF         () const     { return fHwpram->ENSOF; }
  void        SetENSOF    (double e)        { fHwpram->ENSOF = e; }
  int         GetIPRINT        () const     { return fHwpram->IPRINT; }
  void        SetIPRINT   (int i)           { fHwpram->IPRINT = i; }
  int         GetMODPDF   (int i) const     { return fHwpram->MODPDF[i-1];}
  void        SetMODPDF   (int i, int j)  { fHwpram->MODPDF[i-1] = j; }
  int         GetNSTRU         () const     { return fHwpram->NSTRU; }
  void        SetNSTRU    (int i)          { fHwpram->NSTRU = i; }

  // /HWPRCH/
  Hwprch_t*   GetHwprch        ()           { return fHwprch; }
  char*       GetAUTPDF     (int i)         { return fHwprch->AUTPDF[i-1]; }
  void        SetAUTPDF(int i,const char* s){ strncpy(fHwprch->AUTPDF[i-1],s,20);}
  char*       GetBDECAY        ()           { return fHwprch->BDECAY; }

  // /HWPART/
  Hwpart_t*   GetHwpart        ()           { return fHwpart; }

  // /HWPARP/
  Hwparp_t*   GetHwparp        ()           { return fHwparp; }

  // /HWBOSC/
  Hwbosc_t*   GetHwbosc        ()           { return fHwbosc; }

  // /HWPARC/
  Hwparc_t*   GetHwparc        ()           { return fHwparc; }

  // /HWBRCH/
  Hwbrch_t*   GetHwbrch        ()           { return fHwbrch; }

  // /HWEVNT/
  Hwevnt_t*   GetHwevnt        ()           { return fHwevnt; }
  double      GetAVWGT         () const     { return fHwevnt->AVWGT; }
  int         GetMAXPR         () const     { return fHwevnt->MAXPR; }
  void        SetMAXPR    (int i)           { fHwevnt->MAXPR = i; }

  void        SetEV1PR    (int i)           { fHwevnt->EV1PR = i; }
  void        SetEV2PR    (int i)           { fHwevnt->EV2PR = i; }

  int         GetMAXER         () const     { return fHwevnt->MAXER; }
  void        SetMAXER    (int i)           { fHwevnt->MAXER = i; }
  int         GetNRN      (int i) const     { return fHwevnt->NRN[i-1]; }
  void        SetNRN    (int i, int j)      { fHwevnt->NRN[i-1] = j; }
  double      GetEVWGT         () const     { return fHwevnt->EVWGT; }

  int         GetIDHW     (int i) const     { return fHwevnt->IDHW[i]; }

  int         GetIERROR        () const     { return fHwevnt->IERROR; }

  // /HWHARD/
  Hwhard_t*   GetHwhard        ()           { return fHwhard; }
  double      GetPTMIN         () const     { return fHwhard->PTMIN; }
  void        SetPTMIN    (double d)        { fHwhard->PTMIN = d; }
  double      GetPTPOW         () const     { return fHwhard->PTPOW; }
  void        SetPTPOW    (double d)        { fHwhard->PTPOW = d; }
  double      GetYJMIN         () const     { return fHwhard->YJMIN; }
  void        SetYJMIN    (double d)        { fHwhard->YJMIN = d; }
  double      GetYJMAX         () const     { return fHwhard->YJMAX; }
  void        SetYJMAX    (double d)        { fHwhard->YJMAX = d; }
  double      GetQ2MIN         () const     { return fHwhard->Q2MIN; }
  void        SetQ2MIN    (double d)        { fHwhard->Q2MIN = d; }
  double      GetQ2MAX         () const     { return fHwhard->Q2MAX; }
  void        SetQ2MAX    (double d)        { fHwhard->Q2MAX = d; }
  double      GetYBMIN         () const     { return fHwhard->YBMIN; }
  void        SetYBMIN    (double d)        { fHwhard->YBMIN = d; }
  double      GetYBMAX         () const     { return fHwhard->YBMAX; }
  void        SetYBMAX    (double d)        { fHwhard->YBMAX = d; }
  double      GetZJMAX        ()  const     { return fHwhard->ZJMAX; }
  void        SetZJMAX    (double d)        { fHwhard->ZJMAX = d; }

  // /HWPROP/
  Hwprop_t*   GetHwprop        ()           { return fHwprop; }
  double      GetRMASS      (int i) const   { return fHwprop->RMASS[i]; }
  void        SetRMASS    (int i, double r) { fHwprop->RMASS[i] = r; }

  // /HWUNAM/
  Hwunam_t*   GetHwunam        ()           { return fHwunam; }

  void        GetRNAME (int i, char a[9])   { for (int j=0;j<8;j++) a[j] = fHwunam->RNAME[i][j]; a[8] = '\0';}
/*  char*       GetRNAME(int i) { return fHwunam->RNAME[i]; }*/

  // /HWUPDT/
  Hwupdt_t*   GetHwupdt        ()           { return fHwupdt; }

  // /HWUWTS/
  Hwuwts_t*   GetHwuwts        ()           { return fHwuwts; }

  // /HWUCLU/
  Hwuclu_t*   GetHwuclu        ()           { return fHwuclu; }

  // /HWDIST/
  Hwdist_t*   GetHwdist        ()           { return fHwdist; }

  // /HWQDKT/
  Hwqdks_t*   GetHwqdkt        ()           { return fHwqdks; }

  // /HWUSUD/
  Hwusud_t*   GetHwusud        ()           { return fHwusud; }

  // /HWSUSY/
  Hwsusy_t*   GetHwsusy        ()           { return fHwsusy; }

  // /HWRPAR/
  Hwrpar_t*   GetHwrpar        ()           { return fHwrpar; }

  // /HWMINB/
  Hwminb_t*   GetHwminb        ()           { return fHwminb; }

  // /HWCLUS/
  Hwclus_t*   GetHwclus        ()           { return fHwclus; }

  // Herwig6 routines
  // the user would call
  //   Initialize
  //   change by himself the parameters s/he wants
  //   Hwusta to make stable the particles s/he wants
  //   PrepareRun
  //   GenerateEvent as many times as wished
  // An example is given in SetupTest

  void             GenerateEvent();
  void             Initialize(const char *beam, const char *target, double pbeam1, double pbeam2, int iproc);
  void             InitializeJimmy(const char *beam, const char *target, double pbeam1, double pbeam2, int iproc);
  void             PrepareRun();
  void             PrepareRunJimmy();
  void             OpenFortranFile(int lun, char* name);
  void             CloseFortranFile(int lun);
  Int_t            ImportParticles(TClonesArray *particles, Option_t *option="");
  TObjArray       *ImportParticles(Option_t *option="");
  TObjArray       *Particles() { return fParticles; }
  void             Hwigin();
  void             Hwuinc();
  void             Hwusta(char * name);
  void             Hweini();
  void             Hwuine();
  void             Hwepro();
  void             Hwbgen();
  void             Hwdhob();
  void             Hwcfor();
  void             Hwcdec();
  void             Hwdhad();
  void             Hwdhvy();
  void             Hwmevt();
  void             Hwufne();
  void             Hwefin();
  void             Hwiodk(int iopt);
  void             SetupTest();
 // Jimmy subroutines:
  void             Jminit();
  void             Jimmin();
  void             Jmefin();
protected:

  Hepevt_t* fHepevt; // Standard hep common block
  // Herwig6 common-blocks
  Hwbeam_t* fHwbeam; // Beams, process and number of events
  Hwbmch_t* fHwbmch; // Beams, process and number of events
  Hwproc_t* fHwproc; // Beams, process and number of events
  Hwpram_t* fHwpram; // Basic parameters (and quantities derived from them)
  Hwprch_t* fHwprch; // Basic parameters (and quantities derived from them)
  Hwpart_t* fHwpart; // Parton shower common
  Hwparp_t* fHwparp; // Parton polarization common
  Hwbosc_t* fHwbosc; // Electroweak boson common
  Hwparc_t* fHwparc; // Parton colour common
  Hwbrch_t* fHwbrch; // Branching common
  Hwevnt_t* fHwevnt; // Event common
  Hwhard_t* fHwhard; // Hard subprocess common
  Hwprop_t* fHwprop; // Particle properties
  Hwunam_t* fHwunam; // Particle properties
  Hwupdt_t* fHwupdt; // Particle decays
  Hwuwts_t* fHwuwts; // Weights used in cluster decays
  Hwuclu_t* fHwuclu; // Parameters for cluster decays
  Hwdist_t* fHwdist; // Variables controling mixing and vertex information
  Hwqdks_t* fHwqdks; // Arrays for temporarily storing heavy-b,c-hadrons decaying partonicaly
  Hwusud_t* fHwusud; // Parameters for Sudakov form factors
  Hwsusy_t* fHwsusy; // SUSY parameters
  Hwrpar_t* fHwrpar; // R-Parity violating parameters and colours
  Hwminb_t* fHwminb; // Parameters for minimum bias/soft underlying event
  Hwclus_t* fHwclus; // Cluster common used by soft event routines
  Hwgrav_t* fHwgrav;
  Hw6202_t* fHw6202;
  Hw6203_t* fHw6203;
  Hw6300_t* fHw6300;
  Hwpmrs_t* fHwpmrs;
  Hwcirc_t* fHwcirc;
  Hwdspb_t* fHwdspb;
  Hwdsp2_t* fHwdsp2;
  Hwdsp3_t* fHwdsp3;
  Hwdsp4_t* fHwdsp4;
  Hwdspn_t* fHwdspn;
  Hwspin_t* fHwspin;
  Hwstau_t* fHwstau;
  Hwgupr_t* fHwgupr;
  Hw6500_t* fHw6500;
  Hw6504_t* fHw6504;
  Hw6506_t* fHw6506;

  ClassDef(THerwig6,0)  //Interface to Herwig6.1 Event Generator
};

#endif
