#ifndef ROOT_ICommon
#define ROOT_ICommon

#ifndef __CFORTRAN_LOADED
//*KEEP,cfortran.
#include "cfortran.h"
//*KEND.
#endif

extern "C" {

/*========================================================================*/
/* COMMON/DYLIM/QMIN,QMAX,QTMIN,QTMAX,YWMIN,YWMAX,XWMIN,XWMAX,THWMIN,     */
/*     &        THWMAX,PHWMIN,PHWMAX,SETLMQ(12)                           */
/*------------------------------------------------------------------------*/
typedef struct {
    Float_t qmin;
    Float_t qmax;
    Float_t qtmin;
    Float_t qtmax;
    Float_t ywmin;
    Float_t ywmax;
    Float_t xwmin;
    Float_t xwmax;
    Float_t thwmin;
    Float_t thwmax;
    Float_t phwmin;
    Float_t phwmax;
    Bool_t setlmq[12];
    Bool_t ywset;
    Bool_t thwset;
} DylimCommon;
    
#define DYLIM COMMON_BLOCK(DYLIM, dylim)
COMMON_BLOCK_DEF(DylimCommon,DYLIM);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* qmin, qmax : Mass limits for W                                         */
/* qtmin, qtmax : q-subscript-t limits for W                              */
/* ywmin, ywmax : Y limits for W. Cannot be set simultaneously with       */
/*                thwmin, thwmax.                                         */
/* xwmin, xwmax : Feynman x limit for W                                   */
/* thwmin, thwmax : Theta limits for W. Cannot be set simultaneously      */
/*                  with ywmin, ywmax.                                    */
/* phwmin, phwmax : Phi limits for W                                      */
/* setlmq :                                                               */
/* ywset, thwset : Logical flags internal to the interface.               */
/**************************************************************************/


    
/*========================================================================*/
/* COMMON/EEPAR/SGMXEE,PLEP,PLEM,RSHMIN,RSHMAX,UPSLON,SIGZ,IBREM,IBEAM    */
/*------------------------------------------------------------------------*/
typedef struct {
    Float_t sgmxee;
    Float_t plep;
    Float_t plem;
    Float_t rshmin;
    Float_t rshmax;
    Float_t upslon;
    Float_t sigz;
    Bool_t ibrem;
    Bool_t ibeam;
} EeparCommon;
    
#define EEPAR COMMON_BLOCK(EEPAR, eepar)
COMMON_BLOCK_DEF(EeparCommon,EEPAR);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* sgmxee                                                                 */
/* plep : Polarisation of positron beam (P-subscript-L of e-superscript-+ */
/* plem : Polarisation of electron beam (P-subscript-L of e-superscript-- */
/*                                                                        */
/* Brem- and beam-strahlung parameters                                    */
/* rshmin : minimum subprocess energy (root-s-hyphen-subscript-min)       */
/* rshmax : maximum subprocess energy (root-s-hyphen-subscript-max)       */
/* upslon : beamstrahlung parameter (UPSILON)                             */
/* sigz : Longitudinal beam size in mm (sigma-subscript-z)                */
/* ibrem : True if EEBREM used                                            */
/* ibeam : True if EEBEAM used                                            */
/**************************************************************************/


    
/*========================================================================*/
/* COMMON/FORCE/NFORCE,IFORCE(MXFORC),MFORCE(5,MXFORC),LOOK2(2,MXFORC),   */
/*     &        LOOKST(MXFORC),MEFORC(MXFORC)                             */
/*------------------------------------------------------------------------*/
typedef struct {
    static const Int_t mxforc = 40;
    Int_t nforce;
    Int_t iforce[mxforc];
    Int_t mforce[mxforc][5];
    Int_t look2[mxforc][2];
    Int_t lookst[mxforc];
    Int_t meforc[mxforc];
} ForceCommon;
    
#define FORCE COMMON_BLOCK(FORCE, force)
COMMON_BLOCK_DEF(ForceCommon,FORCE);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* mxforce : Size of forced-decay array                                   */
/* nforce : Number of forced decay paths                                  */
/* iforce : ID code for particle of forced decay                          */
/* mforce : Mode(s) of forced decay                                       */
/* look2                                                                  */
/* lookst                                                                 */
/* meforc : Matrix element of the forced decay                            */
/**************************************************************************/ 


    
/*========================================================================*/
/* COMMON/FRGPAR/PUD,PBARY,SIGQT,PEND,XGEN(8),PSPIN1(8),PMIX1(3,2),       */
/*     &         PMIX2(3,2),XGENSS(9)                                     */
/*------------------------------------------------------------------------*/
typedef struct {
    Float_t pud;
    Float_t pbary;
    Float_t sigqt;
    Float_t pend;
    Float_t xgen[8];
    Float_t pspin1[8];
    Float_t pmix1[2][3];
    Float_t pmix2[2][3];
    Float_t xgenss[9];
    Float_t *pmixx1[6];
    Float_t *pmixx2[6];
    Float_t *frpar[32];
} FrgparCommon;
    
#define FRGPAR COMMON_BLOCK(FRGPAR, frgpar)
COMMON_BLOCK_DEF(FrgparCommon,FRGPAR);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* Fragmentation parameters :                                             */
/* pud :                                                                  */
/* pbary :                                                                */
/* sigqt : Internal k-subscript-t parameter for jet fragmentation (sigma) */
/* pend :                                                                 */
/* xgen : Jet fragmentation, Peterson with epsilon = a(n) / m, n = 4-8    */
/* pspin1 :                                                               */
/* pmix1 :                                                                */
/* pmix2 :                                                                */
/* xgenss : Fragmentation of GLSS, UPSS etc with epsilon = a(n)/m-squared */
/**************************************************************************/ 


    
/*========================================================================*/
/* COMMON/HCON/ANWWWW(4,4,4),ADWWWW(2,4),AIWWWW(4),HMASS,HGAM,HGAMS(29),  */
/*     &       ETAHGG,MATCHH(29),ZSTARS(4,2),IHTYPE,HGAMSS(85,85)         */
/*------------------------------------------------------------------------*/
typedef struct {
    Float_t anwwww[4][4][4];
    Float_t adwwww[4][2];
    Float_t aiwwww[4];
    Float_t hmass;
    Float_t hgam;
    Float_t hgams[29];
    Float_t etahgg;
    Int_t matchh[29];
    Float_t zstars[2][4];
    Int_t ihtype;
    Float_t hgamss[85][85];
} HconCommon;
    
#define HCON COMMON_BLOCK(HCON, hcon)
COMMON_BLOCK_DEF(HconCommon,HCON);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* anwwww :                                                               */
/* adwwww :                                                               */
/* aiwwww :                                                               */
/* hmass :                                                                */
/* hgam :                                                                 */
/* hgams :                                                                */
/* etahgg:                                                                */
/* matchh :                                                               */
/* zstars :                                                               */
/* ihtype : MSSM Higgs type, either HL0, HH0 or HA0                       */
/* hgamss :                                                               */
/**************************************************************************/


    
/*========================================================================*/
/* COMMON/JETLIM/PMIN(MXLIM),PMAX(MXLIM),PTMIN(MXLIM),PTMAX(MXLIM),       */
/*     &         YJMIN(MXLIM),YJMAX(MXLIM),PHIMIN(MXLIM),PHIMAX(MXLIM),   */
/*     &         XJMIN(MXLIM),XJMAX(MXLIM),THMIN(MXLIM),THMAX(MXLIM),     */
/*     &         SETLMJ(12*MXLIM)                                         */
/*------------------------------------------------------------------------*/
typedef struct {
    static const Int_t mxlim = 8;
    Float_t pmin[mxlim];
    Float_t pmax[mxlim];
    Float_t ptmin[mxlim];
    Float_t ptmax[mxlim];
    Float_t yjmin[mxlim];
    Float_t yjmax[mxlim];    
    Float_t phimin[mxlim];
    Float_t phimax[mxlim];
    Float_t xjmin[mxlim];
    Float_t xjmax[mxlim];
    Float_t thmin[mxlim];
    Float_t thmax[mxlim];
    Bool_t setlmj[12*mxlim];
} JetlimCommon;
    
#define JETLIM COMMON_BLOCK(JETLIM, jetlim)
COMMON_BLOCK_DEF(JetlimCommon,JETLIM);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* Changes to mxlim should be echoed in MGLIMS and some of its access     */
/* routines.                                                              */
/* pmin, pmax : Momentum limits for jets.                                 */  
/* ptmin, ptmax : p-subscript-t limits for jets.                          */
/* yjmin, yjmax : Y limits for each jet; not simultaneous to thmin, thmax.*/
/* phimin, phimax : Phi limits for jets.                                  */
/* xjmin, xjmax : Feynman x limits for jets.                              */
/* thmin, thmax : Theta limits for jets, not simultaneous with yjmin,     */
/*                yjmax.                                                  */
/* setlmj :                                                               */
/* yjset, thset : Logical flags internal to the interface.                 */
/**************************************************************************/



/*========================================================================*/
/* COMMON/JETPAR/P(3),PT(3),YJ(3),PHI(3),XJ(3),TH(3),CTH(3),STH(3),       */
/*     &         JETTYP(3),SHAT,THAT,QSQ,X1,X2,PBEAM(2),                  */
/*     &         QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,JWTYP,  */
/*     &         INITYP(2),ISIGS,PBEAMS(5)                                */
/*------------------------------------------------------------------------*/
typedef struct {
    Float_t p[3];
    Float_t pt[3];
    Float_t yj[3];
    Float_t phi[3];
    Float_t xj[3];
    Float_t th[3];
    Float_t cth[3];
    Float_t sth[3];
    Int_t jettyp[3];
    Float_t shat;
    Float_t that;
    Float_t uhat;
    Float_t qsq;
    Float_t x1;
    Float_t x2;
    Float_t pbeam[2];
    Float_t qmw;
    Float_t qw;
    Float_t qtw;
    Float_t yw;
    Float_t xw;
    Float_t thw;
    Float_t qtmw;
    Float_t phiw;
    Float_t shat1;
    Float_t that1;
    Float_t uhat1;
    Int_t jwtyp;
    Float_t alfqsq;
    Float_t cthw;
    Float_t sthw;
    Float_t q0w;
    Int_t inityp[2];
    Int_t isigs;
    Float_t pbeams[5];
} JetparCommon;
    
#define JETPAR COMMON_BLOCK(JETPAR, jetpar)
COMMON_BLOCK_DEF(JetparCommon,JETPAR);

    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* jwtyp : W type, either W+-, Z0 or GM.                                  */
/**************************************************************************/
    


/*========================================================================*/
/* COMMON/KEYS/IKEY, KEYON, KEYS(MXKEYS)                                 */
/*------------------------------------------------------------------------*/

typedef struct {
    static const Int_t mxkeys = 20;
    Bool_t keys[mxkeys];
    Bool_t keyon;
    Int_t ikey;
    Char_t *reac;
} KeysCommon;
    
#define KEYS COMMON_BLOCK(KEYS, keys)
COMMON_BLOCK_DEF(KeysCommon,KEYS);

/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* keys : Flag for process type.                                          */
/**************************************************************************/ 


    
/*========================================================================*/
/* COMMON/KKGRAV/NEXTRAD,MASSD,KKGSD,SURFD,UVCUT                          */
/*------------------------------------------------------------------------*/

typedef struct {
    Int_t nextrad;
    Float_t massd;
    Float_t kkgsd;
    Float_t surfd;
    Bool_t uvcut;
} KkgravCommon;
    
#define KKGRAV COMMON_BLOCK(KKGRAV, kkgrav)
COMMON_BLOCK_DEF(KkgravCommon,KKGRAV);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* EXTRADIM process parameters :                                          */
/* nextrad : delta                                                        */
/* massd : M-subscript-D                                                  */
/* uvcut : Logical flag                                                   */
/*                                                                        */   
/* kkgsd :                                                                */
/* surfd :                                                                */
/**************************************************************************/ 



/*========================================================================*/
/* COMMON/MBGEN/POMWT(LIMPOM),POMGEN(LIMPOM),MNPOM,MXPOM,PDIFFR,NPOM,     */
/*     &        XBARY(2),DXBARY(2),XPOM(LIMPOM,2)                         */
/*------------------------------------------------------------------------*/
typedef struct {
    static const Int_t limpom = 20;
    Float_t pomwt[limpom];
    Float_t pomgen[limpom];
    Int_t mnpom;
    Int_t mxpom;
    Float_t pdiffr;
    Int_t npom;
    Float_t xbary[2];
    Float_t dxbary[2];
    Float_t xpom[2][limpom];
} MbgenCommon;
    
#define MBGEN COMMON_BLOCK(MBGEN, mbgen)
COMMON_BLOCK_DEF(MbgenCommon,MBGEN);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* pomwt                                                                  */
/* pomgen                                                                 */  
/* mnpom, mxpom : Min and max number of cut pomerons.                     */
/* pdiffr                                                                 */
/* npom                                                                   */
/* xbary                                                                  */
/* dxbary                                                                 */
/* xpom                                                                   */
/**************************************************************************/


    
/*========================================================================*/
/* COMMON/MGLIMS/EHMGMN,EHMGMX,YHMGMN,YHMGMX,AMIJMN(MXLIM,MXLIM),         */
/*     &         AMIJMX(MXLIM,MXLIM),FIXMIJ(MXLIM,MXLIM)                  */
/*------------------------------------------------------------------------*/
typedef struct {
    static const Int_t mxlim = 8;
    Float_t ehmgmn;
    Float_t ehmgmx;
    Float_t yhmgmn;
    Float_t yhmgmx;
    Float_t amijmn[mxlim][mxlim];
    Float_t amijmx[mxlim][mxlim];
    Bool_t fixmij[mxlim][mxlim];
} MglimsCommon;
    
#define MGLIMS COMMON_BLOCK(MGLIMS, mglims)
COMMON_BLOCK_DEF(MglimsCommon,MGLIMS);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* Changes in mxlim should be echoed in JETLIM and in access routines for */
/* amijmn and amijmx.                                                     */
/* Limits for MadGraph multiparton processes                              */
/* ehmgmn, ehmgmx : Mass range                                            */
/* yhmgmn, yhmgmx :                                                       */
/* amijmn, amijmx : Multimet mass limits                                  */
/* fixmij :                                                               */
/**************************************************************************/ 


    
/*========================================================================*/
/* COMMON/NODCAY/NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOGRAV           */
/*------------------------------------------------------------------------*/
typedef struct {
    Bool_t nodcay;
    Bool_t noeta;
    Bool_t nopi0;
    Bool_t nonunu;
    Bool_t noevol;
    Bool_t nohadr;
    Bool_t nograv;
} NodcayCommon;
    
#define NODCAY COMMON_BLOCK(NODCAY, nodcay)
COMMON_BLOCK_DEF(NodcayCommon,NODCAY);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* Booleans used to forbid decay channels.                                */
/* nodcay : Suppress all decays                                           */
/* noeta  : Suppress eta decays                                           */
/* nopi0  : Suppress pi-zero decays                                       */
/* nonunu : Suppess Z-zero neutrino decays                                */
/* noevol : Suppress QCD evolution and hadronisation                      */
/* nohadr : Suppress hadronisation of jets and beam jets                  */
/* nograv : Suppress gravitino decays in GMSB model                       */
/**************************************************************************/


    
/*========================================================================*/
/* COMMON/PARTCL/NPTCL,PPTCL(5, MXPTCL),IORIG(MXPTCL),IDENT(MXPTCL),      */
/*     &         IDCAY(MXPTCL)                                            */
/*------------------------------------------------------------------------*/
typedef struct {
    static const Int_t mxptcl = 4000;
    static const Int_t ipack = 10000;
    Int_t nptcl;
    Float_t pptcl[mxptcl][5];
    Int_t iorig[mxptcl];
    Int_t ident[mxptcl];
    Int_t idcay[mxptcl];
} PartclCommon;
    
#define PARTCL COMMON_BLOCK(PARTCL, partcl)
COMMON_BLOCK_DEF(PartclCommon,PARTCL);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* Stores HEPEVT data.                                                    */
/* nptcl : Number of particles.                                           */
/* pptcl : x, y, z, and 0 momentum of particles, and mass.                */
/* iorig : origin of particles.                                           */
/* ident : ID code.                                                       */
/* idcay : ID of decay products.                                          */
/* mxptcl : Max number of particles.                                      */
/* ipack : Packing integer, used in decoding idcay.                       */
/**************************************************************************/


    
/*========================================================================*/
/* COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA          */
/*------------------------------------------------------------------------*/
typedef struct {
    Int_t njet;
    Float_t scm;
    Float_t halfe;
    Float_t ecm;
    Int_t idin[2];
    Int_t nevent;
    Int_t ntries;
    Int_t nsigma;
} PrimarCommon;
    
#define PRIMAR COMMON_BLOCK(PRIMAR, primar)
COMMON_BLOCK_DEF(PrimarCommon,PRIMAR);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* njet :                                                                 */
/* scm : Square of the center-of-mass energy.                             */
/* halfe : Half the center-of-mass energy.                                */
/* ecm : Center-of-mass energy.                                           */
/* idin : Beam types 1 and 2.                                             */
/* nevent : Number of events to build.                                    */
/* ntries : Maximum number of tries to find a good event.                 */
/* nsigma : Number of unevolved events for SIGF calculation.              */
/**************************************************************************/



/*========================================================================*/
/* COMMON/PRTOUT/NEVPRT,NJUMP                                             */    
/*------------------------------------------------------------------------*/
typedef struct {
    Int_t nevprt;
    Int_t njump;    
} PrtoutCommon;
    
#define PRTOUT COMMON_BLOCK(PRTOUT, prtout)
COMMON_BLOCK_DEF(PrtoutCommon,PRTOUT);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* nevprt : Number of events to print                                     */
/* njump : Print every njump events.                                      */
/**************************************************************************/
    

    
/*========================================================================*/
/* COMMON/QCDPAR/ALAM,ALAM2,CUTJET,ISTRUC                	          */
/*------------------------------------------------------------------------*/
typedef struct {
    Float_t alam;
    Float_t alam2;
    Float_t cutjet;
    Int_t istruc;
} QcdparCommon;
    
#define QCDPAR COMMON_BLOCK(QCDPAR, qcdpar)
COMMON_BLOCK_DEF(QcdparCommon,QCDPAR);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* alam : QCD scale (LAMBDA)                                              */
/* alam2 : Square of alam                                                 */
/* cutjet : Cutoff mass for QCD jet evolution (mu-subscript-c)            */
/* istruc : Structure functions CTEQ3L, CTEQ2L, EHLQ or DO.               */   
/**************************************************************************/

    

/*========================================================================*/
/* COMMON/QLMASS/AMLEP(100),NQLEP,NMES,NBARY               	          */
/*------------------------------------------------------------------------*/
typedef struct {
    static const Int_t nqlep = 61;
    static const Int_t nmes = 2;
    static const Int_t nbary = 2;
    Float_t amlep[100];
} QlmassCommon;
    
#define QLMASS COMMON_BLOCK(QLMASS, qlmass)
COMMON_BLOCK_DEF(QlmassCommon,QLMASS);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* amlep : (C++ indices, Fortran is one greater)                          */
/* [5-7] : t, y and x quark mass (m-subscript-t, -y, -x)                  */
/* [21-26] : Masses of u, d, s, c, b and t, all -tilde                    */
/* [29] : g-tilde mass                                                    */
/* [30] : gamma-tilde mass                                                */
/* [31-36] : Masses for nu-e, e, nu-mu, mu, nu-tau and tau, all -tilde    */
/* [39] : W-plus-tilde mass                                               */
/* [40] : Z-zero-tilde mass                                               */
/* [63-71] : Higgs meson masses, charges 0,0,0,0,0,1,1,2,2                */
/**************************************************************************/


    
/*========================================================================*/
/* COMMON/SEED/XSEED                                      	          */
/*------------------------------------------------------------------------*/
typedef struct {
    Char_t xseed[24];
} SeedCommon;
    
#define SEED COMMON_BLOCK(SEED, seed)
COMMON_BLOCK_DEF(SeedCommon,SEED);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* Random number seed                                                     */
/**************************************************************************/


    
/*========================================================================*/
/* COMMON/SUGNU/XNUSUG(18)                                	          */
/*------------------------------------------------------------------------*/
typedef struct {
    Float_t xnusug[18];
} SugnuCommon;
    
#define SUGNU COMMON_BLOCK(SUGNU, sugnu)
COMMON_BLOCK_DEF(SugnuCommon,SUGNU);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* Non-universal SUGRA terms : (C++ indexes again)                        */
/* xnusug[0-2] : Gaugino masses                                           */
/* [3-5] : A terms (A-subscript-tau, -b, -t)                              */
/* [6-7] : Higgs masses (M-subscript-H-subscript-d, -u)                   */
/* [8-12] : 1st / 2nd generation masses (M-subscript-eR, -eL, -dR,        */
/*                                                  -uR, -uL)             */
/* [13-17] : 3rd generation masses (M-subscript-tau R, -tau L, -bR,       */
/*                                             -tR, -tL)                  */
/**************************************************************************/


    
/*========================================================================*/
/* COMMON/TCPAR/TCMRHO,TCGRHO                             	          */
/*------------------------------------------------------------------------*/
typedef struct {
    Float_t tcmrho;
    Float_t tcgrho;
} TcparCommon;
    
#define TCPAR COMMON_BLOCK(TCPAR, tcpar)
COMMON_BLOCK_DEF(TcparCommon,TCPAR);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* tcmrho : Technicolor mass                                              */
/* tcgrho : Technicolor width                                             */
/**************************************************************************/


    
/*========================================================================*/
/* COMMON/TYPES/LOC(100),NTYP,NJTTYP,NWWTYP(2),NWMODE(3)                  */
/*------------------------------------------------------------------------*/
typedef struct {
    static const Int_t mxtype = 8;
    Int_t loc[100];
    Int_t ntyp;
    Int_t njttyp[mxtype];
    Int_t nwwtyp[2];
    Int_t nwmode[3];
} TypesCommon;
    
#define TYPES COMMON_BLOCK(TYPES, types)
COMMON_BLOCK_DEF(TypesCommon,TYPES);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* Note : Any change in mxtype should be copied in xtypes.                */
/* loc :                                                                  */
/* ntyp :                                                                 */
/* njttyp : Number of jet types in jet [index].                           */
/* nwwtyp :                                                               */
/* nwmode :                                                               */
/**************************************************************************/



/*========================================================================*/
/* COMMON/XMSSM/GOMSSM,GOSUG,GOGMSB,GOAMSB,AL3UNI,                        */
/*     &        XGLSS,XMUSS,XHASS,XTBSS,                                  */
/*     &        XQ1SS,XDRSS,XURSS,XL1SS,XERSS,                            */
/*     &        XQ2SS,XSRSS,XCRSS,XL2SS,XMRSS,                            */
/*     &        XQ3SS,XBRSS,XTRSS,XL3SS,XTARSS,XATSS,XABSS,XATASS,        */
/*     &        XM1SS,XM2SS,XM0SU,XMHSU,XA0SU,XTGBSU,XSMUSU,              */
/*     &        XLAMGM,XMESGM,XN5GM,XCMGV,XMGVTO,                         */
/*     &        XRSLGM,XDHDGM,XDHUGM,XDYGM,XN51GM,XN52GM,XN53GM,          */
/*     &        XMN3NR,XMAJNR,XANSS,XNRSS,XSBCS                           */
/*------------------------------------------------------------------------*/
    
typedef struct {
    Bool_t gomssm;
    Bool_t gosug;
    Bool_t gogmsb;
    Bool_t goamsb;
    Bool_t al3uni;
    Float_t xglss;
    Float_t xmuss;
    Float_t xhass;
    Float_t xtbss;
    Float_t xq1ss;
    Float_t xdrss;
    Float_t xurss;
    Float_t xl1ss;
    Float_t xerss;
    Float_t xq2ss;
    Float_t xsrss;
    Float_t xcrss;
    Float_t xl2ss;
    Float_t xmrss;
    Float_t xq3ss;
    Float_t xbrss;
    Float_t xtrss;
    Float_t xl3ss;
    Float_t xtarss;
    Float_t xatss;
    Float_t xabss;
    Float_t xatass;
    Float_t xm1ss;
    Float_t xm2ss;
    Float_t xm0su;
    Float_t xmhsu;
    Float_t xa0su;
    Float_t xtgbsu;
    Float_t xsmusu;
    Float_t xlamgm;
    Float_t xmesgm;
    Float_t xn5gm;
    Float_t xcmgv;
    Float_t mgvto;
    Float_t xrslgm;
    Float_t xdhdgm;
    Float_t xdhugm;
    Float_t xdygm;
    Float_t xn51gm;
    Float_t xn52gm;
    Float_t xn53gm;
    Float_t xmn3nr;
    Float_t xmajnr;
    Float_t xanss;
    Float_t xnrss;
    Float_t xsbcs;

} XmssmCommon;
    
#define XMSSM COMMON_BLOCK(XMSSM, xmssm)
COMMON_BLOCK_DEF(XmssmCommon,XMSSM);

/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* Keyword flags :                                                        */
/* gomssm : True if any of the MSSM* keywords or SUGRA, GMSB or AMSB have */
/*          been used.                                                    */
/* gosug  : True if SUGRA has been used.                                  */
/* gogmsb : True if GMSB has been used.                                   */
/* goamsb : True if AMSB has been used.                                   */
/* al3uni : True if AL3UNI has been used.                                 */
/*                                                                        */
/* MSSM parameters :                                                      */
/* xglss : Gluino mass                                                    */
/* xmuss : mu                                                             */
/* xhass : A mass                                                         */
/* xtbss : tan beta                                                       */
/*                                                                        */
/* MSSM first generation :                                                */
/* xq1ss : q-subscript-1 mass                                             */
/* xdrss : d-subscript-r mass                                             */
/* xurss : u-subscript-r mass                                             */
/* xl1ss : l-subscript-1 mass                                             */
/* xerss : e-subscript-1 mass                                             */
/*                                                                        */
/* MSSM second generation :                                               */
/* xq2ss : q-subscript-2 mass                                             */
/* xsrss : s-subscript-r mass                                             */
/* xcrss : c-subscript-r mass                                             */
/* xl2ss : l-subscript-2 mass                                             */
/* xmrss : mu-subscript-r mass                                            */
/*                                                                        */
/* MSSM third generation :                                                */
/* xq3ss : q-subscript-3 mass                                             */
/* xbrss : b-subscript-r mass                                             */
/* xtrss : t-subscript-r mass                                             */
/* xl3ss : l-subscript-3 mass                                             */
/* xtarss : tau-subscript-r mass                                          */
/* xatss, xabss, xatass : A-subscript-t, -b, -tau mixings.                */
/*                                                                        */
/* MSSM gaugino masses                                                    */
/* xm1ss, xm2ss : M-subscript-1, -2                                       */
/*                                                                        */
/* Anomaly-mediated SUSY breaking / Minimal supergravity parameters :     */
/* xm0su : scalar mass (m-subscript-0)                                    */
/* xmhsu : gravitino mass (m-subscript-three halves)                      */
/* xa0su : trilinear breaking term (A-subscript-0)                        */
/* xtgbsu : VEV ratio (tan beta)                                          */
/* xsmusu : sign (sgn mu)                                                 */
/*                                                                        */
/* GMSB messinger SUSY breaking                                           */
/* xlamgm : mass (LAMBDA-subscript-m)                                     */
/* xmesgm :                                                               */
/* xn5gm : number of 5 + 5-bar (M-subscript-m)                            */
/* xcmgv : gravitino scale (C-subscript-gr)                               */
/*                                                                        */
/* Non-minimal GMSB parameters                                            */
/* xrslgm : gaugino mass multiplier (R-slash)                             */
/* xdhdgm : Higgs mass shift (delta M-squared-subscript-Hd)               */
/* xdhugm : Higgs mass shift (delta M-squared-subscript-Hu)               */
/* xdygm : D-term mass-squared (D-subscript-Y of M)                       */
/* xn51gm, xn52gm, xn53gm : independent gauge group messengers            */
/*                          (N-subscript-5-subscript-1, -2, -3)           */
/* SUGRA see-saw nu-effect                                                */
/* xmn3nr : nu-mass (m-subscript-nu-subscript-tau)                        */
/* xmajnr : int. scale (M-subscript-N)                                    */
/* xanss : GUT scale (A-subscript-n)                                      */
/* xnrss : nu SSB terms (m-subscript-nu-tilde-subscript-R)                */
/*                                                                        */
/* xsbcs                                                                  */
/*                                                                        */
/* Gravitino mass                                                         */
/* xmgvto : gravitino mass (M-subscript-gravitino)                        */
/**************************************************************************/  


    
/*========================================================================*/
/* COMMON/XTYPES/PARTYP(40),TITLE(10),JETYP(30,MXTYPE),WWTYP(30,2),       */
/*     &         WMODES(30,3)                                             */
/*------------------------------------------------------------------------*/
typedef struct {
    static const Int_t mxtype = 8;
    Char_t* partyp[40];
    Char_t* title;
    Char_t* jetyp[mxtype][30];
    Char_t* wwtyp[2][30];
    Char_t* wmodes[3][30];
} XtypesCommon;
    
#define XTYPES COMMON_BLOCK(XTYPES, xtypes)
COMMON_BLOCK_DEF(XtypesCommon,XTYPES);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/   
/* Note : Any change in mxtype should be copied in types.                 */
/* partyp                                                                 */
/* title                                                                  */
/* jetyp                                                                  */
/* wwtyp : Decay modes for WPAIR                                          */
/* wmodes : Decay modes for weak force carriers                           */
/**************************************************************************/



/*========================================================================*/
/* COMMON/WCON/SIN2W,WMASS(4),WGAM(4),AQ(12,4),BQ(12,4),COUT(4),          */
/*     &       MATCH(25,4),WCBR(25,4),CUTOFF,CUTPOW,TBRWW(4,2),           */
/*     &       RBRWW(12,4,2),EZ,AQDP(12,4),BQDP(12,4),EZDP,WFUDGE         */
/*------------------------------------------------------------------------*/
typedef struct {
    
    Float_t sin2w;
    Float_t wmass[4];
    Float_t wgam[4];
    Float_t aq[4][12];
    Float_t bq[4][12];
    Float_t cout[4];
    Int_t match;
    Float_t wcbr[4][25];
    Float_t cutoff;
    Float_t cutpow;
    Float_t tbrww[2][4];
    Float_t rbrww[2][4][12];
    Float_t ez;
    Float_t aqdp[4][12];
    Float_t bqdp[4][12];
    Float_t ezdp;
    Float_t wfudge;
} WconCommon;
    
#define WCON COMMON_BLOCK(WCON, wcon)
COMMON_BLOCK_DEF(WconCommon, WCON);
    
/**************************************************************************/
/*           D E S C R I P T I O N :                                      */
/*------------------------------------------------------------------------*/
/* sin2w : Weinberg angle (sin-squared(theta-subscript-W))                */
/* wmass : W and Z masses                                                 */
/* wgam  :                                                                */
/* aq    :                                                                */
/* bq    :                                                                */
/* cout  :                                                                */
/* match                                                                  */
/* wcbr                                                                   */
/* cutoff, cutpow : mu-square and nu, respectively, in formula            */
/*                  cutoff q*t-square = mu-square*Q-power-nu              */
/* 		    for DRELLYAN events.                                  */
/* tbrww                                                                  */
/* rbrww                                                                  */
/* ez                                                                     */
/* aqdp                                                                   */
/* bqdp                                                                   */
/* ezdp                                                                   */
/* wfudge : Fudge factor for DRELLYAN evolution scale                     */
/**************************************************************************/
}

#endif
    
// Endfile.








