#ifndef TGeant3_H 
#define TGeant3_H 
//////////////////////////////////////////////// 
//  C++ interface to Geant3 basic routines    // 
//////////////////////////////////////////////// 
 
#include <AliMC.h> 
  
//______________________________________________________________
//
//       Geant3 prototypes for commons
//
//______________________________________________________________
//

//----------QUEST 
//      COMMON/QUEST/IQUEST(100) 
typedef struct { 
  Int_t    iquest[100]; 
} Quest_t; 
 
//----------GCBANK
//      COMMON/GCBANK/NZEBRA,GVERSN,ZVERSN,IXSTOR,IXDIV,IXCONS,FENDQ(16)
//     +             ,LMAIN,LR1,WS(KWBANK)
typedef struct {
  Int_t nzebra;
  Float_t gversn;
  Float_t zversn;
  Int_t ixstor;
  Int_t ixdiv;
  Int_t ixcons;
  Float_t fendq[16];
  Int_t lmain;
  Int_t lr1;
} Gcbank_t;

//----------GCLINK 
//      COMMON/GCLINK/JDIGI ,JDRAW ,JHEAD ,JHITS ,JKINE ,JMATE ,JPART 
//     +      ,JROTM ,JRUNG ,JSET  ,JSTAK ,JGSTAT,JTMED ,JTRACK,JVERTX 
//     +      ,JVOLUM,JXYZ  ,JGPAR ,JGPAR2,JSKLT 
typedef struct { 
  Int_t    jdigi; 
  Int_t    jdraw; 
  Int_t    jhead; 
  Int_t    jhits; 
  Int_t    jkine; 
  Int_t    jmate; 
  Int_t    jpart; 
  Int_t    jrotm; 
  Int_t    jrung; 
  Int_t    jset; 
  Int_t    jstak; 
  Int_t    jgstat; 
  Int_t    jtmed; 
  Int_t    jtrack; 
  Int_t    jvertx; 
  Int_t    jvolum; 
  Int_t    jxyz; 
  Int_t    jgpar; 
  Int_t    jgpar2; 
  Int_t    jsklt; 
} Gclink_t; 
 
 
//----------GCFLAG 
//      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN 
//     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2) 
typedef struct { 
  Int_t    idebug; 
  Int_t    idemin; 
  Int_t    idemax; 
  Int_t    itest; 
  Int_t    idrun; 
  Int_t    idevt; 
  Int_t    ieorun; 
  Int_t    ieotri; 
  Int_t    ievent; 
  Int_t    iswit[10]; 
  Int_t    ifinit[20]; 
  Int_t    nevent; 
  Int_t    nrndm[2]; 
} Gcflag_t; 
 
//----------GCKINE 
//      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP 
//     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD 
typedef struct { 
  Int_t    ikine; 
  Float_t  pkine[10]; 
  Int_t    itra; 
  Int_t    istak; 
  Int_t    ivert; 
  Int_t    ipart; 
  Int_t    itrtyp; 
  Int_t    napart[5]; 
  Float_t  amass; 
  Float_t  charge; 
  Float_t  tlife; 
  Float_t  vert[3]; 
  Float_t  pvert[4]; 
  Int_t    ipaold; 
} Gckine_t; 
 
//----------GCKING 
//      COMMON/GCKING/KCASE,NGKINE,GKIN(5,MXGKIN), 
//     +                           TOFD(MXGKIN),IFLGK(MXGKIN) 
#define MXGKIN 100 
typedef struct  { 
  Int_t    kcase; 
  Int_t    ngkine; 
  Float_t  gkin[MXGKIN][5]; 
  Int_t    tofd[MXGKIN]; 
  Int_t    iflgk[MXGKIN]; 
} Gcking_t; 

//----------GCKIN2
//      COMMON/GCKIN2/NGPHOT,XPHOT(11,MXPHOT)
#define MXPHOT 800
typedef struct {
  Int_t ngphot;
  Float_t xphot[MXPHOT][11];
} Gckin2_t;

//----------GCKIN3 
//      COMMON/GCKIN3/GPOS(3,MXGKIN)
typedef struct {
  Float_t gpos[MXGKIN][3];
} Gckin3_t;

//----------GCMATE 
//      COMMON/GCMATE/NMAT,NAMATE(5),A,Z,DENS,RADL,ABSL 
typedef struct { 
  Int_t    nmat; 
  Int_t    namate[5]; 
  Float_t  a; 
  Float_t  z; 
  Float_t  dens; 
  Float_t  radl; 
  Float_t  absl; 
} Gcmate_t; 
 
//----------GCTMED 
//      COMMON/GCTMED/NUMED,NATMED(5),ISVOL,IFIELD,FIELDM,TMAXFD,STEMAX 
//     +      ,DEEMAX,EPSIL,STMIN,CFIELD,PREC,IUPD,ISTPAR,NUMOLD 
typedef struct { 
  Int_t    numed; 
  Int_t    natmed[5]; 
  Int_t    isvol; 
  Int_t    ifield; 
  Float_t  fieldm; 
  Float_t  tmaxfd; 
  Float_t  stemax; 
  Float_t  deemax; 
  Float_t  epsil; 
  Float_t  stmin; 
  Float_t  cfield; 
  Float_t  prec; 
  Int_t    iupd; 
  Int_t    istpar; 
  Int_t    numold; 
} Gctmed_t; 
 
//----------GCTRAK 
#define MAXMEC 30 
//      PARAMETER (MAXMEC=30) 
//      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC) 
//     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG 
//     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL 
//     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN 
//     + ,NLVSAV,ISTORY 
typedef struct { 
  Float_t  vect[7]; 
  Float_t  getot; 
  Float_t  gekin; 
  Int_t    vout[7]; 
  Int_t    nmec; 
  Int_t    lmec[MAXMEC]; 
  Int_t    namec[MAXMEC]; 
  Int_t    nstep; 
  Int_t    maxnst; 
  Float_t  destep; 
  Float_t  destel; 
  Float_t  safety; 
  Float_t  sleng; 
  Float_t  step; 
  Float_t  snext; 
  Float_t  sfield; 
  Float_t  tofg; 
  Float_t  gekrat; 
  Float_t  upwght; 
  Int_t    ignext; 
  Int_t    inwvol; 
  Int_t    istop; 
  Int_t    igauto; 
  Int_t    iekbin; 
  Int_t    ilosl; 
  Int_t    imull; 
  Int_t    ingoto; 
  Int_t    nldown; 
  Int_t    nlevin; 
  Int_t    nlsav; 
  Int_t    istory; 
} Gctrak_t; 
 
//----------GCVOLU 
//      COMMON/GCVOLU/NLEVEL,NAMES(15),NUMBER(15), 
//     +LVOLUM(15),LINDEX(15),INFROM,NLEVMX,NLDEV(15),LINMX(15), 
//     +GTRAN(3,15),GRMAT(10,15),GONLY(15),GLX(3) 
typedef struct { 
  Int_t    nlevel; 
  Int_t    names[15]; 
  Int_t    number[15]; 
  Int_t    lvolum[15]; 
  Int_t    lindex[15]; 
  Int_t    infrom; 
  Int_t    nlevmx; 
  Int_t    nldev[15]; 
  Int_t    linmx[15]; 
  Float_t  gtran[15][3]; 
  Float_t  grmat[15][10]; 
  Float_t  gonly[15]; 
  Float_t  glx[3]; 
} Gcvolu_t; 
 
//----------GCSETS 
//  COMMON/GCSETS/IHSET,IHDET,ISET,IDET,IDTYPE,NVNAME,NUMBV(20) 
typedef struct { 
  Int_t    ihset; 
  Int_t    ihdet; 
  Int_t    iset; 
  Int_t    idet; 
  Int_t    idtype; 
  Int_t    nvname; 
  Int_t    numbv[20]; 
} Gcsets_t; 
 
//----------GCNUM 
//   COMMON/GCNUM/NMATE ,NVOLUM,NROTM,NTMED,NTMULT,NTRACK,NPART 
//  +            ,NSTMAX,NVERTX,NHEAD,NBIT 
typedef struct { 
  Int_t    nmate; 
  Int_t    nvolum; 
  Int_t    nrotm; 
  Int_t    ntmed; 
  Int_t    ntmult; 
  Int_t    ntrack; 
  Int_t    npart; 
  Int_t    nstmax; 
  Int_t    nvertx; 
  Int_t    nhead; 
  Int_t    nbit; 
} Gcnum_t; 
 
//----------GCCUTS 
//  COMMON/GCCUTS/CUTGAM,CUTELE,CUTNEU,CUTHAD,CUTMUO,BCUTE,BCUTM 
//   +             ,DCUTE ,DCUTM ,PPCUTM,TOFMAX,GCUTS(5) 
typedef struct { 
  Float_t cutgam; 
  Float_t cutele; 
  Float_t cutneu; 
  Float_t cuthad; 
  Float_t cutmuo; 
  Float_t bcute; 
  Float_t bcutm; 
  Float_t dcute; 
  Float_t dcutm; 
  Float_t ppcutm; 
  Float_t tofmax; 
  Float_t gcuts[5]; 
} Gccuts_t; 

//----------GCPHYS
//      COMMON/GCPHYS/IPAIR,SPAIR,SLPAIR,ZINTPA,STEPPA
//     +             ,ICOMP,SCOMP,SLCOMP,ZINTCO,STEPCO
//     +             ,IPHOT,SPHOT,SLPHOT,ZINTPH,STEPPH
//     +             ,IPFIS,SPFIS,SLPFIS,ZINTPF,STEPPF
//     +             ,IDRAY,SDRAY,SLDRAY,ZINTDR,STEPDR
//     +             ,IANNI,SANNI,SLANNI,ZINTAN,STEPAN
//     +             ,IBREM,SBREM,SLBREM,ZINTBR,STEPBR
//     +             ,IHADR,SHADR,SLHADR,ZINTHA,STEPHA
//     +             ,IMUNU,SMUNU,SLMUNU,ZINTMU,STEPMU
//     +             ,IDCAY,SDCAY,SLIFE ,SUMLIF,DPHYS1
//     +             ,ILOSS,SLOSS,SOLOSS,STLOSS,DPHYS2
//     +             ,IMULS,SMULS,SOMULS,STMULS,DPHYS3
//     +             ,IRAYL,SRAYL,SLRAYL,ZINTRA,STEPRA
//      COMMON/GCPHLT/ILABS,SLABS,SLLABS,ZINTLA,STEPLA
//     +             ,ISYNC
//     +             ,ISTRA
typedef struct { 
  Int_t    ipair;
  Float_t  spair;
  Float_t  slpair;
  Float_t  zintpa;
  Float_t  steppa;
  Int_t    icomp;
  Float_t  scomp;
  Float_t  slcomp;
  Float_t  zintco;
  Float_t  stepco;
  Int_t    iphot;
  Float_t  sphot;
  Float_t  slphot;
  Float_t  zintph;
  Float_t  stepph;
  Int_t    ipfis;
  Float_t  spfis;
  Float_t  slpfis;
  Float_t  zintpf;
  Float_t  steppf;
  Int_t    idray;
  Float_t  sdray;
  Float_t  sldray;
  Float_t  zintdr;
  Float_t  stepdr;
  Int_t    ianni;
  Float_t  sanni;
  Float_t  slanni;
  Float_t  zintan;
  Float_t  stepan;
  Int_t    ibrem;
  Float_t  sbrem;
  Float_t  slbrem;
  Float_t  zintbr;
  Float_t  stepbr;
  Int_t    ihadr;
  Float_t  shadr;
  Float_t  slhadr;
  Float_t  zintha;
  Float_t  stepha;
  Int_t    imunu;
  Float_t  smunu;
  Float_t  slmunu;
  Float_t  zintmu;
  Float_t  stepmu;
  Int_t    idcay;
  Float_t  sdcay;
  Float_t  slife;
  Float_t  sumlif;
  Float_t  dphys1;
  Int_t    iloss;
  Float_t  sloss;
  Float_t  soloss;
  Float_t  stloss;
  Float_t  dphys2;
  Int_t    imuls;
  Float_t  smuls;
  Float_t  somuls;
  Float_t  stmuls;
  Float_t  dphys3;
  Int_t    irayl;
  Float_t  srayl;
  Float_t  slrayl;
  Float_t  zintra;
  Float_t  stepra;
} Gcphys_t; 
 
//----------GCOPTI 
//      COMMON/GCOPTI/IOPTIM
typedef struct { 
  Int_t   ioptim;
} Gcopti_t; 
 
//----------GCTLIT 
//      COMMON/GCTLIT/THRIND,PMIN,DP,DNDL,JMIN,ITCKOV,IMCKOV,NPCKOV
typedef struct { 
  Float_t   thrind;
  Float_t   pmin;
  Float_t   dp;
  Float_t   dndl;
  Int_t     jmin;
  Int_t     itckov;
  Int_t     imckov;
  Int_t     npckov;
} Gctlit_t; 
 
//----------GCVDMA 
//      COMMON/GCVDMA/NVMANY,MANYLE(20),MANYNA(20,15),
//     +MANYNU(20,15),NFMANY,MYCOUN,IMYSE,RAYTRA,VECCOS(3)
typedef struct { 
  Int_t     vdma[624];
  Float_t   raytra;
  Float_t   veccos[3];
} Gcvdma_t; 
 
//----------GCTPOL 
#define MAXME1 30 
//      COMMON/GCTPOL/POLAR(3), NAMEC1(MAXME1) 
typedef struct { 
  Float_t polar[3]; 
  Int_t   namec1[MAXME1]; 
} Gctpol_t; 


class TGeant3 : public AliMC { 

private:
  Int_t fNextVol;    // Iterator for GeomIter

//--------------Declarations for ZEBRA--------------------- 
  Int_t *fZiq, *fZlq; 
  Float_t *fZq; 

  Quest_t  *fQuest; 
  Gcbank_t *fGcbank;
  Gclink_t *fGclink; 
  Gccuts_t *fGccuts; 
  Gcmate_t *fGcmate; 
  Gctpol_t *fGctpol; 
  Gcnum_t  *fGcnum; 
  Gcsets_t *fGcsets; 
  Gcopti_t *fGcopti; 
  Gctlit_t *fGctlit; 
  Gcvdma_t *fGcvdma; 
  Gcvolu_t *fGcvolu; 
  Gckine_t *fGckine; 
  Gcflag_t *fGcflag; 
  Gctmed_t *fGctmed; 
  Gcphys_t *fGcphys; 
  Gcking_t *fGcking; 
  Gckin2_t *fGckin2; 
  Gckin3_t *fGckin3; 
  Gctrak_t *fGctrak; 

  enum {kMaxParticles = 100};

  Int_t fNPDGCodes;

  Int_t fPDGCode[kMaxParticles];

public: 
  TGeant3(); 
  TGeant3(const char *title, Int_t nwgeant=0); 
  virtual ~TGeant3() {} 
  virtual void LoadAddress(); 
 
///////////////////////////////////////////////////////////////////////
//                                                                   //
//                                                                   //
//     Here are the service routines from the geometry               //
//     which could be implemented also in other geometries           //
//                                                                   //
//                                                                   //
///////////////////////////////////////////////////////////////////////

  void  GeomIter();
  Int_t CurrentMaterial(Float_t &a, Float_t &z, Float_t &dens, Float_t &radl, Float_t &absl) const;
  Int_t NextVolUp(Text_t *name, Int_t &copy);
  Int_t CurrentVolID(Int_t &copy) const;
  Int_t CurrentVolOffID(Int_t off, Int_t &copy) const;
  const char* CurrentVolName() const;
  const char *CurrentVolOffName(Int_t off) const;
  Int_t VolId(Text_t *name) const;
  Int_t IdFromPDG(Int_t pdg) const;
  Int_t PDGFromId(Int_t pdg) const;
  void  DefineParticles();
  const char* VolName(Int_t id) const;
  void  TrackPosition(TLorentzVector &xyz) const;
  void  TrackMomentum(TLorentzVector &xyz) const;  
  Int_t NofVolumes() const;
  Float_t TrackTime() const;  
  Float_t TrackCharge() const;
  Float_t TrackMass() const;
  Float_t TrackStep() const;
  Float_t TrackLength() const;
  Int_t   TrackPid() const;
  Bool_t IsTrackInside() const;
  Bool_t IsTrackEntering() const;
  Bool_t IsTrackExiting() const;
  Bool_t IsTrackOut() const;
  Bool_t IsTrackDisappeared() const;
  Bool_t IsTrackStop() const;
  Bool_t IsTrackAlive() const;
  Int_t   NSecondaries() const;
  Int_t   CurrentEvent() const;
  void    ProdProcess(char*) const;
  void    GetSecondary(Int_t, Int_t&, Float_t*, Float_t*);
  void   StopTrack();
  void   StopEvent();
  Float_t MaxStep() const;
  void   SetColors();
  void  SetMaxStep(Float_t maxstep);
  void  SetMaxNStep(Int_t maxnstp);
  Int_t GetMaxNStep() const;
  void GetParticle(const Int_t pdg, char *name, Float_t &mass) const;
  virtual Int_t GetMedium() const;
  virtual Float_t Edep() const;
  virtual Float_t Etot() const;
  virtual void    Rndm(Float_t* r, const Int_t n) const;
  virtual void    Material(Int_t&, const char*, Float_t, Float_t, Float_t, Float_t,
			    Float_t, Float_t* buf=0, Int_t nwbuf=0);
  virtual void    Mixture(Int_t&, const char*, Float_t*, Float_t*, Float_t, Int_t, Float_t*);
  virtual void    Medium(Int_t&, const char*, Int_t, Int_t, Int_t, Float_t, Float_t, 
		   Float_t, Float_t, Float_t, Float_t, Float_t* ubuf=0, Int_t nbuf=0);
  virtual void    Matrix(Int_t&, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);

/////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                         //
//                                                                                         //
//     Here are the interface functions with GEANT3.21                                     //
//                                                                                         //
//                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////

  // access functions to commons
 
  virtual Quest_t* Quest() const {return fQuest;}
  virtual Gcbank_t* Gcbank() const {return fGcbank;}
  virtual Gclink_t* Gclink() const {return fGclink;}
  virtual Gccuts_t* Gccuts() const {return fGccuts;}
  virtual Gcmate_t* Gcmate() const {return fGcmate;}
  virtual Gctpol_t* Gctpol() const {return fGctpol;}
  virtual Gcnum_t* Gcnum() const {return fGcnum;}
  virtual Gcsets_t* Gcsets() const {return fGcsets;}
  virtual Gcopti_t* Gcopti() const {return fGcopti;}
  virtual Gctlit_t* Gctlit() const {return fGctlit;}
  virtual Gcvdma_t* Gcvdma() const {return fGcvdma;}
  virtual Gcvolu_t* Gcvolu() const {return fGcvolu;}
  virtual Gckine_t* Gckine() const {return fGckine;}
  virtual Gcflag_t* Gcflag() const {return fGcflag;}
  virtual Gctmed_t* Gctmed() const {return fGctmed;}
  virtual Gcphys_t* Gcphys() const {return fGcphys;}
  virtual Gcking_t* Gcking() const {return fGcking;}
  virtual Gckin2_t* Gckin2() const {return fGckin2;}
  virtual Gckin3_t* Gckin3() const {return fGckin3;}
  virtual Gctrak_t* Gctrak() const {return fGctrak;}
  virtual Int_t* Iq() const {return fZiq;}
  virtual Int_t* Lq() const {return fZlq;}
  virtual Float_t* Q() const {return fZq;}


      // functions from GBASE 
   virtual  void  Gpcxyz(); 
   virtual  void  Ggclos(); 
   virtual  void  Gfile(const char *filename, const char *option="I"); 
   virtual  void  Glast(); 
   virtual  void  Gprint(const char *name); 
   virtual  void  Grun(); 
   virtual  void  Gtrig(); 
   virtual  void  Gtrigc(); 
   virtual  void  Gtrigi(); 
   virtual  void  Gwork(Int_t nwork); 
   virtual  void  Gzinit(); 
 
      // functions from GCONS 
   virtual  void  Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z, Float_t &dens, 
                         Float_t &radl, Float_t &absl, Float_t* ubuf, Int_t& nbuf); 
   virtual  void  Gfpart(Int_t ipart, char *name, Int_t &itrtyp,  
                         Float_t &amass, Float_t &charge, Float_t &tlife); 
   virtual  void  Gftmed(Int_t numed, char *name, Int_t &nmat, Int_t &isvol,  
                         Int_t &ifield, Float_t &fieldm, Float_t &tmaxfd, 
                         Float_t &stemax, Float_t &deemax, Float_t &epsil, 
                         Float_t &stmin, Float_t *buf=0, Int_t *nbuf=0); 
   virtual  void  Gmate(); 
   virtual  void  Gpart(); 
   virtual  void  Gsckov(Int_t itmed, Int_t npckov, Float_t *ppckov,
			 Float_t *absco, Float_t *effic, Float_t *rindex); 
   virtual  void  Gsdk(Int_t ipart, Float_t *bratio, Int_t *mode); 
   virtual  void  Gsmate(Int_t imat, const char *name, Float_t a, Float_t z,  
                         Float_t dens, Float_t radl, Float_t absl); 
   virtual  void  Gsmixt(Int_t imat, const char *name, Float_t *a, Float_t *z,  
                         Float_t dens, Int_t nlmat, Float_t *wmat); 
   virtual  void  Gspart(Int_t ipart, const char *name, Int_t itrtyp,  
                         Float_t amass, Float_t charge, Float_t tlife); 
   virtual  void  Gstmed(Int_t numed, const char *name, Int_t nmat, Int_t isvol,  
                         Int_t ifield, Float_t fieldm, Float_t tmaxfd, 
                         Float_t stemax, Float_t deemax, Float_t epsil, 
                         Float_t stmin); 
   virtual  void  Gstpar(Int_t itmed, const char *param, Float_t parval); 
 
      // functions from GKINE 
   virtual  void  Gfkine(Int_t itra, Float_t *vert, Float_t *pvert, 
                         Int_t &ipart, Int_t &nvert); 
   virtual  void  Gfvert(Int_t nvtx, Float_t *v, Int_t &ntbeam, Int_t &nttarg, Float_t &tofg); 
   virtual  Int_t Gskine(Float_t *plab, Int_t ipart, Int_t nv, Float_t *ubuf=0, Int_t nwbuf=0); 
   virtual  Int_t Gsvert(Float_t *v, Int_t ntbeam, Int_t nttarg, Float_t *ubuf=0, Int_t nwbuf=0); 
 
      // functions from GPHYS 
   virtual  void  Gphysi(); 
 
      // functions from GTRAK 
   virtual  void  Gdebug(); 
   virtual  void  Gekbin(); 
   virtual  void  Gfinds(); 
   virtual  void  Gsking(Int_t igk); 
   virtual  void  Gskpho(Int_t igk); 
   virtual  void  Gsstak(Int_t iflag); 
   virtual  void  Gsxyz(); 
   virtual  void  Gtrack(); 
   virtual  void  Gtreve(); 
   virtual  void  Gtreve_root(); 
   virtual  void  Grndm(Float_t *rvec, const Int_t len) const; 
   virtual  void  Grndmq(Int_t &is1, Int_t &is2, const Int_t iseq, const Text_t *chopt); 
 
      // functions from GGEOM 
   virtual  void  Gdxyz(Int_t ); 
   virtual  void  Gdcxyz(); 

      // functions from GGEOM 
   virtual  void  Gdtom(Float_t *xd, Float_t *xm, Int_t iflag); 
   virtual  void  Glmoth(const char* iudet, Int_t iunum, Int_t &nlev, 
                         Int_t *lvols, Int_t *lindx); 
   virtual  void  Gmedia(Float_t *x, Int_t &numed); 
   virtual  void  Gmtod(Float_t *xm, Float_t *xd, Int_t iflag); 
   virtual  void  Gsdvn(const char *name, const char *mother, Int_t ndiv, Int_t iaxis); 
   virtual  void  Gsdvn2(const char *name, const char *mother, Int_t ndiv, Int_t iaxis, Float_t c0i, Int_t numed); 
   virtual  void  Gsdvs(const char *name, const char *mother, Float_t step, Int_t iaxis, Int_t numed); 
   virtual  void  Gsdvs2(const char *name, const char *mother, Float_t step, Int_t iaxis, Float_t c0, Int_t numed); 
   virtual  void  Gsdvt(const char *name, const char *mother, Float_t step, Int_t iaxis, Int_t numed, Int_t ndvmx); 
   virtual  void  Gsdvt2(const char *name, const char *mother, Float_t step, Int_t iaxis,
			 Float_t c0, Int_t numed, Int_t ndvmx); 
   virtual  void  Gsord(const char *name, Int_t iax); 
   virtual  void  Gspos(const char *name, Int_t nr, const char *mother,  
                         Float_t x, Float_t y, Float_t z, Int_t irot, const char *konly="ONLY"); 
   virtual  void  Gsposp(const char *name, Int_t nr, const char *mother,  
                         Float_t x, Float_t y, Float_t z, Int_t irot, const char *konly, Float_t *upar, Int_t np); 
   virtual  void  Gsrotm(Int_t nmat, Float_t theta1, Float_t phi1, Float_t theta2, Float_t phi2, 
                         Float_t theta3, Float_t phi3); 
   virtual  void  Gprotm(Int_t nmat=0); 
   virtual  Int_t Gsvolu(const char *name, const char *shape, Int_t nmed,  
                         Float_t *upar, Int_t np); 
   virtual  void  Gsatt(const char *name, const char *att, Int_t val);
   virtual  void  Gfpara(const char *name, Int_t number, Int_t intext, Int_t& npar,
			 Int_t& natt, Float_t* par, Float_t* att);
   virtual  void  Gckpar(Int_t, Int_t, Float_t*);
   virtual  void  Gckmat(Int_t, char*);
    
      // functions from GDRAW 
   virtual  void  DefaultRange();
   virtual  void  InitHIGZ();
   virtual  void  Gdopen(Int_t view);
   virtual  void  Gdclose();
   virtual  void  Gdelete(Int_t view);
   virtual  void  Gdshow(Int_t view);
   virtual  void  Gdopt(const char *name,const char *value);
   virtual  void  Gdraw(const char *name,Float_t theta=30, Float_t phi=30, Float_t psi=0,Float_t u0=10,Float_t v0=10,Float_t ul=0.01,Float_t vl=0.01);
   virtual  void  Gdrawc(const char *name,Int_t axis=1, Float_t cut=0,Float_t u0=10,Float_t v0=10,Float_t ul=0.01,Float_t vl=0.01);
   virtual  void  Gdrawx(const char *name,Float_t cutthe, Float_t cutphi, Float_t cutval,
                         Float_t theta=30, Float_t phi=30,Float_t u0=10,Float_t v0=10,Float_t ul=0.01,Float_t vl=0.01);
   virtual  void  Gdhead(Int_t isel, const char *name, Float_t chrsiz=0.6);   
   virtual  void  Gdman(Float_t u0, Float_t v0, const char *type="MAN");
   virtual  void  Gdspec(const char *name);
   virtual  void  DrawOneSpec(const char *name);
   virtual  void  Gdtree(const char *name,Int_t levmax=15,Int_t ispec=0);
   virtual  void  GdtreeParent(const char *name,Int_t levmax=15,Int_t ispec=0);

   virtual  void  WriteEuclid(const char*, const char*, Int_t, Int_t);

   virtual  void  SetABAN(Int_t par=1);
   virtual  void  SetANNI(Int_t par=1);
   virtual  void  SetAUTO(Int_t par=1);
   virtual  void  SetBOMB(Float_t bomb=1);
   virtual  void  SetBREM(Int_t par=1);
   virtual  void  SetCKOV(Int_t par=1);
   virtual  void  SetClipBox(const char *name,Float_t xmin=-9999,Float_t xmax=0, Float_t ymin=-9999,Float_t ymax=0,Float_t zmin=-9999,Float_t zmax=0);
   virtual  void  SetCOMP(Int_t par=1);
   virtual  void  SetCUTS(Float_t cutgam,Float_t cutele,Float_t cutneu,Float_t cuthad,
                      Float_t cutmuo ,Float_t bcute ,Float_t bcutm ,Float_t dcute ,
                      Float_t dcutm ,Float_t ppcutm, Float_t tofmax);
   virtual  void  SetDCAY(Int_t par=1);
   virtual  void  SetDEBU(Int_t emin=1, Int_t emax=999, Int_t emod=1);
   virtual  void  SetDRAY(Int_t par=1);
   virtual  void  SetHADR(Int_t par=1);
   virtual  void  SetKINE(Int_t kine, Float_t xk1=0, Float_t xk2=0, Float_t xk3=0, Float_t xk4=0,
                         Float_t xk5=0, Float_t xk6=0, Float_t xk7=0, Float_t xk8=0, Float_t xk9=0,
                         Float_t xk10=0);
   virtual  void  SetLOSS(Int_t par=2);
   virtual  void  SetMULS(Int_t par=1);
   virtual  void  SetMUNU(Int_t par=1);
   virtual  void  SetOPTI(Int_t par=2);
   virtual  void  SetPAIR(Int_t par=1);
   virtual  void  SetPFIS(Int_t par=1);
   virtual  void  SetPHOT(Int_t par=1);
   virtual  void  SetRAYL(Int_t par=1);
   virtual  void  SetSWIT(Int_t sw, Int_t val=1);
   virtual  void  SetTRIG(Int_t nevents=1);
   virtual  void  SetUserDecay(Int_t ipart);

   virtual  void  Vname(const char *name, char *vname);

   virtual  void  InitLego();
        
   ClassDef(TGeant3,1)  //C++ interface to Geant basic routines 
}; 

#endif 
