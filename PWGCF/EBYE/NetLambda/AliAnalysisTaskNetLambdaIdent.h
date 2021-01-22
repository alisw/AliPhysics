#ifndef AliAnalysisTaskNetLambdaIdent_h
#define AliAnalysisTaskNetLambdaIdent_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task for net-lambda fluctuations analysis
// Author: Alice Ohlson (alice.ohlson@cern.ch)

#include "AliAnalysisTaskSE.h"
//#include "AliNetLambdaHelper.h"
class TList;
class AliESDEvent;
class AliESDtrack;
class AliAODTrack;
//class AliExternalTrackParam;
class AliAnalysisUtils;
class AliPIDResponse;
class TTree;
class TH1;
class TH2;
class TH3;
class TH3F;
class TOBjArray;
class AliVVertex;
class AliAODv0;
class AliESDv0;
class AliEventPoolManager;
class AliLightV0;
//class AliAnalysisTaskWeakDecayVertexer;
#include "AliEventCuts.h"
#include "AliExternalTrackParam.h"
#include "TArrayF.h"

class AliAnalysisTaskNetLambdaIdent : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNetLambdaIdent(const char* name="AliAnalysisTaskNetLambdaIdent");
  virtual ~AliAnalysisTaskNetLambdaIdent(){};
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void Tracks2V0vertices(TClonesArray* fv0s, TObjArray* ev1, TObjArray* ev2, Bool_t mixing, AliVVertex *vtxT3D, Double_t b);
  Bool_t TrackCutsForTreeAOD(AliAODTrack* trk);
  Bool_t TrackCutsForTreeESD(AliESDtrack* trk);
  Bool_t V0CutsForTreeAOD(AliAODv0* v0, Double_t* vt);
  Bool_t V0CutsForTreeESD(AliESDv0* v0, Double_t* vt, AliExternalTrackParam* ptrk, AliExternalTrackParam* ntrk, Double_t b);
  Int_t IsGenLambda(Int_t poslabel, Int_t neglabel, Float_t pt, Float_t &mpt, Float_t &meta, Float_t &cascpt, Float_t &casceta);
  Double_t GetDCAV0Dau( AliExternalTrackParam *pt, AliExternalTrackParam *nt, Double_t &xp, Double_t &xn, Double_t b, Double_t lNegMassForTracking, Double_t lPosMassForTracking);
  void GetHelixCenter(const AliExternalTrackParam *track,Double_t center[2], Double_t b);
  void Evaluate(const Double_t *h, Double_t t, Double_t r[3], Double_t g[3], Double_t gg[3]);
  // event cuts
  void SetCentCut(Float_t cut){centcut = cut;};
  // daughter track cuts
  void SetCrossedRowsCut(Float_t rows, Float_t ratio){ncrossedrowscut = rows; crossedrowsclustercut = ratio;};
  void SetImpParamCut(Double_t cut){fDmin = cut;};
  // v0 cuts
  void SetPtMinLambda(Float_t pt){ptminlambda = pt;};
  void SetPtMaxLambda(Float_t pt){ptmaxlambda = pt;};
  void SetEtaCutLambda(Float_t eta){etacutlambda = eta;};
  void SetCosPACut(Double_t cut){fCPAmin = cut;};
  void SetRadiusCuts(Double_t mincut, Double_t maxcut){fRmin = mincut; fRmax = maxcut;};
  void SetDCADaughterCut(Double_t cut){fDCAmax = cut;};
  void SetChi2Cut(Double_t cut){fChi2max = cut;};
  void SetMinPionDCA(Float_t cut){minpiondca = cut;};
  void SetGammaCuts(Float_t cutmass, Float_t cutpt){gammamasscut = cutmass; gammaptcut = cutpt;};
  
  void SetIsMC(Bool_t val){fIsMC = val;};
  Bool_t GetIsMC(){return fIsMC;};
  void SetIsAOD(Bool_t val){fIsAOD = val;};
  Bool_t GetIsAOD(){return fIsAOD;};
  void SetEventSelection(UInt_t val) {fEvSel = val;}
  void MakeLambdaTree(Bool_t val){fLambdaTree = val;};
  void MakeEventMixingTree(Bool_t val){fEventMixingTree = val;};
  void DoRevertexing(Bool_t val){fRevertex = val;};
  void SetMaxMixTracks(Int_t val){nmaxmixtracks = val;};
  void SetMaxMixEvents(Int_t val){nmaxmixevents = val;};

 protected:
  AliAnalysisTaskNetLambdaIdent(const  AliAnalysisTaskNetLambdaIdent &task);
  AliAnalysisTaskNetLambdaIdent& operator=(const  AliAnalysisTaskNetLambdaIdent &task);
  
  AliESDEvent* fESD;
  AliAODEvent* fAOD;
  //AliAnalysisUtils* fUtils;   //! analysis utils to detect pileup
  AliPIDResponse* fPIDResponse; // points to class for PID
  AliEventCuts fEventCuts;      /// Event cuts
  TList* fListOfHistos;         //! list of output histograms
  TTree* fTree;                 //! output tree

  TH1I*  hEventStatistics;      //! cut-by-cut counter of events

  AliEventPoolManager*     fPoolMgr;         // event pool manager
  
  TH1F* hGenPt;
  TH1F* hGenPhi;
  TH1F* hGenEta;
  TH1F* hTrackPt;
  TH1F* hTrackPhi;
  TH1F* hTrackEta;
  TH1F* hNTracksEsd;
  TH2F* hNSigmaProton;
  TH2F* hArmPod;

  TH2F* hVxVy;
      
  TH1F*  hLambdaPtGen;
  TH1F*  hAntiLambdaPtGen;
  TH1F*  hLambdaPtReco;
  TH1F*  hAntiLambdaPtReco;
  
  TH2F*  hInvMassLambda;
  TH2F*  hInvMassAntiLambda;
  TH2F*  hInvMassLambdaOnTheFly;
  TH2F*  hInvMassAntiLambdaOnTheFly;
  TH2F*  hInvMassLambdaReco;
  TH2F*  hInvMassAntiLambdaReco;
  TH2F*  hInvMassLambdaSecFromMaterial;
  TH2F*  hInvMassAntiLambdaSecFromMaterial;
  TH2F*  hInvMassLambdaSecFromWeakDecay;
  TH2F*  hInvMassAntiLambdaSecFromWeakDecay;
  TH2F*  hInvMassLike;
  TH2F*  hInvMassLambdaK0S;
  TH2F*  hInvMassAntiLambdaK0S;
  TH2F*  hInvMassLambdaConversions;
  TH2F*  hInvMassAntiLambdaConversions;
  //TH3F*  hInvMassPtPidLambda;
  //TH3F*  hInvMassPtPidAntiLambda;

  TH3F* hXiPlus;
  TH3F* hXiMinus;
  TH3F* hXiZero;
  TH3F* hXiZeroAnti;

  TH2F* hPtResLambda;
  TH2F* hPtResAntiLambda;
  TH2F* hPtResLambdaPrim;
  TH2F* hPtResAntiLambdaPrim;
  
  // kinematic cuts
  Float_t centcut;
  Float_t ptminlambda;
  Float_t ptmaxlambda;
  Float_t etacutlambda;
  Float_t ncrossedrowscut;
  Float_t crossedrowsclustercut;

  Float_t fCentV0M;
  Float_t fCentCL1;
  Float_t fMultV0M;
  Float_t fNtracksTPCout;
  Double_t fVtxZ;
  Int_t fRunNumber;
  TClonesArray *fAcceptV0; //!
  //TClonesArray *fAcceptV0test; //!
  TClonesArray *fGenLambda; //!
  TClonesArray *fGenCascade; //!
  TClonesArray *fMixV0; //!

  Bool_t fIsMC;
  Bool_t fIsAOD;
  UInt_t fEvSel;
  Bool_t fLambdaTree;
  Bool_t fEventMixingTree;
  Bool_t fRevertex;
  Int_t nmaxmixtracks;
  Int_t nmaxmixevents;

  Double_t fChi2max;      //max chi2
  Double_t fDmin;//0.05;  //min imp parameter for the 1st daughter
  Double_t fDCAmax;//1.5; //max DCA between the daughter tracks
  Double_t fCPAmin;//0.9; //min cosine of V0's pointing angle
  Double_t fRmin;         //min radius of the fiducial volume
  Double_t fRmax;//200.;  //max radius of the fiducial volume
  Float_t minpiondca;
  Float_t gammamasscut;
  Float_t gammaptcut;
  
  const Float_t massPi = 0.139570;
  const Float_t massP = 0.938272;
  const Float_t massLambda = 1.115683;

  //AliMCEvent*              fMcEvent;    //! MC event
  //AliInputEventHandler*    fMcHandler;  //! MCEventHandler 
 
  ClassDef(AliAnalysisTaskNetLambdaIdent,12);
};

//_____________________________________________________________________________
class AliLightV0 : public TObject, public TArrayF
{
 public:
 AliLightV0() : TObject(), TArrayF(), pt(-999), eta(-999), invmass(-999), invmassK0S(-999), invmassGamma(-999),
    cospt(-999), decayr(-999), proplife(-999), dcav0pv(-999), dcadaughters(-999), mcstatus(0),
    ppt(-999), peta(-999), pnsigmapr(-999), pdca(-999), npt(-999), neta(-999), nnsigmapr(-999), ndca(-999) {};//genpt(-999), geneta(-999), cascpt(-999), casceta(-999) {};
 AliLightV0(Float_t ptin, Float_t etain, UInt_t ndim = 0) : TObject(), TArrayF(ndim), pt(ptin), eta(etain), invmass(-999), invmassK0S(-999), invmassGamma(-999),
    cospt(-999), decayr(-999), proplife(-999), dcav0pv(-999), dcadaughters(-999), mcstatus(0), 
    ppt(-999), peta(-999), pnsigmapr(-999), pdca(-999), npt(-999), neta(-999), nnsigmapr(-999), ndca(-999) {};//genpt(-999), geneta(-999), cascpt(-999), casceta(-999) {};
 AliLightV0(Float_t ptin, Float_t etain, Float_t invmassin, Float_t invmassK0Sin, Float_t invmassGammain,  
	       Float_t cosptin, Float_t decayrin, Float_t proplifein, UInt_t ndim = 0) : TObject(), TArrayF(ndim), pt(ptin), eta(etain), invmass(invmassin), invmassK0S(invmassK0Sin), invmassGamma(invmassGammain),
    cospt(cosptin), decayr(decayrin), proplife(proplifein), dcav0pv(-999), dcadaughters(-999), mcstatus(0),
    ppt(-999), peta(-999), pnsigmapr(-999), pdca(-999), npt(-999), neta(-999), nnsigmapr(-999), ndca(-999) {};//genpt(-999), geneta(-999), cascpt(-999), casceta(-999) {};
  virtual ~AliLightV0(){};
  void SetPt(Float_t val){pt = val;};
  void SetEta(Float_t val){eta = val;};
  void SetInvMass(Float_t val){invmass = val;};
  void SetInvMassK0S(Float_t val){invmassK0S = val;};
  void SetInvMassGamma(Float_t val){invmassGamma = val;};
  void SetCosPointingAngle(Float_t val){cospt = val;};
  void SetDecayR(Float_t val){decayr = val;};
  void SetProperLifetime(Float_t val){proplife = val;};
  void SetDCAV0(Float_t val){dcav0pv = val;};
  void SetDCADaughters(Float_t val){dcadaughters = val;};
  void SetMcStatus(Int_t val){mcstatus = val;};
  void SetPosDaughter(Float_t ptin, Float_t etain, Float_t nsigma, Float_t dca){ppt = ptin; peta = etain; pnsigmapr = nsigma; pdca = dca;};
  void SetNegDaughter(Float_t ptin, Float_t etain, Float_t nsigma, Float_t dca){npt = ptin; neta = etain; nnsigmapr = nsigma; ndca = dca;};
  //void SetGenPtEta(Float_t ptin, Float_t etain){genpt = ptin; geneta = etain;};
  //void SetCascadePtEta(Float_t ptin, Float_t etain){cascpt = ptin; casceta = etain;};
  
  Float_t GetPt(){return pt;};
  Float_t GetEta(){return eta;};
  Float_t GetInvMass(){return invmass;};
  Float_t GetInvMassK0S(){return invmassK0S;};
  Float_t GetInvMassGamma(){return invmassGamma;};
  Float_t GetCosPointingAngle(){return cospt;};
  Float_t GetDecayR(){return decayr;};
  Float_t GetProperLifetime(){return proplife;};
  Float_t GetDCAV0(){return dcav0pv;};
  Float_t GetDCADaughters(){return dcadaughters;};
  Int_t   GetMcStatus(){return mcstatus;};
  void    GetPosDaughter(Float_t& ptout, Float_t& etaout, Float_t& nsigma, Float_t& dca){ptout = ppt; etaout = peta; nsigma = pnsigmapr; dca = pdca;};
  void    GetNegDaughter(Float_t& ptout, Float_t& etaout, Float_t& nsigma, Float_t& dca){ptout = npt; etaout = neta; nsigma = nnsigmapr; dca = ndca;};
  //Float_t GetGenPt(){return genpt;};
  //Float_t GetGenEta(){return geneta;};
  //Float_t GetCascadePt(){return cascpt;};
  //Float_t GetCascadeEta(){return casceta;};

 private:
  Float_t   pt;
  Float_t   eta;
  Float_t   invmass;
  Float_t   invmassK0S;
  Float_t   invmassGamma;
  Float_t   cospt;
  Float_t   decayr;
  Float_t   proplife;
  Float_t   dcav0pv;
  Float_t   dcadaughters;
  Int_t    mcstatus;
  Float_t   ppt; // positive daughter properties
  Float_t   peta;
  Float_t   pnsigmapr;
  Float_t   pdca;
  Float_t   npt; // negative daughter properties
  Float_t   neta;
  Float_t   nnsigmapr;
  Float_t   ndca;
  //Float_t   genpt;
  //Float_t   geneta;
  //Float_t   cascpt;
  //Float_t   casceta;
  
  ClassDef(AliLightV0, 9);
};

//_____________________________________________________________________________
class AliLightGenV0 : public TObject
{
 public:
 AliLightGenV0() : pt(-999), eta(-999), id(-999), ppt(-999), peta(-999), npt(-999), neta(-999) {};
 AliLightGenV0(Float_t ptin, Float_t etain, Int_t idin) : pt(ptin), eta(etain), id(idin), ppt(-999), peta(-999), npt(-999), neta(-999) {};
  virtual ~AliLightGenV0(){};
  void SetPt(Float_t val){pt = val;};
  void SetEta(Float_t val){eta = val;};
  void SetId(Int_t val){id = val;}
  void SetPosDaughter(Float_t ptin, Float_t etain){ppt = ptin; peta = etain;};
  void SetNegDaughter(Float_t ptin, Float_t etain){npt = ptin; neta = etain;};
  
  Float_t GetPt(){return pt;};
  Float_t GetEta(){return eta;};
  Int_t   GetId(){return id;};
  void    GetPosDaughter(Float_t& ptout, Float_t& etaout){ptout = ppt; etaout = peta;};
  void    GetNegDaughter(Float_t& ptout, Float_t& etaout){ptout = npt; etaout = neta;};

 private:
  Float_t   pt;
  Float_t   eta;
  Int_t     id;
  Float_t   ppt; // positive daughter properties
  Float_t   peta;
  Float_t   npt; // negative daughter properties
  Float_t   neta;
  
  ClassDef(AliLightGenV0, 3);
};

//_____________________________________________________________________________
class AliLightV0track : public TObject
{
 public:
 AliLightV0track() : extparam(), prpid(-999), vtxx(-999), vtxy(-999), vtxz(-999), massForTracking(-999), mclabel(-999) {};
 AliLightV0track(AliExternalTrackParam& inparam, Float_t inpid, Double_t* vertin, Double_t inm) : extparam(inparam), prpid(inpid),
    vtxx(vertin[0]), vtxy(vertin[1]), vtxz(vertin[2]), massForTracking(inm), mclabel(-999){};
  virtual ~AliLightV0track(){};

  void SetMcLabel(Int_t labelin){mclabel = labelin;};

  //void Set(AliVTrack* trk, Float_t val){extparam = AliExternalParam(trk); prpid = val;};
  AliExternalTrackParam *GetExtParam(){return &extparam;};
  Float_t GetProtonPID(){return prpid;};
  void GetPrimaryVertex(Double_t &xin, Double_t &yin, Double_t &zin){xin = vtxx; yin = vtxy; zin = vtxz;};
  Double_t GetMassForTracking(){return massForTracking;};
  Int_t GetMcLabel(){return mclabel;};
  
 private:
  AliExternalTrackParam extparam;
  Float_t prpid;
  Double_t vtxx;
  Double_t vtxy;
  Double_t vtxz;
  Double_t massForTracking;
  Int_t mclabel;
  
  ClassDef(AliLightV0track, 4);
};
#endif
