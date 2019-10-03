#ifndef AliAnalysisTaskFemtoESE_cxx
#define AliAnalysisTaskFemtoESE_cxx


class TH1F;
class TH2F;
class TH3F;
class TH1D;
class TH2D;
class TH3D;

class TProfile;
class TProfile2D;

class AliESDEvent;
class AliAODEvent;
class AliESDtrackCuts;
class AliESDpid;
class AliHelperPID;
class AliEventPoolManager;
class AliSpectraAODEventCuts;
class AliSpectraAODTrackCuts;
class AliAODTrack;
class AliEventPool;

//class TStopwatch;

class AliFemtoESEBasicParticle;

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDpid.h"
#include "AliAODPid.h"

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliSpectraAODEventCuts.h"
#include "AliSpectraAODTrackCuts.h"

class AliAnalysisTaskFemtoESE : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskFemtoESE();
  AliAnalysisTaskFemtoESE(const char* name);
  virtual ~AliAnalysisTaskFemtoESE();
  AliAnalysisTaskFemtoESE(const AliAnalysisTaskFemtoESE &/*obj*/); 
  AliAnalysisTaskFemtoESE &operator=(const AliAnalysisTaskFemtoESE &/*obj*/);

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  void TrackLoop(TObjArray *tracks, AliEventPool *pool, Int_t z, Double_t psiEP, Float_t centralityPercentile);
  virtual void   Terminate(Option_t *);

  AliHelperPID* GetHelperPID() { return fHelperPID; }
  void SetHelperPID(AliHelperPID* pid){ fHelperPID = pid; }
  AliSpectraAODEventCuts* GetEventCuts() {return fEventCuts;}
  void SetEventCuts(AliSpectraAODEventCuts* cuts) {fEventCuts = cuts;}
  AliSpectraAODTrackCuts* GetTrackCuts() {return fTrackCuts;}
  void SetTrackCuts(AliSpectraAODTrackCuts* cuts) {fTrackCuts = cuts;}

  Int_t GetTrackFilterBit(){return fFilterBit;}
  void SetTrackFilterBit(Int_t bit){fFilterBit=bit;}
  void SetEventSelectionBit( UInt_t val ) {fSelectBit = val;}
  void SetIsLHC10h(Bool_t val) {bIsLHC10h = val;}
  void SetMinSepPair(Double_t eta,Double_t phi){fMinSepPairEta = eta; fMinSepPairPhi = phi;}
  void SetMinQinv(Double_t min){fQinvMin = min;}
  void SetMaxDCA(Double_t maxxy, Double_t maxz){fMaxDcaXY = maxxy; fMaxDcaZ = maxz;}
  void SetPtCuts(Double_t min, Double_t max){fPtMin = min; fPtMax = max;}
  void SetMaxEta(Double_t eta){fEtaMax = eta;}
  void SetShareFraction(Double_t val) {fShareFraction = val;}
  Double_t GetShareFraction() {return fShareFraction;}
  void SetShareQuality(Double_t val) {fShareQuality = val;}
  Double_t GetShareQuality() {return fShareQuality;}
  void SetQPercCuts(Double_t min, Double_t max) {fMinQPerc = min; fMaxQPerc = max;}
  Double_t GetMinQPerc(){return fMinQPerc;}
  Double_t GetMaxQPerc(){return fMaxQPerc;}
  void SetQPercDetector(Int_t det){fQPercDet = det;};
  void SetEPDetector(Int_t det){fEPDet = det;};
  void SetNMixingTracks(Int_t n){fMixingTracks = n;};
  void SetQBinning(Int_t n, Double_t q){qbins = n; qlimit = q;};

  void SetKtBins(Int_t n, Double_t* bins);
  void SetEPBins(Int_t n);
  void SetCentBins(Int_t n, Double_t* bins);
  void SetVzBins(Int_t n, Double_t* bins);

  //Double_t GetQPercLHC11h(Double_t qvec);

  Double_t GetCentralityWeight(Double_t cent);

 private:
  Double_t GetQinv(Double_t[], Double_t[]);
  void GetQosl(Double_t[], Double_t[], Double_t&, Double_t&, Double_t&);
  Bool_t TrackCut(AliAODTrack* ftrack);
  Bool_t EventCut(/*AliAODEvent* fevent*/);
  Bool_t PairCut(AliFemtoESEBasicParticle* ftrack1, AliFemtoESEBasicParticle* ftrack2, Bool_t mix);
  Double_t DeltaPhiStar(AliAODTrack* ftrack1, AliAODTrack* ftrack2, Double_t r);
  TObjArray* CloneAndReduceTrackList(TObjArray* tracks, Double_t psi);
  Double_t GetDeltaPhiEP(Double_t px1, Double_t py1, Double_t px2, Double_t py2, Double_t psi);
  Bool_t FindBin(Double_t kt, Double_t phi, Double_t cent, Int_t& a, Int_t& b, Int_t&c);

  AliAODEvent            *fAOD; //!    // AOD object
  TList                  *fOutputList; //! Compact Output list
  //AliPIDResponse         *fPIDResponse; //! PID response object; equivalent to AliAODpidUtil
  AliHelperPID*     fHelperPID;      // points to class for PID
  AliEventPoolManager**     fPoolMgr;         //![2] event pool manager
  AliSpectraAODEventCuts* fEventCuts;
  AliSpectraAODTrackCuts* fTrackCuts;

  Int_t          fFilterBit;         // track selection cuts
  UInt_t         fSelectBit;            // Select events according to AliAnalysisTaskJetServices bit maps 
  Bool_t         bIsLHC10h;
  Int_t fEventCounter;
  Int_t fMixingTracks;
  Double_t fBfield;
  Double_t fMinSepPairEta;
  Double_t fMinSepPairPhi;
  Double_t fQinvMin;
  Double_t fMaxDcaXY;
  Double_t fMaxDcaZ;
  Double_t fPtMin;
  Double_t fPtMax;
  Double_t fEtaMax;
  Double_t fShareQuality;
  Double_t fShareFraction;

  Int_t nCountSamePairs;
  Int_t nCountMixedPairs;
  Int_t nCountTracks;

  Double_t fMinQPerc;
  Double_t fMaxQPerc;

  Double_t qlimit;
  Int_t qbins;

  Int_t fQPercDet; // detector used for q-vector (0-V0A, 1-V0C)
  Int_t fEPDet; // detector used for event plane (0-V0A, 1-V0C)

  // binning for histograms
  Int_t nKtBins;
  Int_t nKtBins1;
  Double_t* ktBins; //[nKtBins1]
  Int_t nEPBins;
  Int_t nEPBins1;
  Double_t* epBins; //[nEPBins1]
  Int_t nEPBinsMix;
  Int_t nEPBinsMix1;
  Double_t* epBinsMix; //[nEPBinsMix1]
  Int_t nCentBins;
  Int_t nCentBins1;
  Double_t* centBins; //[nCentBins1]
  Int_t nVzBins;
  Int_t nVzBins1;
  Double_t* vzBins; //[nVzBins1]

  Double_t vertex[3];

  TH3F***** hq;
  TH3F***** hqmix;
  TH3F***** hqinv;

  Int_t nqPercBinsLHC11h;
  Double_t* qPercBinsLHC11h; //[nqPercBinsLHC11h]

  //TStopwatch *stopwatch;

  // histograms
  TH1D *hpx;
  TH1D *hpy;
  TH1D *hpz;
  TH1D *hpt;
  TH1D *hE;
  TH2D *hphieta;
  TH2D *hphieta_pid;
  TH1D *hpt_pid;
  TH2D *hvzcent;
  TH1D *hcent;
  TH1D* hcentUnweighted;
  TH2D *hcentn;
  TH3D *hphistaretapair10;
  TH3D *hphistaretapair16;
  TH3D *hphistaretapair10a;
  TH3D *hphistaretapair16a;
  TH3D *hphistaretapair10b;
  TH3D *hphistaretapair16b;
  TH3D *hphietapair;
  TH3D *hphietapair2;
  TH1D *hpidid;
  TH2D *hkt;
  TH1D *hktcheck;
  TH3D *hkt3;
  TH2D* hdcaxy;
  TH1D* hdcaz;
  TH1D* hsharequal;
  TH1D* hsharefrac;
  TH1D* hsharequalmix;
  TH1D* hsharefracmix;
  TH2D *hPsiTPC;
  TH2D *hPsiV0A;
  TH2D *hPsiV0C;
  TH2D* hShiftTPC;
  TH2D* hShiftV0A;
  TH2D* hShiftV0C;
  TH2D *hPsiMix;
  TH2D *hCheckEPA;
  TH2D *hCheckEPC;
  TH2D* hCheckEPmix;
  TH3D* hAvDphi;
  TH3D* hNpairs;
  TH1D* hPairDphi;
  TH1D* hPairDphiMix;
  TH2D *hcentq;
  TH2D* hMixedDistTracks;
  TH2D* hMixedDistEvents;
  TH2D *hQvecV0A;
  TH2D *hQvecV0C;
  TH1D *hresV0ATPC;
  TH1D *hresV0CTPC;
  TH1D *hresV0AV0C;
  TH3F* hqinvcheck;

  TH1F* hktbins;
  TH1F* hcentbins;
  TH1F* hepbins;

  ClassDef(AliAnalysisTaskFemtoESE, 1);
};


class AliFemtoESEBasicParticle : public AliAODTrack
{
 public:
 AliFemtoESEBasicParticle(Double_t En, Double_t px, Double_t py, Double_t pz, Short_t charge, Double_t phi, Double_t eta)
   : fE(En), fPx(px), fPy(py), fPz(pz), fCharge(charge), fPhi(phi), fEta(eta), fPsiEP(0.0)
  {
  }
  ~AliFemtoESEBasicParticle() {}

  // kinematics
  virtual Double_t Px() const { return fPx;}
  virtual Double_t Py() const { return fPy; }
  virtual Double_t Pz() const { return fPz; }
  virtual Double_t Pt() const { return sqrt(fPx*fPx+fPy*fPy); }
  virtual Double_t P()  const { return sqrt(fPx*fPx+fPy*fPy+fPz*fPz);; }
  virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
  
  virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
  virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
  
  virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
  virtual Double_t Phi()        const { return fPhi; }
  virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }
  
  virtual Double_t E()          const { return fE; }
  virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }
  
  virtual Double_t Eta()        const { return fEta; }
  virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }
  
  virtual Short_t Charge()      const { return fCharge; }
  virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
  // PID
  virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }      
  virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }

  void SetPsiEP(Double_t psi) {fPsiEP = psi;}
  Double_t GetPsiEP() {return fPsiEP;}

  /*void SetTPCClusterMap(TBits b){fClusterMap = TBits(b);};
  void SetTPCSharedMap(TBits b){fSharedMap = TBits(b);};
  TBits GetTPCClusterMap(){return fClusterMap;};
  TBits GetTPCSharedMap(){return fSharedMap;};*/
  
private:
  Double_t fE;
  Double_t fPx;
  Double_t fPy;
  Double_t fPz;
  Short_t fCharge;
  Double_t fPhi;
  Double_t fEta;
  Double_t fPsiEP;
  //TBits fClusterMap;
  //TBits fSharedMap;

  
  ClassDef( AliFemtoESEBasicParticle, 1);
};

#endif
