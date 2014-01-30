#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "AliLog.h"

#include "AliAnalysisTaskSE.h"

#define ID(x) x, #x
#define LAB(x) x + 1, #x

class TList;
class TClonesArray;
class AliVParticle;
class AliVTrack;
class AliPIDResponse;
class AliEventPoolManager;
class AliEventPool;
class AliESDtrackCuts;

class AliAnalysisTaskJetProtonCorr :
  public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskJetProtonCorr(const char *name = "jets_trg_trd");
  ~AliAnalysisTaskJetProtonCorr();

  // analysis operations
  virtual void   UserCreateOutputObjects();
  virtual Bool_t Notify();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t *option);

  // task configuration
  void SetJetBranchName(const char* const branchName) { strncpy(fJetBranchName, branchName, fgkStringLength-1); }
  const char* GetJetBranchName() const { return fJetBranchName; }

  void SetPtThrPart(Float_t minPt, Float_t maxPt) { fTrgPartPtMin = minPt; fTrgPartPtMax = maxPt; }
  Float_t GetPtMinPart() const { return fTrgPartPtMin; }
  Float_t GetPtMaxPart() const { return fTrgPartPtMax; }
  void SetPtThrJet(Float_t minPt, Float_t maxPt) { fTrgJetPtMin = minPt; fTrgJetPtMax = maxPt; }
  Float_t GetPtMinJet() const { return fTrgJetPtMin; }
  Float_t GetPtMaxJet() const { return fTrgJetPtMax; }
  void SetPtThrAss(Float_t minPt, Float_t maxPt) { fAssPartPtMin = minPt; fAssPartPtMax = maxPt; }
  Float_t GetPtMinAss() const { return fAssPartPtMin; }
  Float_t GetPtMaxAss() const { return fAssPartPtMax; }

  void PrintTask(Option_t *option, Int_t indent) const;

  static Double_t TOFsignal(Double_t *x, Double_t *par)
  {
    Double_t norm = par[0];
    Double_t mean = par[1];
    Double_t sigma = par[2];
    Double_t tail = par[3];

    if (x[0] <= (tail + mean))
      return norm * TMath::Gaus(x[0], mean, sigma);
    else
      return norm * TMath::Gaus(tail + mean, mean, sigma) * TMath::Exp(-tail * (x[0] - tail - mean) / (sigma * sigma));
  }

  // histograms
  enum Hist_t {
      kHistStat = 0,
      kHistCentrality,
      kHistCentralityUsed,
      kHistCentralityCheck,
      kHistCentralityCheckUsed,
      kHistSignalTPC,
      kHistSignalTOF,
      kHistBetaTOF,
      kHistDeltaTPC,
      kHistDeltaTPCSemi,
      kHistDeltaTOF,
      kHistDeltaTOFSemi,
      kHistExpSigmaTOFe,
      kHistExpSigmaTOFmu,
      kHistExpSigmaTOFpi,
      kHistExpSigmaTOFk,
      kHistExpSigmaTOFp,
      kHistExpSigmaTOFd,
      kHistExpSigmaTOFeSemi,
      kHistExpSigmaTOFmuSemi,
      kHistExpSigmaTOFpiSemi,
      kHistExpSigmaTOFkSemi,
      kHistExpSigmaTOFpSemi,
      kHistExpSigmaTOFdSemi,
      kHistCmpSigmaTOFe,
      kHistCmpSigmaTOFmu,
      kHistCmpSigmaTOFpi,
      kHistCmpSigmaTOFk,
      kHistCmpSigmaTOFp,
      kHistCmpSigmaTOFd,
      kHistCmpSigmaTOFeSemi,
      kHistCmpSigmaTOFmuSemi,
      kHistCmpSigmaTOFpiSemi,
      kHistCmpSigmaTOFkSemi,
      kHistCmpSigmaTOFpSemi,
      kHistCmpSigmaTOFdSemi,
      kHistNsigmaTPCe,
      kHistNsigmaTPCmu,
      kHistNsigmaTPCpi,
      kHistNsigmaTPCk,
      kHistNsigmaTPCp,
      kHistNsigmaTPCd,
      kHistNsigmaTPCe_e,
      kHistNsigmaTOFe,
      kHistNsigmaTOFmu,
      kHistNsigmaTOFpi,
      kHistNsigmaTOFk,
      kHistNsigmaTOFp,
      kHistNsigmaTOFd,
      kHistNsigmaTOFmismatch,
      kHistNsigmaTOFmismatch2,
      kHistDeltaTOFe,
      kHistDeltaTOFmu,
      kHistDeltaTOFpi,
      kHistDeltaTOFk,
      kHistDeltaTOFp,
      kHistDeltaTOFd,
      kHistNsigmaTPCeSemi,
      kHistNsigmaTPCmuSemi,
      kHistNsigmaTPCpiSemi,
      kHistNsigmaTPCkSemi,
      kHistNsigmaTPCpSemi,
      kHistNsigmaTPCdSemi,
      kHistNsigmaTPCe_eSemi,
      kHistNsigmaTOFeSemi,
      kHistNsigmaTOFmuSemi,
      kHistNsigmaTOFpiSemi,
      kHistNsigmaTOFkSemi,
      kHistNsigmaTOFpSemi,
      kHistNsigmaTOFdSemi,
      kHistNsigmaTOFmismatchSemi,
      kHistNsigmaTOFmismatch2Semi,
      kHistDeltaTOFeSemi,
      kHistDeltaTOFmuSemi,
      kHistDeltaTOFpiSemi,
      kHistDeltaTOFkSemi,
      kHistDeltaTOFpSemi,
      kHistDeltaTOFdSemi,
      kHistNsigmaTPCTOF,
      kHistNsigmaTPCTOFPt,
      kHistNsigmaTPCTOFUsed,
      kHistNsigmaTPCTOFUsedCentral,
      kHistNsigmaTPCTOFUsedSemiCentral,
      kHistNsigmaTPCTOFUsedPt,
      kHistNsigmaTPCTOFUsedPtCentral,
      kHistNsigmaTPCTOFUsedPtSemiCentral,
      kHistEvPlane,
      kHistEvPlaneUsed,
      kHistEvPlaneCheck,
      kHistEvPlaneCheckUsed,
      kHistEvPlaneCorr,
      kHistJetPtCentral,
      kHistJetPtSemi,
      kHistEtaPhiTrgHad,
      kHistEtaPhiTrgJet,
      kHistEtaPhiAssHad,
      kHistEtaPhiAssProt,
      kHistLast
  };

  // statistics
  enum Stat_t {
      kStatSeen = 1,
      kStatTrg,
      kStatCent,
      kStatEvPlane,
      kStatPID,
      kStatUsed,
      kStatEvCuts,
      kStatCentral,
      kStatSemiCentral,
      kStatLast
  };

  // trigger conditions
  enum Trigger_t {
      kTriggerMB = 0,
      kTriggerInt,
      kTriggerLast
  };

  // classification
  enum CorrType_t {
    kCorrHadHad = 0,
    kCorrHadProt,
    kCorrJetHad,
    kCorrJetProt,
    kCorrLast
  };

  enum Class_t {
    kClCentral = 0,
    kClSemiCentral,
    kClDijet,
    kClLast
  };

  enum Trg_t {
    kTrgHad = 0,
    kTrgJet,
    kTrgLast
  };

  enum Ass_t {
    kAssHad,
    kAssProt,
    kAssLast
  };

  enum Ev_t {
    kEvSame = 0,
    kEvMix,
    kEvLast
  };

  class AliHistCorr : public TNamed {
  public:
    AliHistCorr(TString name, TList *outputList = 0x0);
    ~AliHistCorr();

    void Trigger(Float_t weight = 1.) { fHistStat->Fill(1., weight); }
    void Fill(AliVParticle *trgPart, AliVParticle *assPart, Float_t weight = 1.);

  protected:
    TList *fOutputList;

    TH1F *fHistStat;

    TH1F *fHistCorrPhi;
    TH2F *fHistCorrPhi2;
    TH2F *fHistCorrEtaPhi;

    AliHistCorr(const AliHistCorr &rhs);
    AliHistCorr& operator=(const AliHistCorr &rhs);

    ClassDef(AliHistCorr, 1);
  };

  AliHistCorr*& GetHistCorr(CorrType_t corr, Class_t cl, Ev_t ev) { return fHistCorr[kEvLast*(kClLast*corr + cl) + ev]; }
  AliEventPoolManager*& GetPoolMgr(Trg_t trg, Ass_t ass) { return fPoolMgr[kTrgLast * ass + trg]; }
  AliEventPool*& GetPool(Class_t cls, Trg_t trg, Ass_t ass) { return fPool[kClLast * (kTrgLast * ass + trg) + cls]; }

protected:
  AliMCEvent  *fMCEvent; //!
  AliESDEvent *fESDEvent; //!
  AliAODEvent *fAODEvent; //!

  UInt_t fTriggerMask;		//! internal representation of trigger conditions
  UInt_t fClassMask;		//! internal representation of event classes
  Float_t fCentrality; //!
  Float_t fCentralityCheck; //!
  Float_t fZvtx; //!
  AliPIDResponse *fPIDResponse; //!
  Float_t fEventPlane; //!
  Float_t fEventPlaneCheck; //!
  TObjArray *fPrimTrackArray; //!
  TClonesArray *fJetArray; //!

  AliEventPoolManager *fPoolMgr[kTrgLast * kAssLast]; //!
  AliEventPool *fPool[kClLast * kTrgLast * kAssLast]; //!

  AliHistCorr **fHistCorr; //! [kCorrLast*kEvLast*kClLast]; //!

  Bool_t DetectTriggers();
  void   MarkTrigger(Trigger_t trg) { fTriggerMask |= (1 << trg); }
  Bool_t IsTrigger(Trigger_t trg) const { return (fTriggerMask & (1 << trg)); }

  Bool_t DetectClasses();
  void   MarkClass(Class_t cl) { fClassMask |= (1 << cl); }
  Bool_t IsClass(Class_t cl) const { return (fClassMask & (1 << cl)); }

  Bool_t PrepareEvent();
  Bool_t CleanUpEvent();

  Float_t GetCentrality() const { return fCentrality; }
  Float_t GetEventPlane() const { return fEventPlane; }
  AliPIDResponse* GetPID() const { return fPIDResponse; }
  Bool_t IsCentral() { return ((fCentrality >= 0.) && (fCentrality <= 10.)); }
  Bool_t IsSemiCentral() { return ((fCentrality >= 30.) && (fCentrality <= 50.)); }

  AliVTrack* GetLeadingTrack(AliAODJet *jet) const;

  Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1,
		      Float_t phi2, Float_t pt2, Float_t charge2,
		      Float_t radius, Float_t bSign);

  Bool_t AcceptTrigger(AliVTrack *trg);
  Bool_t AcceptTrigger(AliAODJet *trg);
  Bool_t AcceptAssoc(AliVTrack *trk);
  Bool_t IsProton(AliVTrack *trk);
  Bool_t AcceptAngleToEvPlane(Float_t phi, Float_t psi);
  Bool_t AcceptTwoTracks(AliVParticle *trgPart, AliVParticle *assPart);

  TObjArray* CloneTracks(TObjArray *tracks) const;

  Bool_t Correlate(CorrType_t corr, Class_t cl, Ev_t ev,
		   TCollection *trgArray, TCollection *assArray, Float_t weight = 1.);

  Bool_t Correlate(Trg_t trg, Ass_t ass, Class_t cl, Ev_t ev,
		   TCollection *trgArray, TCollection *assArray, Float_t weight = 1.);

  // output objects
  TList *fOutputList;		//! list of output objects

  // histogram management
  TH1  *fHist[kHistLast];	//! pointers to histogram
  const char *fShortTaskId;	//! short identifier for the task

  TH1*&  GetHistogram(Hist_t hist, const Int_t idx = 0) { return fHist[hist + idx]; }

  TH1*   AddHistogram(Hist_t hist, const char *hid, TString title,
                      Int_t xbins, Float_t xmin, Float_t xmax, Int_t binType = 1);
  TH2*   AddHistogram(Hist_t hist, const char *hid, TString title,
                      Int_t xbins, Float_t xmin, Float_t xmax,
                      Int_t ybins, Float_t ymin, Float_t ymax, Int_t binType = 1);
  TH3*   AddHistogram(Hist_t hist, const char *hid, TString title,
                      Int_t xbins, Float_t xmin, Float_t xmax,
                      Int_t ybins, Float_t ymin, Float_t ymax,
                      Int_t zbins, Float_t zmin, Float_t zmax, Int_t binType = 1);

  void    FillH1(Hist_t hist, Float_t x, Float_t weight = 1., Int_t idx = 0)
  { GetHistogram(hist, idx)->Fill(x, weight); }
  void    FillH2(Hist_t hist, Float_t x, Float_t y, Float_t weight = 1., Int_t idx = 0)
  { ((TH2*) GetHistogram(hist, idx))->Fill(x, y, weight); }
  void    FillH3(Hist_t hist, Float_t x, Float_t y, Float_t z, Float_t weight = 1., Int_t idx = 0)
  { ((TH3*) GetHistogram(hist, idx))->Fill(x, y, z, weight); }

  const char* fkCorrTypeName[kCorrLast]; //!
  const char* fkClassName[kClLast]; //!
  const char* fkEvName[kEvLast]; //!

  // task configuration
  static const Int_t fgkStringLength = 100; // max length for the jet branch name
  char fJetBranchName[fgkStringLength];     // jet branch name

  Bool_t fUseStandardCuts;

  AliESDtrackCuts *fCutsPrim;	// track cuts for primary particles
  Float_t fCutsTwoTrackEff;

  Float_t fTrgPartPtMin;
  Float_t fTrgPartPtMax;
  Float_t fTrgJetPtMin;
  Float_t fTrgJetPtMax;
  Float_t fTrgJetLeadTrkPtMin;
  Float_t fAssPartPtMin;
  Float_t fAssPartPtMax;
  Float_t fTrgAngleToEvPlane;

  // not implemented
  AliAnalysisTaskJetProtonCorr(const AliAnalysisTaskJetProtonCorr &rhs);
  AliAnalysisTaskJetProtonCorr& operator=(const AliAnalysisTaskJetProtonCorr &rhs);

  ClassDef(AliAnalysisTaskJetProtonCorr, 1);
};
