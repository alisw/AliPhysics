#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

#include "AliLog.h"

#include "AliVParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"

#define ID(x) x, #x
#define LAB(x) x + 1, #x

class TList;
class TClonesArray;
class AliOADBContainer;
class AliTOFPIDParams;
class AliVTrack;
class AliPIDResponse;
class AliEventPoolManager;
class AliEventPool;
class AliEventplane;
class TSpline;

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

  void SetParamsTOF();

  // task configuration
  void SetJetBranchName(const char* branchName) { strncpy(fJetBranchName, branchName, fgkStringLength-1); }
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

  void SetTwoTrackCut(Float_t cut) { fCutsTwoTrackEff = cut; }
  Float_t GetTwoTrackCut() const { return fCutsTwoTrackEff; }

  void SetTrackCutsAss(const AliESDtrackCuts &cuts) { *fCutsPrimAss = cuts; }
  void SetTrackCutsTrg(const AliESDtrackCuts &cuts) { *fCutsPrimTrg = cuts; }
  void SetTrackCutsTrgConstrain(const AliESDtrackCuts &cuts) { *fCutsPrimTrgConstrain = cuts; }

  void SetTrgJetEtaMax(Float_t etamax) { fTrgJetEtaMax = etamax; }
  Float_t GetTrgJetEtaMax() const { return fTrgJetEtaMax; }
  void SetHadEtaMax(Float_t etamax) { fHadEtaMax = etamax; }
  Float_t GetHadEtaMax() const { return fHadEtaMax; }

  void SetUseEvplaneV0(Bool_t usev0 = kTRUE) { fUseEvplaneV0 = usev0; }
  Bool_t GetUseEvplaneV0() const { return fUseEvplaneV0; }

  void SetJetV2(Float_t v2Cent, Float_t v2Semi) {
    fTrgJetPhiModCent->SetParameter(0, v2Cent);
    fTrgJetPhiModSemi->SetParameter(0, v2Semi);
  }
  void SetHadV2(Float_t v2Cent, Float_t v2Semi) {
    fTrgHadPhiModCent->SetParameter(0, v2Cent);
    fTrgHadPhiModSemi->SetParameter(0, v2Semi);
  }

  void SetLeadingTrackPtMin(Float_t pt) { fTrgJetLeadTrkPtMin = pt; }
  void SetLeadingTrackPtMax(Float_t pt) { fTrgJetLeadTrkPtMax = pt; }

  void SetJetAreaMin(Float_t area) { fTrgJetAreaMin = area; }
  Float_t GetJetAreaMin() const { return fTrgJetAreaMin; }

  void SetRequirePID(Bool_t req = kTRUE) { fRequirePID = req; }
  Bool_t GetRequirePID() const { return fRequirePID; }

  void SetFilterMask(Int_t mask) { fAssFilterMask = mask; }
  Int_t GetFilterMask() const { return fAssFilterMask; }

  void SetErrorCount(Int_t cnt) { fErrorMsg = cnt; }
  Int_t GetErrorCount() const { return fErrorMsg; }

  void SetTrgAngleToEvPlane(Float_t angle) { fTrgAngleToEvPlane = angle; }
  Float_t GetTrgAngleToEvPlane() const { return fTrgAngleToEvPlane; }

  void SetToyMeanNoPart(Float_t mean) { fToyMeanNoPart = mean; }
  Float_t GetToyMeanNoPart() const { return fToyMeanNoPart; }
  void SetToyRadius(Float_t radius) { fToyRadius = radius; }
  Float_t GetToyRadius() const { return fToyRadius; }
  void SetToySmearPhi(Float_t sigma) { fToySmearPhi = sigma; }
  Float_t GetToySmearPhi() const { return fToySmearPhi; }

  void SetEventPlaneResSpline(TSpline *spline) { fSplineEventPlaneRes = spline; }
  const TSpline *GetEventPlaneResSpline() const { return fSplineEventPlaneRes; }

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
      kHistVertexNctb,
      kHistVertexZ,
      kHistCentrality,
      kHistCentralityUsed,
      kHistCentralityCheck,
      kHistCentralityCheckUsed,
      kHistCentralityVsMult,
      kHistSignalTPC,
      kHistSignalTOF,
      kHistBetaTOF,
      kHistDeltaTPC,
      kHistDeltaTPCSemi,
      kHistDeltaTOF,
      kHistDeltaTOFSemi,
      // kHistExpSigmaTOFe,
      // kHistExpSigmaTOFmu,
      // kHistExpSigmaTOFpi,
      // kHistExpSigmaTOFk,
      // kHistExpSigmaTOFp,
      // kHistExpSigmaTOFd,
      // kHistExpSigmaTOFeSemi,
      // kHistExpSigmaTOFmuSemi,
      // kHistExpSigmaTOFpiSemi,
      // kHistExpSigmaTOFkSemi,
      // kHistExpSigmaTOFpSemi,
      // kHistExpSigmaTOFdSemi,
      // kHistCmpSigmaTOFe,
      // kHistCmpSigmaTOFmu,
      // kHistCmpSigmaTOFpi,
      // kHistCmpSigmaTOFk,
      // kHistCmpSigmaTOFp,
      // kHistCmpSigmaTOFd,
      // kHistCmpSigmaTOFeSemi,
      // kHistCmpSigmaTOFmuSemi,
      // kHistCmpSigmaTOFpiSemi,
      // kHistCmpSigmaTOFkSemi,
      // kHistCmpSigmaTOFpSemi,
      // kHistCmpSigmaTOFdSemi,
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

      kHistNsigmaTPCTOFUsedCentralMCe,
      kHistNsigmaTPCTOFUsedCentralMCmu,
      kHistNsigmaTPCTOFUsedCentralMCpi,
      kHistNsigmaTPCTOFUsedCentralMCk,
      kHistNsigmaTPCTOFUsedCentralMCp,
      kHistNsigmaTPCTOFUsedCentralMCd,

      kHistNsigmaTPCTOFUsedSemiCentralMCe,
      kHistNsigmaTPCTOFUsedSemiCentralMCmu,
      kHistNsigmaTPCTOFUsedSemiCentralMCpi,
      kHistNsigmaTPCTOFUsedSemiCentralMCk,
      kHistNsigmaTPCTOFUsedSemiCentralMCp,
      kHistNsigmaTPCTOFUsedSemiCentralMCd,

      kHistNevMix,

      kHistEvPlane,
      kHistEvPlaneRes,
      kHistEvPlaneUsed,
      kHistEvPlaneCheck,
      kHistEvPlaneCheckUsed,
      kHistEvPlane3,
      kHistEvPlaneCorr,
      kHistEvPlaneCross,
      kHistEvPlaneCorrNoTrgJets,
      kHistEvPlaneCorrNoTrgJetsTrgd,
      kHistJetPtCentral,
      kHistJetPtSemi,
      kHistEtaPhiTrgHad,
      kHistEtaPhiTrgJet,
      kHistEtaPhiAssHad,
      kHistEtaPhiAssProt,
      kHistPhiTrgJetEvPlane,
      kHistPhiTrgHadEvPlane,
      kHistPhiRndTrgJetEvPlane,
      kHistPhiRndTrgHadEvPlane,
      kHistPhiAssHadEvPlane,
      kHistPhiAssHadVsEvPlane,
      kHistPhiAssProtEvPlane,
      kHistPhiTrgJetEvPlane3,
      kHistPhiTrgHadEvPlane3,
      kHistPhiAssHadEvPlane3,
      kHistPhiAssProtEvPlane3,
      kHistLast
  };

  // statistics
  enum Stat_t {
      kStatSeen = 1,
      kStatTrg,
      kStatVtx,
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

    kCorrRndHadHad,
    kCorrRndHadProt,
    kCorrRndJetHad,
    kCorrRndJetProt,

    kCorrRndHadExcHad,
    kCorrRndHadExcProt,
    kCorrRndJetExcHad,
    kCorrRndJetExcProt,

    kCorrLast
  };

  enum Class_t {
    kClCentral = 0,
    kClSemiCentral,
    // kClDijet,
    kClLast
  };

  enum Trg_t {
    kTrgHad = 0,
    kTrgJet,
    kTrgHadRnd,
    kTrgJetRnd,
    kTrgLast
  };

  enum Ass_t {
    kAssHad = 0,
    kAssProt,
    kAssHadJetExc,
    kAssProtJetExc,
    kAssHadHadExc,
    kAssProtHadExc,
    kAssLast
  };

  enum Ev_t {
    kEvSame = 0,
    kEvMix,
    kEvLast
  };

  class AliPartCorr : public AliVParticle {
  public:
    AliPartCorr(Float_t eta = 0., Float_t phi = 0., Float_t pt = 0., Float_t charge = 0) :
      fPt(pt), fEta(eta), fPhi(phi), fCharge(charge) {}
    AliPartCorr(const AliVParticle &rhs) :
      fPt(rhs.Pt()), fEta(rhs.Eta()), fPhi(rhs.Phi()), fCharge(rhs.Charge()) {}
    virtual ~AliPartCorr() {}
    
    // kinematics
    virtual Double_t Px() const { AliFatal("not implemented"); return 0; }
    virtual Double_t Py() const { AliFatal("not implemented"); return 0; }
    virtual Double_t Pz() const { AliFatal("not implemented"); return 0; }
    virtual Double_t Pt() const { return fPt; }
    virtual Double_t P()  const { AliFatal("not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("not implemented"); return 0; }

    virtual Double_t Xv() const { AliFatal("not implemented"); return 0; }
    virtual Double_t Yv() const { AliFatal("not implemented"); return 0; }
    virtual Double_t Zv() const { AliFatal("not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("not implemented"); return 0; }

    virtual Double_t OneOverPt()  const { AliFatal("not implemented"); return 0; }
    virtual Double_t Phi()        const { return fPhi; }
    virtual Double_t Theta()      const { AliFatal("not implemented"); return 0; }

    virtual Double_t E()          const { AliFatal("not implemented"); return 0; }
    virtual Double_t M()          const { AliFatal("not implemented"); return 0; }

    virtual Double_t Eta()        const { return fEta; }
    virtual Double_t Y()          const { AliFatal("not implemented"); return 0; }

    virtual Short_t Charge()      const { return fCharge; }
    virtual Int_t   GetLabel()    const { AliFatal("not implemented"); return 0; }

    virtual Int_t   PdgCode()     const { AliFatal("not implemented"); return 0; }
    virtual const Double_t *PID() const { AliFatal("not implemented"); return 0; }

  protected:
    Float_t fPt;
    Float_t fEta;
    Float_t fPhi;
    Short_t fCharge;

    ClassDef(AliPartCorr, 1);
  };

  class AliHistCorr : public TNamed {
  public:
    AliHistCorr(TString name, TList *outputList = 0x0);
    ~AliHistCorr();

    void Trigger(Float_t phi, Float_t eta, Float_t qpt, Float_t weight) {
      fHistStat->Fill(1., weight);
      if (fHistCorrTrgEtaPhi)
	fHistCorrTrgEtaPhi->Fill(phi, eta, weight);
      if (fHistCorrTrgEtaPhiQpt)
	fHistCorrTrgEtaPhiQpt->Fill(phi, eta, qpt, weight);
    }
    void Ass(Float_t phi, Float_t eta, Float_t qpt, Float_t weight) {
      if (fHistCorrAssEtaPhi)
	fHistCorrAssEtaPhi->Fill(phi, eta, weight);
      if (fHistCorrAssEtaPhiQpt)
	fHistCorrAssEtaPhiQpt->Fill(phi, eta, qpt, weight);
    }
    void Fill(AliVParticle *trgPart, AliVParticle *assPart, Float_t weight = 1.);
    void Fill(TLorentzVector *trgPart, AliVParticle *assPart, Float_t weight = 1.);
    void Fill(TLorentzVector *trgPart, TLorentzVector *assPart, Float_t weight = 1.);

  protected:
    TList *fOutputList;

    TH1F *fHistStat;

    TH1F *fHistCorrPhi;
    TH2F *fHistCorrPhi2;
    TH2F *fHistCorrEtaPhi;
    TH2F *fHistCorrAvgEtaPhi;
    TH2F *fHistCorrTrgEtaPhi;
    TH2F *fHistCorrAssEtaPhi;
    TH3F *fHistCorrTrgEtaPhiQpt;
    TH3F *fHistCorrAssEtaPhiQpt;

    const Float_t fHistDphiLo;
    const Int_t   fHistDphiNbins;
    const Int_t   fHistDetaNbins;

    AliHistCorr(const AliHistCorr &rhs);
    AliHistCorr& operator=(const AliHistCorr &rhs);

    ClassDef(AliHistCorr, 1);
  };

  AliHistCorr*& GetHistCorr(CorrType_t corr, Class_t cl, Ev_t ev) { return fHistCorr[kEvLast*(kClLast*corr + cl) + ev]; }
  AliEventPoolManager*& GetPoolMgr(Ass_t ass) { return fPoolMgr[ass]; }
  AliEventPool*& GetPool(Ass_t ass) { return fPool[ass]; }

protected:
  AliMCEvent  *fMCEvent; //!
  AliESDEvent *fESDEvent; //!
  AliAODEvent *fAODEvent; //!

  Int_t fRunNumber; //! current run number
  AliOADBContainer *fOADBContainerTOF; //! container for OADB entry with TOF parameters
  AliTOFPIDParams *fParamsTOF; //! TOF parametrization

  AliEventplane *fEventplane; //! pointer to the event plane

  UInt_t fTriggerMask;		//! internal representation of trigger conditions
  UInt_t fClassMask;		//! internal representation of event classes
  Float_t fCentrality; //!
  Float_t fCentralityCheck; //!
  Float_t fZvtx; //!
  AliPIDResponse *fPIDResponse; //!
  Float_t fEventPlaneAngle; //!
  Float_t fEventPlaneRes; //!
  Float_t fEventPlaneAngleCheck; //!
  Float_t fEventPlaneAngle3; //!

  TObjArray *fPrimTrackArrayAss; //!
  TObjArray *fPrimTrackArrayTrg; //!
  TClonesArray *fPrimConstrainedTrackArray; //!
  TClonesArray *fJetArray; //!

  AliEventPoolManager *fPoolMgr[kAssProt + 1]; //!
  AliEventPool *fPool[kAssProt + 1]; //!

  AliHistCorr **fHistCorr; //! [kCorrLast*kEvLast*kClLast]; //!

  Int_t fErrorMsg; //! remaining number of error messages to be printed

  Bool_t DetectTriggers();
  void   MarkTrigger(Trigger_t trg) { fTriggerMask |= (1 << trg); }
  Bool_t IsTrigger(Trigger_t trg) const { return (fTriggerMask & (1 << trg)); }

  Bool_t DetectClasses();
  void   MarkClass(Class_t cl) { fClassMask |= (1 << cl); }
  Bool_t IsClass(Class_t cl) const { return (fClassMask & (1 << cl)); }

  Bool_t PrepareEvent();
  Bool_t CleanUpEvent();

  Float_t GetCentrality() const { return fCentrality; }
  Float_t GetEventPlaneAngle() const { return fEventPlaneAngle; }
  AliPIDResponse* GetPID() const { return fPIDResponse; }
  Bool_t IsCentral() const { return ((fCentrality >= 0.) && (fCentrality <= 10.)); }
  Bool_t IsSemiCentral() const { return ((fCentrality >= 30.) && (fCentrality <= 50.)); }

  AliVTrack* GetLeadingTrack(const AliAODJet *jet) const;

  Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1,
		      Float_t phi2, Float_t pt2, Float_t charge2,
		      Float_t radius, Float_t bSign) const;

  Bool_t AcceptTrigger(AliVTrack *trg);
  Bool_t AcceptTrigger(AliAODJet *trg);
  Bool_t AcceptAssoc(const AliVTrack *trk) const;
  Bool_t IsProton(const AliVTrack *trk) const;
  Bool_t AcceptAngleToEvPlane(Float_t phi, Float_t psi) const;
  Bool_t AcceptTwoTracks(const AliVParticle *trgPart, const AliVParticle *assPart) const;

  Float_t GetPhiRel2(AliVParticle *part) const;

  TObjArray* CloneTracks(TObjArray *tracks) const;

  Bool_t GenerateRandom(TCollection *trgJetArray, TCollection *trgHadArray,
			TCollection *assHadJetArray, TCollection *assProtJetArray,
			TCollection *assHadHadArray, TCollection *assProtHadArray,
			Float_t pFraction = .5);

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

  const Bool_t fUseStandardCuts;
  Bool_t fUseEvplaneV0;

  AliESDtrackCuts *fCutsPrimTrg;	// track cuts for primary particles (trigger)
  AliESDtrackCuts *fCutsPrimTrgConstrain;	// track cuts for primary particles (trigger)
  AliESDtrackCuts *fCutsPrimAss;	// track cuts for primary particles (associate)
  Float_t fCutsTwoTrackEff;

  UInt_t  fAssFilterMask;
  Bool_t  fRequirePID;
  Float_t fTrgJetEtaMax;
  Float_t fHadEtaMax;
  Float_t fTrgPartPtMin;
  Float_t fTrgPartPtMax;
  Float_t fTrgJetPtMin;
  Float_t fTrgJetPtMax;
  Float_t fTrgJetLeadTrkPtMin;
  Float_t fTrgJetLeadTrkPtMax;
  Float_t fTrgJetAreaMin;
  Float_t fAssPartPtMin;
  Float_t fAssPartPtMax;
  Float_t fTrgAngleToEvPlane;

  Float_t fToyMeanNoPart;
  Float_t fToyRadius;
  Float_t fToySmearPhi;

  TF1 *fTrgJetPhiModCent;
  TF1 *fTrgJetPhiModSemi;
  TF1 *fTrgHadPhiModCent;
  TF1 *fTrgHadPhiModSemi;

  Float_t fTrgJetV2Cent;
  Float_t fTrgJetV2Semi;
  Float_t fTrgHadV2Cent;
  Float_t fTrgHadV2Semi;

  TSpline *fSplineEventPlaneRes;

  // not implemented
  AliAnalysisTaskJetProtonCorr(const AliAnalysisTaskJetProtonCorr &rhs);
  AliAnalysisTaskJetProtonCorr& operator=(const AliAnalysisTaskJetProtonCorr &rhs);

  ClassDef(AliAnalysisTaskJetProtonCorr, 1);
};
