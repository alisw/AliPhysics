#ifndef ALIANALYSISTASKEFFCONTPIDBF_cxx
#define ALIANALYSISTASKEFFCONTPIDBF_cxx

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency and contamination of the Balance 
// Function for single particles and pairs
// Modified By Noor Alam(VECC ,Kolkata)
//[ Special thanks to Michael Weber]
// ---------------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH2D;
class TH3D;
class TString;
class AliAODEvent;
class AliAODInputHandler;
class AliAODMCParticle;
class AliAODTrack;

#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEffContPIDBF : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEffContPIDBF() : AliAnalysisTaskSE(), 
    fAOD(0), 
    fPIDType(kNSigmaTPCTOF),
    fNSigmaPID(3.0),
    fArrayMC(0),  
    fQAList(0), 
    fOutputList(0),
    fHistEventStats(0),
    fHistCentrality(0), 
    fHistVz(0),
    fHistNSigmaTPCvsPtbeforePID(0),
    fHistNSigmaTPCvsPtafterPID(0),  
    fHistTruthPionPlus(0),
    fHistTruthKaonPlus(0),
    fHistTruthProtonPlus(0),
    fHistTruthPionMinus(0),
    fHistTruthKaonMinus(0),
    fHistTruthProtonMinus(0),
    fHistTruthPion(0),
    fHistTruthKaon(0),
    fHistTruthProton(0),
    HistMCTruthPtAll(0),
    HistMCTruthEtaAll(0),
    HistMCTruthPhiAll(0),
    HistMCTruthPtPion(0),
    HistMCTruthEtaPion(0),
    HistMCTruthPhiPion(0),
    HistMCTruthPtKaon(0),
    HistMCTruthEtaKaon(0),
    HistMCTruthPhiKaon(0),
    HistMCTruthPtProton(0),
    HistMCTruthEtaProton(0),
    HistMCTruthPhiProton(0),
    fHistSigmaTPCVsTOFPionForPionAfterCut(0),
    fHistSigmaTPCVsTOFProtonForPionAfterCut(0),
    fHistSigmaTPCVsTOFKaonForPionAfterCut(0),
    fHistSigmaTPCVsTOFPionForKaonAfterCut(0),
    fHistSigmaTPCVsTOFProtonForKaonAfterCut(0),
    fHistSigmaTPCVsTOFKaonForKaonAfterCut(0),
    fHistSigmaTPCVsTOFPionForProtonAfterCut(0),
    fHistSigmaTPCVsTOFProtonForProtonAfterCut(0),
    fHistSigmaTPCVsTOFKaonForProtonAfterCut(0),
    h1PionAfterCut(0),
    h1KaonAfterCut(0),
    h1ProtonAfterCut(0),
    h1PionAsNonPion(0),
    h1KaonAsNonKaon(0),
    h1ProtonAsNonProton(0),
    h1PionAsKaon(0),
    h1PionAsProton(0),
    h1KaonAsPion(0),
    h1KaonAsProton(0),
    h1ProtonAsPion(0),
    h1ProtonAsKaon(0),
    fHistMCRecoPionPlus(0),
    fHistMCRecoKaonPlus(0),
    fHistMCRecoProtonPlus(0),
    fHistMCRecoPionMinus(0),
    fHistMCRecoKaonMinus(0),
    fHistMCRecoProtonMinus(0),
    fHistMCRecoPion(0),
    fHistMCRecoKaon(0),
    fHistMCRecoProton(0),
    fHistMCRecoPionAsKaon(0),
    fHistMCRecoPionAsProton(0),
    fHistMCRecoProtonAsKaon(0),
    fHistMCRecoProtonAsPion(0),
    fHistMCRecoKaonAsPion(0),
    fHistMCRecoKaonAsProton(0),
    fUseCentrality(kFALSE), 
    fCentralityEstimator("V0M"), 
    fCentralityPercentileMin(0.0), 
    fCentralityPercentileMax(5.0), 
    fInjectedSignals(kFALSE),
    fPIDResponse(0),
    fElectronRejection(kFALSE),
    fElectronOnlyRejection(kFALSE),
    fElectronRejectionNSigma(-1.),
    fElectronRejectionMinPt(0.),
    fElectronRejectionMaxPt(1000.),
    fVxMax(3.0), 
    fVyMax(3.0), 
    fVzMax(10.), 
    fAODTrackCutBit(128),
    fMinNumberOfTPCClusters(80), 
    fMaxChi2PerTPCCluster(4.0),
    fMaxDCAxy(3.0), 
    fMaxDCAz(3.0),
    fMinPt(0.0), 
    fMaxPt(20.0), 
    fPtTPCMax(0.6),
    fMinEta(-0.8),
    fMaxEta(0.8), 
    fEtaRangeMin(0.0), 
    fEtaRangeMax(1.6), 
    fPtRangeMin(0.0), 
    fPtRangeMax(20.0), 
    fEtaBin(100),
    fdEtaBin(64),
    fPtBin(100),
    fHistdEdxTPC(0),
    fHistBetaTOF(0),
    fParticleType_(kPion){
    for(Int_t ipart=0;ipart<3;ipart++)
     for(Int_t ipid=0;ipid<3;ipid++)
      fnsigmas[ipart][ipid]=999.;
}

    AliAnalysisTaskEffContPIDBF(const char *name);
    virtual ~AliAnalysisTaskEffContPIDBF() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  Bool_t   IsLabelUsed(TArrayI array, Int_t label);

  void SetAODtrackCutBit(Int_t bit){
    fAODTrackCutBit = bit;
  }
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;
    fVyMax = vy;
    fVzMax = vz;
  }

  //Centrality
  void UseCentrality() { fUseCentrality = kTRUE;}
  void SetCentralityEstimator(const char* centralityEstimator) {
    fCentralityEstimator = centralityEstimator;}
  void SetCentralityPercentileRange(Float_t min, Float_t max) { 
    fCentralityPercentileMin=min;
    fCentralityPercentileMax=max;
  }

  //Injected signals
  void SetRejectInjectedSignals() {fInjectedSignals = kTRUE;}

  // electron rejection
  void SetElectronRejection(Double_t gMaxNSigma){
    fElectronRejection = kTRUE;
    fElectronRejectionNSigma = gMaxNSigma;
  }
  
  void SetElectronOnlyRejection(Double_t gMaxNSigma){
    fElectronRejection       = kTRUE;
    fElectronOnlyRejection   = kTRUE;
    fElectronRejectionNSigma = gMaxNSigma;
  }
  
  void SetElectronRejectionPt(Double_t minPt,Double_t maxPt){
    fElectronRejectionMinPt  = minPt;
    fElectronRejectionMaxPt  = maxPt;
  }

   enum ParticleSies_{kSpPion = 0,kSpKaon,kSpProton,kNSpecies,kSpUndefined=999};
   enum AliParticleSpecies_{kPion = 0,kKaon,kProton,kUndefined=999};
   enum PIDType_t{ kNSigmaTPC = 0, kNSigmaTPCTOF,kNSigmaTOF};


   void SetParticleType(AliParticleSpecies_ particletype_) {
    fParticleType_ = particletype_; }
    void SetPIDType(PIDType_t PIDType) { fPIDType = PIDType; }

  
  //Track cuts
  void SetMinNumberOfTPCClusters(Double_t min) {
    fMinNumberOfTPCClusters = min;}
  void SetMaxChi2PerTPCCluster(Double_t max) {
    fMaxChi2PerTPCCluster = max;}
  void SetMaxDCAxy(Double_t max) {
    fMaxDCAxy = max;}
  void SetMaxDCAz(Double_t max) {
    fMaxDCAz = max;}
  void SetMinPt(Double_t minPt) {
    fMinPt = minPt;}
  void SetMaxPt(Double_t maxPt) {
    fMaxPt = maxPt;}
 
  void SetEtaRange(Double_t minEta, Double_t maxEta, Int_t binEta, Double_t minRangeEta, Double_t maxRangeEta, Int_t bindEta){
    fMinEta = minEta;
    fMaxEta = maxEta;
    fEtaBin = binEta;
    fEtaRangeMax = maxRangeEta;
    fEtaRangeMin = minRangeEta;
    fdEtaBin = bindEta;
  }
  void SetPtRange(Double_t minRangePt, Double_t maxRangePt,Int_t binPt){
    fPtRangeMin = minRangePt;
    fPtRangeMax = maxRangePt;
    fPtBin = binPt;} 
  void SetTPCPtMax(Double_t ptmax){
     fPtTPCMax = ptmax;
      }
 void SetNSigmaCut(Double_t nsigma) { fNSigmaPID = nsigma; }

  Bool_t IsMCParticleCut(AliAODMCParticle* particle);
 
  Int_t GetParticleSpecies(AliAODTrack* trk );
  Double_t Beta(AliAODTrack *track);
  Bool_t IsTPCPID(AliAODTrack* track) const; //Check Track ,Is Track TPC?
  Bool_t IsTOFPID(AliAODTrack* track) const; // Check Track , Is Track TOF?
  void SigmaCalculate(AliAODTrack *trk );
  Int_t SigmaCutForParticleSpecies(AliAODTrack *trk );

 
 private:
  PIDType_t fPIDType; // PID type
  Double_t fNSigmaPID;
  AliAODEvent* fAOD; //! AOD object  
  TClonesArray *fArrayMC; //! array of MC particles  
  TList       *fQAList; //! QA list
  TList       *fOutputList; //! Output list
  
  // Truth MC Particle Pt , Eta abd Phi distribution Histogram 

   TH1F *HistMCTruthPtAll; // Pt All
   TH1F *HistMCTruthEtaAll; // Eta All
   TH1F *HistMCTruthPhiAll; // Phi All

   TH1F *HistMCTruthPtPion; //Pt Pion
   TH1F *HistMCTruthEtaPion; // Eta Pion
   TH1F *HistMCTruthPhiPion; // Phi Pion

   TH1F *HistMCTruthPtKaon; // Pt Kaon
   TH1F *HistMCTruthEtaKaon; // Eta Kaon
   TH1F *HistMCTruthPhiKaon; // Phi Kaon

   TH1F *HistMCTruthPtProton; // Pt Proton
   TH1F *HistMCTruthEtaProton; // Eta Proton
   TH1F *HistMCTruthPhiProton; // Phi Proton
   
   // For MC Reco 

   TH2F *fHistSigmaTPCVsTOFPionForPionAfterCut;
   TH2F *fHistSigmaTPCVsTOFProtonForPionAfterCut;
   TH2F *fHistSigmaTPCVsTOFKaonForPionAfterCut;

   TH2F *fHistSigmaTPCVsTOFPionForKaonAfterCut;
   TH2F *fHistSigmaTPCVsTOFProtonForKaonAfterCut;
   TH2F *fHistSigmaTPCVsTOFKaonForKaonAfterCut;

   TH2F *fHistSigmaTPCVsTOFPionForProtonAfterCut;
   TH2F *fHistSigmaTPCVsTOFProtonForProtonAfterCut;
   TH2F *fHistSigmaTPCVsTOFKaonForProtonAfterCut;

   TH1F *h1PionAfterCut; TH1F *h1KaonAfterCut; TH1F *h1ProtonAfterCut;
   TH1F *h1PionAsNonPion; TH1F *h1KaonAsNonKaon; TH1F *h1ProtonAsNonProton;
   TH1F *h1PionAsKaon;TH1F *h1PionAsProton;TH1F *h1KaonAsPion;TH1F *h1KaonAsProton; TH1F *h1ProtonAsPion;TH1F *h1ProtonAsKaon;

   






  // QA histograms
  TH1F        *fHistEventStats; //!event stats
  TH1F        *fHistCentrality; //!centrality
  TH1F        *fHistVz;//!
  TH2F        *fHistNSigmaTPCvsPtbeforePID;//TPC nsigma vs pT before PID cuts (QA histogram)
  TH2F        *fHistNSigmaTPCvsPtafterPID;//TPC nsigma vs pT after PID cuts (QA histogram)

  TH2F *fHistdEdxTPC;
  TH2F *fHistBetaTOF;

  // output histograms
  TH3D        *fHistTruthPionPlus;//!
  TH3D        *fHistTruthKaonPlus;//!
  TH3D        *fHistTruthProtonPlus;//!
  
  TH3D        *fHistTruthPionMinus;//!
  TH3D        *fHistTruthKaonMinus;//!
  TH3D        *fHistTruthProtonMinus;//!

  TH3D        *fHistTruthPion;
  TH3D        *fHistTruthKaon;
  TH3D        *fHistTruthProton;

  TH3D        *fHistMCRecoPionPlus;//!
  TH3D        *fHistMCRecoKaonPlus;//!
  TH3D        *fHistMCRecoProtonPlus;//!

  TH3D        *fHistMCRecoPionMinus;//!
  TH3D        *fHistMCRecoKaonMinus;//!
  TH3D        *fHistMCRecoProtonMinus;//!

  TH3D        *fHistMCRecoPion;
  TH3D        *fHistMCRecoKaon;
  TH3D        *fHistMCRecoProton;

  TH3D        *fHistMCRecoPionAsKaon;
  TH3D        *fHistMCRecoPionAsProton;
 
  TH3D        *fHistMCRecoProtonAsKaon;
  TH3D        *fHistMCRecoProtonAsPion;
  TH3D        *fHistMCRecoKaonAsPion;
  TH3D        *fHistMCRecoKaonAsProton;

  Bool_t  fUseCentrality;// Bool_t use centrality or not
  TString fCentralityEstimator;//"V0M","TRK","TKL","ZDC","FMD"
  Float_t fCentralityPercentileMin, fCentralityPercentileMax; //min-max centrality percentile

  Bool_t fInjectedSignals;//Flag for using the rejection of injected signals

  AliPIDResponse *fPIDResponse;     //! PID response object
  Bool_t   fElectronRejection;//flag to use electron rejection
  Bool_t   fElectronOnlyRejection;//flag to use electron rejection with exclusive electron PID (no other particle in nsigma range)
  Double_t fElectronRejectionNSigma;//nsigma cut for electron rejection
  Double_t fElectronRejectionMinPt;//minimum pt for electron rejection (default = 0.)
  Double_t fElectronRejectionMaxPt;//maximum pt for electron rejection (default = 1000.)

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax
  
  Int_t fAODTrackCutBit;//track cut bit from track selection (only used for AODs)

  Double_t fMinNumberOfTPCClusters;//!
  Double_t fMaxChi2PerTPCCluster;//!
  Double_t fMaxDCAxy, fMaxDCAz;//!
  Double_t fMinPt, fMaxPt,fPtTPCMax; //!
  Double_t fMinEta, fMaxEta;//!
  Double_t fEtaRangeMin;// acceptance cuts 
  Double_t fEtaRangeMax; // acceptance cuts
  Double_t fPtRangeMin;  //acceptance cuts
  Double_t fPtRangeMax;  //acceptance cuts
  
  Int_t fEtaBin;  //acceptance cuts
  Int_t fdEtaBin;  //acceptance cuts
  Int_t fPtBin; //acceptance cuts
  Double_t fnsigmas[3][3];
 
  AliParticleSpecies_ fParticleType_;

  AliAnalysisTaskEffContPIDBF(const AliAnalysisTaskEffContPIDBF&); // not implemented
  AliAnalysisTaskEffContPIDBF& operator=(const AliAnalysisTaskEffContPIDBF&); // not implemented
  
  ClassDef(AliAnalysisTaskEffContPIDBF, 1); // example of analysis
};

#endif
