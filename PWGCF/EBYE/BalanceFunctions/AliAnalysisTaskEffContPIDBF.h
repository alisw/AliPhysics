#ifndef ALIANALYSISTASKEFFCONTPIDBF_cxx
#define ALIANALYSISTASKEFFCONTPIDBF_cxx

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency and contamination of the Balance 
// Function for single particles (for Identified particles) 
// Modified By Noor Alam (VECC ,Kolkata) sk.noor.alam@cern.ch
//[ Special thanks to Michael Weber ]
// ---------------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH2D;
class TH3F;
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
    HistPionContaminationInPt(0),
    HistPionPlusContaminationInPt(0),
    HistPionMinusContaminationInPt(0),
    Hist3dPionContamination(0),
    Hist3dPionPlusContamination(0),
    Hist3dPionMinusContamination(0),
    HistKaonContaminationInPt(0),
    HistKaonPlusContaminationInPt(0),
    HistKaonMinusContaminationInPt(0),
    Hist3dKaonContamination(0),
    Hist3dKaonPlusContamination(0),
    Hist3dKaonMinusContamination(0),
    HistProtonContaminationInPt(0),
    HistProtonPlusContaminationInPt(0),
    HistProtonMinusContaminationInPt(0),
    Hist3dProtonContamination(0),
    Hist3dProtonPlusContamination(0),
    Hist3dProtonMinusContamination(0), 
    HistPionPurityInPt(0),
    HistPionPlusPurityInPt(0),
    HistPionMinusPurityInPt(0),
    Hist3dPionPurity(0),
    Hist3dPionPlusPurity(0),
    Hist3dPionMinusPurity(0),
    HistKaonPurityInPt(0),
    HistKaonPlusPurityInPt(0),
    HistKaonMinusPurityInPt(0),
    Hist3dKaonPurity(0),
    Hist3dKaonPlusPurity(0),
    Hist3dKaonMinusPurity(0),
    HistProtonPurityInPt(0),
    HistProtonPlusPurityInPt(0),
    HistProtonMinusPurityInPt(0),
    Hist3dProtonPurity(0),
    Hist3dProtonPlusPurity(0),
    Hist3dProtonMinusPurity(0),
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
    fParticleType_(kPion),
    fHistNsigmaTPCPionBeforePIDCut(0),
 fHistNsigmaTPCKaonBeforePIDCut(0),
 fHistNsigmaTPCProtonBeforePIDCut(0),
 fHistNsigmaTOFPionBeforePIDCut(0),
 fHistNsigmaTOFKaonBeforePIDCut(0),
 fHistNsigmaTOFProtonBeforePIDCut(0),

 fHistNsigmaTPCPionAfterPIDCut(0),
 fHistNsigmaTPCKaonAfterPIDCut(0),
 fHistNsigmaTPCProtonAfterPIDCut(0),
 fHistNsigmaTOFPionAfterPIDCut(0),
 fHistNsigmaTOFKaonAfterPIDCut(0),
 fHistNsigmaTOFProtonAfterPIDCut(0),

 fHistNsigmaTPCTOFPionBeforePIDCut(0),
 fHistNsigmaTPCTOFKaonBeforePIDCut(0),
 fHistNsigmaTPCTOFProtonBeforePIDCut(0),
 fHistNsigmaTPCTOFPionAfterPIDCut(0),
 fHistNsigmaTPCTOFKaonAfterPIDCut(0),
 fHistNsigmaTPCTOFProtonAfterPIDCut(0),
 fHistdEdxTPCPionAfterPIDCut(0),
 fHistdEdxTPCKaonAfterPIDCut(0),
 fHistdEdxTPCProtonAfterPIDCut(0),
 fHistBetaTOFPionAfterPIDCut(0),
 fHistBetaTOFKaonAfterPIDCut(0),
 fHistBetaTOFProtonAfterPIDCut(0),
 fSigmaIndividually(kFALSE),
 fSigmaCutMethodOne(kFALSE)
{
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
 void SetSigmaIndividually(Bool_t sigmaIndividually){fSigmaIndividually=sigmaIndividually;}

 void SetSigmaCutMethodOne(Bool_t SigmaCutMethodOne){fSigmaCutMethodOne=SigmaCutMethodOne;}


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


// Contamination and Purity Histogram 

   TH1F *HistPionContaminationInPt; // Contamination Pion in Pt 
   TH1F *HistPionPlusContaminationInPt; // Contamination Pion + in Pt 
   TH1F *HistPionMinusContaminationInPt; // Contamination Pion - in Pt 
   TH3F *Hist3dPionContamination; // Contamination Pion in Eta , Pt , Phi
   TH3F *Hist3dPionPlusContamination; // Contamination Pion + in Eta , Pt , Phi
   TH3F *Hist3dPionMinusContamination; // Contamination Pion - in Eta , Pt , Phi

   TH1F *HistKaonContaminationInPt; // Contamination Kaon in Pt
   TH1F *HistKaonPlusContaminationInPt; // Contamination Kaon +  in Pt
   TH1F *HistKaonMinusContaminationInPt; // Contamination Kaon - in Pt
   TH3F *Hist3dKaonContamination; // Contamination Kaon in Eta , Pt , Phi
   TH3F *Hist3dKaonPlusContamination; // Contamination Kaon + in Eta , Pt , Phi
   TH3F *Hist3dKaonMinusContamination; // Contamination Kaon - in Eta , Pt , Phi

   TH1F *HistProtonContaminationInPt; // Contamination Proton in Pt
   TH1F *HistProtonPlusContaminationInPt; // Contamination Proton + in Pt 
   TH1F *HistProtonMinusContaminationInPt; // Contamination Proton - in Pt 
   TH3F *Hist3dProtonContamination; // Contamination Proton in Eta , Pt , Phi
   TH3F *Hist3dProtonPlusContamination; // Contamination Proton + in Eta , Pt , Phi
   TH3F *Hist3dProtonMinusContamination; // Contamination Proton - in Eta , Pt , Phi


   TH1F *HistPionPurityInPt; // Purity Pion in Pt 
   TH1F *HistPionPlusPurityInPt; // Purity Pion + in Pt 
   TH1F *HistPionMinusPurityInPt; // Purity Pion - in Pt 
   TH3F *Hist3dPionPurity; // Purity Pion in Eta , Pt , Phi
   TH3F *Hist3dPionPlusPurity; // Purity Pion + in Eta , Pt , Phi
   TH3F *Hist3dPionMinusPurity; // Purity Pion - in Eta , Pt , Phi

   TH1F *HistKaonPurityInPt; // Purity Kaon in Pt
   TH1F *HistKaonPlusPurityInPt; // Purity Kaon +  in Pt
   TH1F *HistKaonMinusPurityInPt; // Purity Kaon - in Pt
   TH3F *Hist3dKaonPurity; // Purity Kaon in Eta , Pt , Phi
   TH3F *Hist3dKaonPlusPurity; // Purity Kaon + in Eta , Pt , Phi
   TH3F *Hist3dKaonMinusPurity; // Purity Kaon - in Eta , Pt , Phi

   TH1F *HistProtonPurityInPt; // Purity Proton in Pt
   TH1F *HistProtonPlusPurityInPt; // Purity Proton + in Pt 
   TH1F *HistProtonMinusPurityInPt; // Purity Proton - in Pt 
   TH3F *Hist3dProtonPurity; // Purity Proton in Eta , Pt , Phi
   TH3F *Hist3dProtonPlusPurity; // Purity Proton + in Eta , Pt , Phi
   TH3F *Hist3dProtonMinusPurity; // Purity Proton - in Eta , Pt , Phi
// Contamination and Purity Histigram 


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

   
   // Nsigma Plot for Pion, Kaon and Proton

TH2F *fHistNsigmaTPCPionBeforePIDCut;
TH2F *fHistNsigmaTPCKaonBeforePIDCut;
TH2F *fHistNsigmaTPCProtonBeforePIDCut;
TH2F *fHistNsigmaTOFPionBeforePIDCut;
TH2F *fHistNsigmaTOFKaonBeforePIDCut;
TH2F *fHistNsigmaTOFProtonBeforePIDCut;

TH2F *fHistNsigmaTPCPionAfterPIDCut;
TH2F *fHistNsigmaTPCKaonAfterPIDCut;
TH2F *fHistNsigmaTPCProtonAfterPIDCut;
TH2F *fHistNsigmaTOFPionAfterPIDCut;
TH2F *fHistNsigmaTOFKaonAfterPIDCut;
TH2F *fHistNsigmaTOFProtonAfterPIDCut;

TH2F *fHistNsigmaTPCTOFPionBeforePIDCut;
TH2F *fHistNsigmaTPCTOFKaonBeforePIDCut;
TH2F *fHistNsigmaTPCTOFProtonBeforePIDCut;
TH2F *fHistNsigmaTPCTOFPionAfterPIDCut;
TH2F *fHistNsigmaTPCTOFKaonAfterPIDCut;
TH2F *fHistNsigmaTPCTOFProtonAfterPIDCut;

TH2F *fHistdEdxTPCPionAfterPIDCut;
TH2F *fHistdEdxTPCKaonAfterPIDCut;
TH2F *fHistdEdxTPCProtonAfterPIDCut;
TH2F *fHistBetaTOFPionAfterPIDCut;
TH2F *fHistBetaTOFKaonAfterPIDCut;
TH2F *fHistBetaTOFProtonAfterPIDCut;

  // QA histograms
  TH1F        *fHistEventStats; //!event stats
  TH1F        *fHistCentrality; //!centrality
  TH1F        *fHistVz;//!
  TH2F        *fHistNSigmaTPCvsPtbeforePID;//TPC nsigma vs pT before PID cuts (QA histogram)
  TH2F        *fHistNSigmaTPCvsPtafterPID;//TPC nsigma vs pT after PID cuts (QA histogram)

  TH2F *fHistdEdxTPC;
  TH2F *fHistBetaTOF;

  // output histograms
  TH3F        *fHistTruthPionPlus;//!
  TH3F        *fHistTruthKaonPlus;//!
  TH3F        *fHistTruthProtonPlus;//!
  
  TH3F        *fHistTruthPionMinus;//!
  TH3F        *fHistTruthKaonMinus;//!
  TH3F        *fHistTruthProtonMinus;//!

  TH3F        *fHistTruthPion;
  TH3F        *fHistTruthKaon;
  TH3F        *fHistTruthProton;

  TH3F        *fHistMCRecoPionPlus;//!
  TH3F        *fHistMCRecoKaonPlus;//!
  TH3F        *fHistMCRecoProtonPlus;//!

  TH3F        *fHistMCRecoPionMinus;//!
  TH3F        *fHistMCRecoKaonMinus;//!
  TH3F        *fHistMCRecoProtonMinus;//!

  TH3F        *fHistMCRecoPion;
  TH3F        *fHistMCRecoKaon;
  TH3F        *fHistMCRecoProton;

  TH3F        *fHistMCRecoPionAsKaon;
  TH3F        *fHistMCRecoPionAsProton;
 
  TH3F        *fHistMCRecoProtonAsKaon;
  TH3F        *fHistMCRecoProtonAsPion;
  TH3F        *fHistMCRecoKaonAsPion;
  TH3F        *fHistMCRecoKaonAsProton;

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

  Bool_t fSigmaIndividually;
  Bool_t fSigmaCutMethodOne;
 
  AliParticleSpecies_ fParticleType_;

  AliAnalysisTaskEffContPIDBF(const AliAnalysisTaskEffContPIDBF&); // not implemented
  AliAnalysisTaskEffContPIDBF& operator=(const AliAnalysisTaskEffContPIDBF&); // not implemented
  
  ClassDef(AliAnalysisTaskEffContPIDBF, 1); // example of analysis
};

#endif
