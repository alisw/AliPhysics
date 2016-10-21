#ifndef ALIANALYSISTASKEFFCONTPIDBF_cxx
#define ALIANALYSISTASKEFFCONTPIDBF_cxx

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency and contamination of the Balance 
// Function for single particles (for Identified particles) 
// By Noor Alam (VECC ,Kolkata) sk.noor.alam@cern.ch and Subhasis Chattopadhyay(sub.chattopadhyay@gmail.com)
//[ Special thanks to Michael Weber  ]
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
   AliAnalysisTaskEffContPIDBF(): AliAnalysisTaskSE(),
    fAOD(0),
    fNSigmaPID(3.0),
    fArrayMC(0), 
    fQAList(0), 
    fOutputList(0), 
    fHistEventStats(0), 
    fHistCentrality(0),
    fHistVz(0), 
    fHistNSigmaTPCvsPtbeforeElectronPID(0),
    fHistNSigmaTPCvsPtafterElectronPID(0),  
    fHistTruthPionPlus(0),
    fHistTruthKaonPlus(0),
    fHistTruthProtonPlus(0),
    fHistTruthPionMinus(0),
    fHistTruthKaonMinus(0),
    fHistTruthProtonMinus(0),
    fHistTruthPion(0),
    fHistTruthKaon(0),
    fHistTruthProton(0),
    Hist3dPionContamination(0),
    Hist3dPionPlusContamination(0),
    Hist3dPionMinusContamination(0),
    Hist3dKaonContamination(0),
    Hist3dKaonPlusContamination(0),
    Hist3dKaonMinusContamination(0),
    Hist3dProtonContamination(0),
    Hist3dProtonPlusContamination(0),
    Hist3dProtonMinusContamination(0), 
    fHistMCRecoPionPlus(0),
    fHistMCRecoKaonPlus(0),
    fHistMCRecoProtonPlus(0),
    fHistMCRecoPionMinus(0),
    fHistMCRecoKaonMinus(0),
    fHistMCRecoProtonMinus(0),
    fHistMCRecoPion(0),
    fHistMCRecoKaon(0),
    fHistMCRecoProton(0),
    fHistNsigmaTPCPionBeforePIDCut(0),
    fHistNsigmaTPCKaonBeforePIDCut(0),
    fHistNsigmaTPCProtonBeforePIDCut(0),
    fHistNsigmaTOFPionBeforePIDCut(0),
    fHistNsigmaTOFKaonBeforePIDCut(0),
    fHistNsigmaTOFProtonBeforePIDCut(0),
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
    fEtaBin(100), //=100 (BF) 16
    fdEtaBin(64), //=64 (BF)  16
    fPtBin(100), //=100 (BF)  36
    fHistdEdxTPC(0),
    fHistBetaTOF(0),
    fParticleType_(kPion_),
    fDetectorPID_(kTPCTOFpid_),
    fRapidityInsteadOfEta(kFALSE),
    fMistMatchTOFProb(0.0),
    fTOFMisMatch(kFALSE),
    fZvertexTruthPion(0),
    fZvertexTruthPionPlus(0),
    fZvertexTruthPionMinus(0),
    fZvertexContaminationPion(0),
    fZvertexContaminationPionPlus(0),
    fZvertexContaminationPionMinus(0),
    fZvertexRecoPion(0),
    fZvertexRecoPionPlus(0),
    fZvertexRecoPionMinus(0),
    fZvertexTruthKaon(0),
    fZvertexTruthKaonPlus(0),
    fZvertexTruthKaonMinus(0),
    fZvertexContaminationKaon(0),
    fZvertexContaminationKaonPlus(0),
    fZvertexContaminationKaonMinus(0),
    fZvertexRecoKaon(0),
    fZvertexRecoKaonPlus(0),
    fZvertexRecoKaonMinus(0),
    fZvertexTruthProton(0),
    fZvertexTruthProtonPlus(0),
    fZvertexTruthProtonMinus(0),
    fZvertexContaminationProton(0),
    fZvertexContaminationProtonPlus(0),
    fZvertexContaminationProtonMinus(0),
    fZvertexRecoProton(0),
    fZvertexRecoProtonPlus(0),
    fZvertexRecoProtonMinus(0),
    fHistEtaVertexzTruthPion(0),
    fHistEtaVertexzTruthPionPlus(0),
    fHistEtaVertexzTruthPionMinus(0),
    fHistEtaVertexzTruthKaon(0),
    fHistEtaVertexzTruthKaonPlus(0),
    fHistEtaVertexzTruthKaonMinus(0),
    fHistEtaVertexzTruthProton(0),
    fHistEtaVertexzTruthProtonPlus(0),
    fHistEtaVertexzTruthProtonMinus(0),
    fHistEtaVertexzRecoPion(0),
    fHistEtaVertexzRecoPionPlus(0),
    fHistEtaVertexzRecoPionMinus(0),
    fHistEtaVertexzRecoKaon(0),
    fHistEtaVertexzRecoKaonPlus(0),
    fHistEtaVertexzRecoKaonMinus(0),
    fHistEtaVertexzRecoProton(0),
    fHistEtaVertexzRecoProtonPlus(0),
    fHistEtaVertexzRecoProtonMinus(0),
    fHistEtaVertexzContaminationPion(0),
    fHistEtaVertexzContaminationPionPlus(0),
    fHistEtaVertexzContaminationPionMinus(0),
    fHistEtaVertexzContaminationKaon(0),
    fHistEtaVertexzContaminationKaonPlus(0),
    fHistEtaVertexzContaminationKaonMinus(0),
    fHistEtaVertexzContaminationProton(0),
    fHistEtaVertexzContaminationProtonPlus(0),
    fHistEtaVertexzContaminationProtonMinus(0),
    fHistPhiVertexzTruthPion(0),
    fHistPhiVertexzTruthPionPlus(0),
    fHistPhiVertexzTruthPionMinus(0),
    fHistPhiVertexzTruthKaon(0),
    fHistPhiVertexzTruthKaonPlus(0),
    fHistPhiVertexzTruthKaonMinus(0),
    fHistPhiVertexzTruthProton(0),
    fHistPhiVertexzTruthProtonPlus(0),
    fHistPhiVertexzTruthProtonMinus(0),
    fHistPhiVertexzRecoPion(0),
    fHistPhiVertexzRecoPionPlus(0),
    fHistPhiVertexzRecoPionMinus(0),
    fHistPhiVertexzRecoKaon(0),
    fHistPhiVertexzRecoKaonPlus(0),
    fHistPhiVertexzRecoKaonMinus(0),
    fHistPhiVertexzRecoProton(0),
    fHistPhiVertexzRecoProtonPlus(0),
    fHistPhiVertexzRecoProtonMinus(0),
    fHistPhiVertexzContaminationPion(0),
    fHistPhiVertexzContaminationPionPlus(0),
    fHistPhiVertexzContaminationPionMinus(0),
    fHistPhiVertexzContaminationKaon(0),
    fHistPhiVertexzContaminationKaonPlus(0),
    fHistPhiVertexzContaminationKaonMinus(0),
    fHistPhiVertexzContaminationProton(0),
    fHistPhiVertexzContaminationProtonPlus(0),
    fHistPhiVertexzContaminationProtonMinus(0)
{ 
    for(Int_t ipart=0;ipart<6;ipart++){
   fNsigmaTPC[ipart]=999.0;
   fNsigmaTOF[ipart]=999.0;
   }

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


  enum kDetectorPID_ { kTPCTOFpid_, kTogether_ };
  enum kParticleType_ { kPion_, kKaon_, kProton_ };


    void SetParticleType(kParticleType_ particletype_) {
    fParticleType_ = particletype_; }
    void SetDetectorPID(kDetectorPID_ detectorpid_) {
    fDetectorPID_ = detectorpid_; }



// SET TOF Miss Match : 

 void SetMisMatchTOFProb(Double_t fmistmatchTOF,Bool_t ftofMisMatch){
       fMistMatchTOFProb=fmistmatchTOF;
       fTOFMisMatch=ftofMisMatch;
    }

void SetRapidityUse(Bool_t rapidityUse) {fRapidityInsteadOfEta=rapidityUse;}
  
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

   Double_t Beta(AliAODTrack *track);

   Bool_t IsTPCPID(AliAODTrack* track) const; //Check Track ,Is Track TPC?
   Bool_t IsTOFPID(AliAODTrack* track) const; // Check Track , Is Track TOF?

   Double_t GetNsigmas(AliPIDResponse* fPIDResponse , AliAODTrack* track , Int_t specie);


 
 private:
  AliAODEvent* fAOD; //! AOD object  
  TClonesArray *fArrayMC; //! array of MC particles  
  TList       *fQAList; //! QA list
  TList       *fOutputList; //! Output list
 
   Double_t fNSigmaPID;
   Bool_t fRapidityInsteadOfEta;

  // Truth MC Particle Pt , Eta abd Phi distribution Histogram 

// Contamination and Purity Histogram 



   TH3F *Hist3dPionContamination; // Contamination Pion in Eta , Pt , Phi
   TH3F *Hist3dPionPlusContamination; // Contamination Pion + in Eta , Pt , Phi
   TH3F *Hist3dPionMinusContamination; // Contamination Pion - in Eta , Pt , Phi

   TH3F *Hist3dKaonContamination; // Contamination Kaon in Eta , Pt , Phi
   TH3F *Hist3dKaonPlusContamination; // Contamination Kaon + in Eta , Pt , Phi
   TH3F *Hist3dKaonMinusContamination; // Contamination Kaon - in Eta , Pt , Phi

   TH3F *Hist3dProtonContamination; // Contamination Proton in Eta , Pt , Phi
   TH3F *Hist3dProtonPlusContamination; // Contamination Proton + in Eta , Pt , Phi
   TH3F *Hist3dProtonMinusContamination; // Contamination Proton - in Eta , Pt , Phi

// Contamination and Purity Histigram 

// Nsigma Plot for Pion, Kaon and Proton

TH2F *fHistNsigmaTPCPionBeforePIDCut;
TH2F *fHistNsigmaTPCKaonBeforePIDCut;
TH2F *fHistNsigmaTPCProtonBeforePIDCut;

TH2F *fHistNsigmaTOFPionBeforePIDCut;
TH2F *fHistNsigmaTOFKaonBeforePIDCut;
TH2F *fHistNsigmaTOFProtonBeforePIDCut;

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
  TH2F        *fHistNSigmaTPCvsPtbeforeElectronPID;
  TH2F        *fHistNSigmaTPCvsPtafterElectronPID;

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

  TH1F    *fZvertexTruthPion;
  TH1F    *fZvertexTruthPionPlus;
  TH1F    *fZvertexTruthPionMinus;
   
  TH1F    *fZvertexTruthKaon;
  TH1F    *fZvertexTruthKaonPlus;
  TH1F    *fZvertexTruthKaonMinus;

  TH1F    *fZvertexTruthProton;
  TH1F    *fZvertexTruthProtonPlus;
  TH1F    *fZvertexTruthProtonMinus;

// Efficieny in Eta for Different vertex cut; 

  TH2F     *fHistEtaVertexzTruthPion;
  TH2F     *fHistEtaVertexzTruthPionPlus;
  TH2F     *fHistEtaVertexzTruthPionMinus;

  TH2F     *fHistEtaVertexzTruthKaon;
  TH2F     *fHistEtaVertexzTruthKaonPlus;
  TH2F     *fHistEtaVertexzTruthKaonMinus;

  TH2F     *fHistEtaVertexzTruthProton;
  TH2F     *fHistEtaVertexzTruthProtonPlus;
  TH2F     *fHistEtaVertexzTruthProtonMinus;

  TH2F     *fHistEtaVertexzRecoPion;
  TH2F     *fHistEtaVertexzRecoPionPlus;
  TH2F     *fHistEtaVertexzRecoPionMinus;

  TH2F     *fHistEtaVertexzRecoKaon;
  TH2F     *fHistEtaVertexzRecoKaonPlus;
  TH2F     *fHistEtaVertexzRecoKaonMinus;

  TH2F     *fHistEtaVertexzRecoProton;
  TH2F     *fHistEtaVertexzRecoProtonPlus;
  TH2F     *fHistEtaVertexzRecoProtonMinus;

  TH2F     *fHistEtaVertexzContaminationPion;
  TH2F     *fHistEtaVertexzContaminationPionPlus;
  TH2F     *fHistEtaVertexzContaminationPionMinus;

  TH2F     *fHistEtaVertexzContaminationKaon;
  TH2F     *fHistEtaVertexzContaminationKaonPlus;
  TH2F     *fHistEtaVertexzContaminationKaonMinus;

  TH2F     *fHistEtaVertexzContaminationProton;
  TH2F     *fHistEtaVertexzContaminationProtonPlus;
  TH2F     *fHistEtaVertexzContaminationProtonMinus;

  TH2F     *fHistPhiVertexzTruthPion;
  TH2F     *fHistPhiVertexzTruthPionPlus;
  TH2F     *fHistPhiVertexzTruthPionMinus;

  TH2F     *fHistPhiVertexzTruthKaon;
  TH2F     *fHistPhiVertexzTruthKaonPlus;
  TH2F     *fHistPhiVertexzTruthKaonMinus;

  TH2F     *fHistPhiVertexzTruthProton;
  TH2F     *fHistPhiVertexzTruthProtonPlus;
  TH2F     *fHistPhiVertexzTruthProtonMinus;

  TH2F     *fHistPhiVertexzRecoPion;
  TH2F     *fHistPhiVertexzRecoPionPlus;
  TH2F     *fHistPhiVertexzRecoPionMinus;

  TH2F     *fHistPhiVertexzRecoKaon;
  TH2F     *fHistPhiVertexzRecoKaonPlus;
  TH2F     *fHistPhiVertexzRecoKaonMinus;

  TH2F     *fHistPhiVertexzRecoProton;
  TH2F     *fHistPhiVertexzRecoProtonPlus;
  TH2F     *fHistPhiVertexzRecoProtonMinus;

  TH2F     *fHistPhiVertexzContaminationPion;
  TH2F     *fHistPhiVertexzContaminationPionPlus;
  TH2F     *fHistPhiVertexzContaminationPionMinus;

  TH2F     *fHistPhiVertexzContaminationKaon;
  TH2F     *fHistPhiVertexzContaminationKaonPlus;
  TH2F     *fHistPhiVertexzContaminationKaonMinus;

  TH2F     *fHistPhiVertexzContaminationProton;
  TH2F     *fHistPhiVertexzContaminationProtonPlus;
  TH2F     *fHistPhiVertexzContaminationProtonMinus;
// Efficieny in Eta and Phi for Different vertex cut; 

  TH1F    *fZvertexContaminationPion;
  TH1F    *fZvertexContaminationPionPlus;
  TH1F    *fZvertexContaminationPionMinus;
  TH1F    *fZvertexRecoPion;
  TH1F    *fZvertexRecoPionPlus;
  TH1F    *fZvertexRecoPionMinus;

  TH1F    *fZvertexContaminationKaon;
  TH1F    *fZvertexContaminationKaonPlus;
  TH1F    *fZvertexContaminationKaonMinus;

  TH1F    *fZvertexRecoKaon;
  TH1F    *fZvertexRecoKaonPlus;
  TH1F    *fZvertexRecoKaonMinus;

  TH1F    *fZvertexContaminationProton;
  TH1F    *fZvertexContaminationProtonPlus;
  TH1F    *fZvertexContaminationProtonMinus;

  TH1F    *fZvertexRecoProton;
  TH1F    *fZvertexRecoProtonPlus;
  TH1F    *fZvertexRecoProtonMinus; 

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

  Double_t fMistMatchTOFProb;
  Bool_t fTOFMisMatch;


  Double_t fNsigmaTPC[6];
  Double_t fNsigmaTOF[6];
 
  kParticleType_ fParticleType_; // particle type for analysis
  kDetectorPID_ fDetectorPID_; // particle type for analysis


  AliAnalysisTaskEffContPIDBF(const AliAnalysisTaskEffContPIDBF&); // not implemented
  AliAnalysisTaskEffContPIDBF& operator=(const AliAnalysisTaskEffContPIDBF&); // not implemented
  
  ClassDef(AliAnalysisTaskEffContPIDBF, 1); // example of analysis
};

#endif
