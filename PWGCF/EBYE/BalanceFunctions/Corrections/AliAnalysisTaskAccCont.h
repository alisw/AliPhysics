#ifndef ALIANALYSISTASKACCCONT_CXX
#define ALIANALYSISTASKACCCONT_CXX

class TList;
class TH1F;
class TH2F;
class TH3F;
class TObjArray;
class AliAODEvent;
class AliAODTrack;
class AliVEvent;

class AliAnalysisUtils;

#include "AliAnalysisTaskSE.h"
#include "TDatabasePDG.h"

class AliAnalysisTaskAccCont : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskAccCont(const char *name = "AliAnalysisTaskAccCont");
  virtual ~AliAnalysisTaskAccCont();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);

  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;
    fVyMax = vy;
    fVzMax = vz;
  }
  void SetAODtrackCutBit(Int_t bit){
    fAODtrackCutBit = bit;
  }

  void SetKinematicsCutsAOD(Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax) {
    fEtaMin  = etamin; fEtaMax  = etamax;
    fPtMin  = ptmin; fPtMax  = ptmax;

  }

  void SetMCRec() {fMCrec = kTRUE;}

  void SetExcludeSecondariesInMC() {fExcludeSecondariesInMCrec = kTRUE;}
 
  void SetExcludeElectronsInMC()  {fExcludeElectronsInMCrec = kTRUE;}  


  void UsePileUpCutsPbPb() {fPbPb = kTRUE;}

  void CheckPileUp() {fCheckPileUp = kTRUE;}

  void SetPileUpCutsParamsLHC15o(Float_t slope, Float_t offset){
    fUseOutOfBunchPileUpCutsLHC15o=kTRUE;
    fPileupLHC15oSlope = slope;
    fPileupLHC15oOffset = offset;
  }

  void SetPileUpCutsJpsigroup(){
    fUseOutOfBunchPileUpCutsLHC15oJpsi=kTRUE;
  }
  
  void UsePileUpCutspPb() {fpPb = kTRUE;}
  
  void USEextendedDCA() {fDCAext = kTRUE;}
  
  void UsePID() {fUsePID = kTRUE;}
  void SetUsePIDnSigmaComb() {fUsePIDnSigmaComb = kTRUE;} //not that if UsePID true and nSigmaComb not activated, Bayesian PID is used!!
  void SetPIDBayesThreshold(Float_t bayesThresh) {fBayesPIDThr = bayesThresh;}
  void SetPIDMomCut(Float_t pidMomCut)  {fPIDMomCut = pidMomCut;} // momentum threshold to move from TPC only and TPC+TOF for both methods: Bayes and nSigma Combined. usually 0.7 for pi and p and o.4 for K.
  void SetUseRapidity() {fUseRapidity = kTRUE;}
  void SetUseNSigmaPIDNewTrial() {
      fUsePIDNewTrial = kTRUE; fUsePIDnSigma = kTRUE;
  }
  
  void SetRejectInjectedSignals() {fExcludeInjectedSignals = kTRUE;}

  void SetRejectInjectedSignalsGenName(TString genToBeKept) {
    fGenToBeKept = genToBeKept; 
    fRejectCheckGenName=kTRUE;
    fExcludeInjectedSignals = kTRUE;
  }
  
  void SetNSigmaPID(Int_t nsigma) {
    fPIDNSigma = nsigma;
  }
 
  void UseOfflineTrigger() {fUseOfflineTrigger = kTRUE;}

  void SetCentralityEstimator(const char* centralityEstimator) {
    fCentralityEstimator = centralityEstimator;}
  void SetCentralityPercentileRange(Double_t min, Double_t max) { 
    fCentralityPercentileMin=min;
    fCentralityPercentileMax=max;
  }

  // enum kParticleOfInterest { kMuon, kElectron, kPion, kKaon, kProton };
  enum kCentralityBinning { kFull, kBins, kMCgen };
  enum kSystem { kPbPb, kpPb };

  void setParticleType(AliPID::EParticleType ptype){
  fParticleOfInterest = ptype;
  
  if(fParticleOfInterest == AliPID::kElectron){
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(11)->Mass();
  }  
 
  else if(fParticleOfInterest == AliPID::kMuon){
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(13)->Mass();
  }
  else if(fParticleOfInterest == AliPID::kPion){
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  }
  else if(fParticleOfInterest == AliPID::kKaon){
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  }
  else if(fParticleOfInterest == AliPID::kProton){
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  }
  else{
    AliWarning("Particle type not known, set fMassParticleOfInterest to pion mass.");
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  }  
  }

 private:
  AliVEvent* gAOD;
  TList *fListQA; //fList object (QA)
  TList *fListResults; //fList object (Results)

  TH2F *fHistEventStats; //event stats
  TH2F *fHistTrackStats; //Track filter bit stats
  TH2F *fHistVx; //x coordinate of the primary vertex
  TH2F *fHistVy; //y coordinate of the primary vertex
  TH2F *fHistVz; //z coordinate of the primary vertex

  TString fCentralityEstimator;      //"V0M","TRK","TKL","ZDC","FMD"
  Double_t fCentralityPercentileMin;//centrality percentile min
  Double_t fCentralityPercentileMax;//centrality percentile max

  TH2F *fHistMultiplicity; //multiplicity of accepted tracks vs centrality
  TH2F *fHistPtCen; //pT spectra vs centrality
  TH2F *fHistPhiCen;
  TH2F *fHistEtaCen;
  TH2F *fHistDCAToVertex2D;
  
  TH1F *fHistPt;
  TH1F *fHistPtbin;
  TH1F *fHistCent;
  TH1F *fHistCentbin;
  TH1F *fHistPhi;
  TH1F *fHistEta;
  TH1F *fHistNClustersTPC;
  TH1F *fHistChi2PerClusterTPC;
  TH1F *fHistDCAToVertexZ;
  TH1F *fHistDCAToVertexXY;
  TH1F *fHistPdg;
     
  TH3F *fHistEtaPhiCent;
  TH3F *fHistPtEtaCent;
  TH3F *fHistPtPhiCent;

  TH3D *fHistEtaPhiVertexPlus;
  TH3D *fHistEtaPhiVertexMinus;   
  TH3D *fHistYPhiVertexPlus;
  TH3D *fHistYPhiVertexMinus;
  TH3F *fHistDCAXYptchargedminus;
  TH3F *fHistDCAXYptchargedplus;
  TH3F *fHistDCAXYptchargedminus_ext;
  TH3F *fHistDCAXYptchargedplus_ext;
    
  TH2D *fHistdEdxVsPTPCbeforePID;
  TH2D *fHistBetavsPTOFbeforePID;
  TH2D *fHistProbTPCvsPtbeforePID;
  TH2D *fHistProbTPCTOFvsPtbeforePID;
  TH2D *fHistNSigmaTPCvsPtbeforePID;
  TH2D *fHistNSigmaTOFvsPtbeforePID;
  TH2D *fHistBetaVsdEdXbeforePID;
  TH2D *fHistNSigmaTPCTOFvsPtbeforePID;
  TH3D *fHistNSigmaTPCTOFPbefPID;
  TH2D *fHistBetavsPTOFafterPID;
  TH2D *fHistdEdxVsPTPCafterPID;
  TH2D *fHistBetaVsdEdXafterPID;
  TH2D *fHistNSigmaTOFvsPtafterPID;
  TH2D *fHistNSigmaTPCvsPtafterPID;
  TH2D *fHistNSigmaTPCTOFvsPtafterPID;
  TH3D *fHistNSigmaTPCTOFPafterPID;

  TH2F *fHistGlobalvsESDBeforePileUpCuts;
  TH2F *fHistGlobalvsESDAfterPileUpCuts;

  TH2F *fHistV0MvsTPCoutBeforePileUpCuts; //histos to monitor pile up cuts J/psi
  TH2F *fHistV0MvsTPCoutAfterPileUpCuts;
 
  TH3F* hNSigmaCutApplied;
  TH3F* hBayesProbab;
  
  Bool_t fUseOfflineTrigger;//Usage of the offline trigger selection
  Bool_t fPbPb;
  Bool_t fpPb;
  Bool_t fCheckPileUp;
  Bool_t fMCrec;
  Bool_t fExcludeSecondariesInMCrec;
  Bool_t fExcludeElectronsInMCrec;
  Bool_t fExcludeInjectedSignals; //Flag to reject MC injected signals from MC analysis
  Bool_t fRejectCheckGenName; // Flag to activate the injected signal rejection based on the name of the MC generator 
  TString fGenToBeKept; //String to select the generator name that has to be kept for analysis
  

  TClonesArray* fArrayMC;

  Float_t fPileupLHC15oSlope; //parameters for LHC15o pile-up rejection  default: slope=3.35, offset 15000
  Float_t fPileupLHC15oOffset;
  Bool_t fUseOutOfBunchPileUpCutsLHC15o;//usage of correlation cuts to exclude out of bunche pile up. To be used for 2015 PbPb data.

  Bool_t fUseOutOfBunchPileUpCutsLHC15oJpsi;//
  
  Bool_t fUsePID;
  Bool_t fDCAext;
  Bool_t fUseRapidity;
  Bool_t fUsePIDnSigmaComb;
  Bool_t fUsePIDNewTrial;
  Bool_t fUsePIDnSigma;
 
  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax
  
  Int_t fAODtrackCutBit;//AOD track cut bit from track selection
  
  Int_t fPIDNSigma;
  Double_t fMassParticleOfInterest;
  AliPID::EParticleType fParticleOfInterest;
 
  Double_t fEtaMin; 
  Double_t fEtaMax; 
  Double_t fPtMin; 
  Double_t fPtMax;

  Float_t fBayesPIDThr;
  Float_t fPIDMomCut;
  
  //AliAnalysisUtils
  AliAnalysisUtils *fUtils;//AliAnalysisUtils
  AliPIDResponse *fPIDResponse;
  AliPIDCombined *fPIDCombined;     //! combined PID object
 
  AliAnalysisTaskAccCont(const AliAnalysisTaskAccCont&); // not implemented
  AliAnalysisTaskAccCont& operator=(const AliAnalysisTaskAccCont&); // not implemented
  
  ClassDef(AliAnalysisTaskAccCont, 4); // example of analysis
};

#endif
