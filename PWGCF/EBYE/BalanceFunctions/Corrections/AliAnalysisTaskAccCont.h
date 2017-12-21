#ifndef ALIANALYSISTASKACCCONT_CXX
#define ALIANALYSISTASKACCCONT_CXX

class TList;
class TH1F;
class TH2F;
class TH3F;
class TObjArray;

class AliAODEvent;
class AliAODTrack;

class AliAnalysisUtils;

#include "AliAnalysisTaskSE.h"

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

  void UsePileUpCuts() {fUseOutOfBunchPileUpCutsLHC15o = kTRUE;}
  
  void USEextendedDCA() {fDCAext = kTRUE;}
  
  void UsePID() {fUsePID = kTRUE;}

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

  enum kParticleOfInterest { kMuon, kElectron, kPion, kKaon, kProton };
  
  void setParticleType(kParticleOfInterest ptype){
  fParticleOfInterest = ptype;
  }

 private:
  AliAODEvent* gAOD;
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
     
  TH3F *fHistEtaPhiCent;
  TH3F *fHistPtEtaCent;
  TH3F *fHistPtPhiCent;

  TH3D *fHistEtaPhiVertexPlus;
  TH3D *fHistEtaPhiVertexMinus;   
  TH3F *fHistDCAXYptchargedminus;
  TH3F *fHistDCAXYptchargedplus;
  TH3F *fHistDCAXYptchargedminus_ext;
  TH3F *fHistDCAXYptchargedplus_ext;
  
  Bool_t fUseOfflineTrigger;//Usage of the offline trigger selection
  Bool_t fUseOutOfBunchPileUpCutsLHC15o;
  Bool_t fUsePID;
  Bool_t fDCAext;

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax

  Int_t fAODtrackCutBit;//AOD track cut bit from track selection
  
  Int_t fPIDNSigma;
  kParticleOfInterest fParticleOfInterest;
 
  Double_t fEtaMin; 
  Double_t fEtaMax; 
  Double_t fPtMin; 
  Double_t fPtMax;

  //AliAnalysisUtils
  AliAnalysisUtils *fUtils;//AliAnalysisUtils
  AliPIDResponse *fPIDResponse;  
 
  AliAnalysisTaskAccCont(const AliAnalysisTaskAccCont&); // not implemented
  AliAnalysisTaskAccCont& operator=(const AliAnalysisTaskAccCont&); // not implemented
  
  ClassDef(AliAnalysisTaskAccCont, 1); // example of analysis
};

#endif
