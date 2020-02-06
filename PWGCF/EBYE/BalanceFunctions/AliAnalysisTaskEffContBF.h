#ifndef ALIANALYSISTASKEFFCONTBF_H
#define ALIANALYSISTASKEFFCONTBF_H


// ---------------------------------------------------------------------
//
// Task for calculating the efficiency and contamination of the Balance 
// Function for single particles and pairs
// 
// ---------------------------------------------------------------------

class TList;
class TH1F;
class TH3D;
class TH2F;
class TString;
class AliAODEvent;
class AliAODInputHandler;
class TH2D;

#include <AliPID.h>
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEffContBF : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskEffContBF();
  AliAnalysisTaskEffContBF(const char *name);
  virtual ~AliAnalysisTaskEffContBF() {}
  
  enum etriggerSel{kMB, kCentral, kINT7, kppHighMult};
  
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

  void SetRejectInjectedSignalsLabelAboveThreshold() {
    fRejectLabelAboveThreshold=kTRUE;
    fInjectedSignals = kTRUE; 
  }
  
  void SetRejectInjectedSignalsGenName(TString genToBeKept) {
    fGenToBeKept = genToBeKept; 
    fRejectCheckGenName=kTRUE;
    fInjectedSignals = kTRUE;
  }

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

  void SetExcludeElectronsInMC()  {fExcludeElectronsInMC = kTRUE;}
  
  void SetUseParticleID(Bool_t usePID=kFALSE, AliPID::EParticleType partOfInterest = AliPID::kPion) {
    fUsePIDstrategy = usePID;
    fUsePIDFromPDG = usePID;
    fpartOfInterest = partOfInterest;
  }

 void SetUsePIDfromPDG(Bool_t usePID=kFALSE, AliPID::EParticleType partOfInterest = AliPID::kPion) {
   fUsePIDFromPDG = usePID;
   fpartOfInterest = partOfInterest;
   fMassParticleOfInterest = AliPID::ParticleMass(partOfInterest);
  }
  
  void SetUsePIDnSigmaComb(Bool_t UsePIDnSigmaComb){
     fUsePIDnSigmaComb = UsePIDnSigmaComb;
  }
  void SetBayesPIDThr(Double_t Thr){
      fBayesPIDThr = Thr;
  }

  void SetTOFBCPileUpCut(){
    fUseTOFBCPileUpCut = kTRUE;
  }

  void SetTPCInOutRowsCut(Int_t innermostRows = 2,  Int_t outermostRows = 20){
    fUseTPCInOutRowsCut = kTRUE;
    fInRows = innermostRows;
    fOutRows = outermostRows; 
  }
  
  //Track cuts
  void SetMinNumberOfTPCClusters(Double_t min) {
    fMinNumberOfTPCClusters = min;}
  void SetMaxChi2PerTPCCluster(Double_t max) {
    fMaxChi2PerTPCCluster = max;}  

  void SetExtra2DDCACutsAOD(Double_t DCAxy, Double_t DCAz){
    fDCAxyCut  = DCAxy;
    fDCAzCut = DCAz;
  }

  void SetExtraTPCCutsAOD(Double_t maxTPCchi2, Int_t minNClustersTPC, Int_t minNTPCCrossedRows, Float_t minNTPCFindableCls){
    fTPCchi2Cut      = maxTPCchi2;
    fNClustersTPCCut = minNClustersTPC;
    fMinTPCCrossedRows = minNTPCCrossedRows;
    fMinTPCRowsOverFindableCls =  minNTPCFindableCls;
  }
  
  void SetMaxDCAxy(Double_t max) {
    fMaxDCAxy = max;}
  void SetMaxDCAz(Double_t max) {
    fMaxDCAz = max;}
  void SetMinPt(Double_t minPt) {
    fMinPt = minPt;}
  void SetMaxPt(Double_t maxPt) {
    fMaxPt = maxPt;}
    
  void SetUseY(){
    fUseY = kTRUE;
  }  
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


  void SetUseRaaGeoCut(Float_t deadZoneWidth = 3, Float_t cutGeoNcrNclLength = 130, Float_t cutGeoNcrNclGeom1Pt = 1.5, Float_t cutGeoNcrNclFractionNcr = 0.85, Float_t cutGeoNcrNclFractionNcl = 0.7){
    fUseRaaGeoCut=kTRUE;
    fDeadZoneWidth = deadZoneWidth; 
    fCutGeoNcrNclLength = cutGeoNcrNclLength;
    fCutGeoNcrNclGeom1Pt = cutGeoNcrNclGeom1Pt;
    fCutGeoNcrNclFractionNcr = cutGeoNcrNclFractionNcr;
    fCutGeoNcrNclFractionNcl = cutGeoNcrNclFractionNcl;
  }

   
 private:
  AliAODEvent* fAOD; //! AOD object  
  TClonesArray *fArrayMC; //! array of MC particles  
  TList       *fQAList; //! QA list
  TList       *fOutputList; //! Output list
  
  // QA histograms
  TH1F        *fHistEventStats; //!event stats
  TH1F        *fHistCentrality; //!centrality
  TH1F        *fHistNMult; //! nmult  
  TH1F        *fHistVz;//!
  TH2F        *fHistDCA;//DCA z vs DCA xy
  TH2F        *fHistNSigmaTPCvsPtbeforePID;//TPC nsigma vs pT before PID cuts (QA histogram)
  TH2F        *fHistNSigmaTPCvsPtafterPID;//TPC nsigma vs pT after PID cuts (QA histogram)

  // output histograms
  TH3D        *fHistContaminationSecondariesPlus;//!
  TH3D        *fHistContaminationSecondariesMinus;//!

  TH3D        *fHistContaminationSecondariesMaterialPlus;//!
  TH3D        *fHistContaminationSecondariesMaterialMinus; //!
  TH3D        *fHistContaminationSecondariesWeakDecPlus; //!
  TH3D        *fHistContaminationSecondariesWeakDecMinus; //!
  
  TH3D        *fHistContaminationPrimariesPlus;//!
  TH3D        *fHistContaminationPrimariesMinus;//!
  
  // output histograms (single particles)
  TH3D        *fHistGeneratedEtaPtPhiPlus;//!correction map for positives (generated)
  TH3D        *fHistSurvivedEtaPtPhiPlus;//!correction map positives (survived)
 
  TH3D        *fHistGeneratedEtaPtPhiMinus;//!correction map for negatives (generated)
  TH3D        *fHistSurvivedEtaPtPhiMinus;//!correction map negatives (survived)
 
  TH2F        *fHistGeneratedEtaPtPlusControl;//!correction map for positives (generated)
  TH2F        *fHistSurvivedEtaPtPlusControl;//!correction map positives (survived)
 
  TH2F        *fHistGeneratedEtaPtMinusControl;//!correction map for negatives (generated)
  TH2F        *fHistSurvivedEtaPtMinusControl;//!correction map negatives (survived)
 
  // output histograms (pairs)
  TH2F        *fHistGeneratedEtaPtPlusPlus;//!correction map for ++ (generated)
  TH2F        *fHistSurvivedEtaPtPlusPlus;//!correction map ++ (survived)
 
  TH2F        *fHistGeneratedEtaPtMinusMinus;//!correction map for -- (generated)
  TH2F        *fHistSurvivedEtaPtMinusMinus;//!correction map -- (survived)
 
  TH2F        *fHistGeneratedEtaPtPlusMinus;//!correction map for +- (generated)
  TH2F        *fHistSurvivedEtaPtPlusMinus;//!correction map +- (survived)

  TH2F        *fHistGeneratedPhiEtaPlusPlus;//!correction map for ++ (generated)
  TH2F        *fHistSurvivedPhiEtaPlusPlus;//!correction map ++ (survived)
 
  TH2F        *fHistGeneratedPhiEtaMinusMinus;//!correction map for -- (generated)
  TH2F        *fHistSurvivedPhiEtaMinusMinus;//!correction map -- (survived)
 
  TH2F        *fHistGeneratedPhiEtaPlusMinus;//!correction map for +- (generated)
  TH2F        *fHistSurvivedPhiEtaPlusMinus;//!correction map +- (survived)

  // check pdg    
  TH1F        *fHistPdgGen;
  TH1F        *fHistPdgSurv;

  Bool_t  fUseCentrality;// Bool_t use centrality or not
  TString fCentralityEstimator;// "V0M","TRK","TKL","ZDC","FMD"
  Float_t fCentralityPercentileMin; // min centrality percentile 
  Float_t fCentralityPercentileMax; // max centrality percentile

  Bool_t fInjectedSignals;//Flag for using the rejection of injected signals
  Bool_t fRejectLabelAboveThreshold;// 
  TString fGenToBeKept; // name of the generator that should be kept in the analysis (in case of rejection of injected signals)
  Bool_t fRejectCheckGenName; // Flag for using the rejection of injected signals on a track by track base (different cocktails with respect to fInjectedSignals) 
  Bool_t fExcludeElectronsInMC;
  
  AliPIDResponse *fPIDResponse;     //! PID response object
  Bool_t   fElectronRejection;//flag to use electron rejection
  Bool_t   fElectronOnlyRejection;//flag to use electron rejection with exclusive electron PID (no other particle in nsigma range)
  Double_t fElectronRejectionNSigma;//nsigma cut for electron rejection
  Double_t fElectronRejectionMinPt;//minimum pt for electron rejection (default = 0.)
  Double_t fElectronRejectionMaxPt;//maximum pt for electron rejection (default = 1000.)
  
  AliPIDCombined* fPIDCombined;  //! PID combined

  Bool_t  fUsePIDnSigmaComb;//
  Double_t fBayesPIDThr;//
    
  Bool_t  fUseY;//
    
  Bool_t fUsePIDstrategy; // flag to switch on PID
  Bool_t fUsePIDFromPDG; //flag to switch on MC PID (used for PID tracking eff) 
  AliPID::EParticleType fpartOfInterest; //
  Int_t fPDGCodeWanted;//
  Float_t fMassParticleOfInterest;//
    
  Double_t fVxMax;// vxmax
  Double_t fVyMax;// vymax
  Double_t fVzMax;// vzmax
  
  Int_t fAODTrackCutBit;// track cut bit from track selection (only used for AODs)

  Double_t fMinNumberOfTPCClusters;//!
  Double_t fMaxChi2PerTPCCluster;//!
  Double_t fMaxDCAxy, fMaxDCAz;//!
  Double_t fMinPt, fMaxPt;//!
  Double_t fMinEta, fMaxEta;//!
  Double_t fEtaRangeMin;// acceptance cuts 
  Double_t fEtaRangeMax; // acceptance cuts
  Double_t fPtRangeMin;  // acceptance cuts
  Double_t fPtRangeMax;  // acceptance cuts

  Bool_t fUseTOFBCPileUpCut;
  
  Bool_t fUseTPCInOutRowsCut;
  Int_t fInRows;
  Int_t fOutRows; 
  
  Double_t fDCAxyCut;//2D DCA cut
  Double_t fDCAzCut;//2D DCA cut
  Double_t fTPCchi2Cut;//Chi2 Per Cluster TPC cut
  Int_t fNClustersTPCCut;//Minimum number of TPC clusters cut
  Int_t fMinTPCCrossedRows; //Minimum number of TPC crossed rows cut
  Float_t fMinTPCRowsOverFindableCls;  //Minimum number of TPC findable cls

  Int_t fEtaBin;  // acceptance cuts
  Int_t fdEtaBin;  // acceptance cuts
  Int_t fPtBin; // acceptance cuts

  TH3F        *fHistSurvived4EtaPtPhiPlus;//!
  TH3F        *fHistSurvived8EtaPtPhiPlus;//!

  AliESDtrackCuts *fESDtrackCuts; //ESD track cuts

  Bool_t fUseRaaGeoCut; //flag to switch on GeoCut for 2018PbPb data pass1
  Float_t fDeadZoneWidth; //parameters of the cut as implemented in AliESDtrackCuts.h, default values implemented as suggested by DPG and D mesons analysis
  Float_t fCutGeoNcrNclLength;
  Float_t fCutGeoNcrNclGeom1Pt;
  Float_t fCutGeoNcrNclFractionNcr;
  Float_t fCutGeoNcrNclFractionNcl;

  AliAnalysisTaskEffContBF(const AliAnalysisTaskEffContBF&); // not implemented
  AliAnalysisTaskEffContBF& operator=(const AliAnalysisTaskEffContBF&); // not implemented
  
  ClassDef(AliAnalysisTaskEffContBF, 9); // example of analysis
};

#endif
