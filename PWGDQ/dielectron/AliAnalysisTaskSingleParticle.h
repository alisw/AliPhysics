#ifndef ALIANALYSISTASKSINGLEPARTICLE_H
#define ALIANALYSISTASKSINGLEPARTICLE_H

//===============================================================
//
// Analysis task for constructing MC or data driven single particle efficiencies
//
// Ionut C. Arsene, EMMI-GSI, 2011/12/07 
//
//===============================================================

#include "TList.h"

#include "AliAnalysisTaskSE.h"

class TH1D;
class TList;
class AliAnalysisCuts;
class AliCFContainer;
class AliVEvent;
class AliDielectronHistos;
class AliESDv0Cuts;
class AliESDv0KineCuts;

class AliAnalysisTaskSingleParticle : public AliAnalysisTaskSE {
  
public:
  enum { kNMaxDimensions = 20 };
    
  AliAnalysisTaskSingleParticle();
  AliAnalysisTaskSingleParticle(const char *name);
  virtual ~AliAnalysisTaskSingleParticle();

  virtual void UserExec(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
    
  void SetFillTRDfriendPH(Bool_t fill=kTRUE) {fFillTRDfriendPH = fill;}
  void UsePhysicsSelection(Bool_t phy=kTRUE) {fSelectPhysics=phy;}
  void SetTriggerMask(UInt_t mask) {fTriggerMask=mask;}
  UInt_t GetTriggerMask() const { return fTriggerMask; }
  void SetRejectPileup(Bool_t pileup=kTRUE)     { fRejectPileup=pileup;     }
  
  void SetEventFilter(AliAnalysisCuts * const filter) {fEventFilter=filter;}
  void SetTrackFilter(AliAnalysisCuts * const filter) {fTrackFilter=filter;}
  void SetPairFilter(AliAnalysisCuts* const filter) {fPairFilter=filter;}
  void SetV0Cuts(AliESDv0Cuts* const cuts) {fV0Cuts = cuts;}
  void SetLambdaFilter(AliAnalysisCuts* const filter) {fLambdaFilter = filter;}
  void SetK0sFilter(AliAnalysisCuts* const filter) {fK0sFilter = filter;}
  void SetHistogramManager(AliDielectronHistos * const histos) { fHistos=histos; }
  void SetV0KineCuts(AliESDv0KineCuts* const v0cuts) {fV0KineCuts = v0cuts;}

  void AddCFVar(Int_t var, Int_t nbins, Double_t lowLim, Double_t highLim);
  void AddCFVar(Int_t var, const Char_t* bins);

protected:
//  enum {kAllEvents=0, kPhysicsSelectionEvents, kFilteredEvents, kEventStatBins};
//  enum {kEventVtxZ=0, kNTracklets, kCentrality, kPt, kPin, kPhi, kEta, kDeltaPt, kDeltaPhi, kDeltaEta, kTPCnCls, kTPCnSigEle, kTPCnSigPio, kTPCnSigPro, kNVariables}; 

  AliCFContainer*  fCfContainer;      //  CF container
  AliDielectronHistos *fHistos;       // Histogram manager
  TList            fHistogramList;    // histogram list from the manager

  Bool_t           fSelectPhysics;    // Whether to use physics selection
  UInt_t           fTriggerMask;      // Event trigger mask
  Bool_t           fRejectPileup;     // pileup rejection wanted
  Bool_t           fFillTRDfriendPH;  // use the task to fill a TRD tracklet PH container

  AliAnalysisCuts* fEventFilter;      // event filter
  AliAnalysisCuts* fTrackFilter;      // tracking filter
  AliAnalysisCuts* fPairFilter;       // pair filter
  AliESDv0Cuts*    fV0Cuts;           // v0 standard filter
  AliAnalysisCuts* fLambdaFilter;     // additional cuts for lambda v0 exclusion
  AliAnalysisCuts* fK0sFilter;        // additional cuts for K0s v0 exclusion
  AliESDv0KineCuts* fV0KineCuts;      // V0 kine cuts
  
  Int_t            fCFNVars;          // number of CF variables
  Int_t            fCFVarsEnabled[kNMaxDimensions];    // id for every dimension
  Int_t            fCFVarsNbins[kNMaxDimensions];      // number of bins for every CF dimension
  Double_t         fCFVarRanges[kNMaxDimensions][2];   // range for CF dimensions
  TString          fCFVarBins[kNMaxDimensions];        // bin limits for CF dimensions

  TH1D *fEventStat;                  //! Histogram with event statistics
  
  void FillContainer(Int_t step, const Double_t* values, Int_t trackId, Int_t pdg);
  
  AliAnalysisTaskSingleParticle(const AliAnalysisTaskSingleParticle &c);
  AliAnalysisTaskSingleParticle& operator= (const AliAnalysisTaskSingleParticle &c);
  
  ClassDef(AliAnalysisTaskSingleParticle, 2); //Analysis Task handling single particle cuts
};
#endif
