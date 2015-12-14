#ifndef ALIANALYSISTASKPROTONABSORBTION_H
#define ALIANALYSISTASKPROTONABSORBTION_H


/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskProtonAbsorbtion class
//            This task is for absorption efficiency of protons and antiprotons from ESD only
//-----------------------------------------------------------------
#include "TObject.h"
#include "TH1I.h"
//#include "AliCFContainer.h"
//#include "AliCFManager.h"

class TString;
class TList;
class TH1F;
class TH2F;
class TH3F;
class TF1;

#include "AliAnalysisTaskSE.h"

class AliESDEvent;
class AliESDVertex;
class AliAODEvent;
class AliPhysicsSelection;
class AliTPCPIDResponse; 
class AliPIDResponse;
class AliCFContainer;
class AliCFManager;


class AliAnalysisTaskProtonAbsorbtion : public AliAnalysisTaskSE {

 public:
 enum {
    kStepGenerated       = 0,
    kStepReconstructible = 1,
    kStepReconstructed   = 2,
    kStepSurvived        = 3,
    kNSteps = 4
  };

  enum PIDMode { kRatio = 0, kSigma};

  AliAnalysisTaskProtonAbsorbtion();
  AliAnalysisTaskProtonAbsorbtion(const char *name);
 ~AliAnalysisTaskProtonAbsorbtion();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);


  void SetPhysicsSelection(AliPhysicsSelection* physicsSelection) { fPhysicsSelection = physicsSelection; }
  AliPhysicsSelection* GetPhysicsSelection() const { return fPhysicsSelection; }
  void   SetUsePhysicsSelection(Bool_t usePhysicsSelection = 0) {fUsePhysicsSelection = usePhysicsSelection;}
  void	 SetCentralityWindow(Int_t minCent, Int_t maxCent) {MINCent = minCent, MAXCent = maxCent;}
  void SetDebugMode() {fDebugMode = kTRUE;}

  void   SetPtSpace(Int_t nBinsPt, Double_t gMinPt, Double_t gMaxPt){
	nbinsPtPrim = nBinsPt; fLowPtPrim = gMinPt; fHighPtPrim = gMaxPt;
  }

  void   SetCollidingSystems(Short_t collidingSystems = 0) {fCollidingSystems = collidingSystems;}
  void   SetMultiplicityMode(Bool_t flagMulti = kTRUE) {fMultiplicityMode = flagMulti;}
  void   SetMaxPrimaryVtxPosZ(const Float_t maxPrimaryVtxPosZ = 100.) {fMaxPrimaryVtxPosZ = maxPrimaryVtxPosZ;}
  void   SetMinTPCClusters(Int_t minTPCClusters) {fMinTPCClusters = minTPCClusters;}
  void   SetMinITSClusters(Int_t minITSClusters) {fMinITSClusters = minITSClusters;}
  void   SetMaxChi2PerTPCCluster(Double_t maxChi2PerTPCCluster) {fMaxChi2PerTPCCluster = maxChi2PerTPCCluster;}
  void   SetMaxChi2PerITSCluster(Double_t maxChi2PerITSCluster) {fMaxChi2PerITSCluster = maxChi2PerITSCluster;}
  void   SetMaxDCAXY(Double_t maxDCAXY) {
    fMaxDCAXY = maxDCAXY;
    fMaxDCAXYFlag = kTRUE;
  }

 void    SetMaxDCAZ(Double_t maxDCAZ) {
    fMaxDCAZ = maxDCAZ;
    fMaxDCAZFlag = kTRUE;
      }



//PID related functions
  void SetPIDMode(PIDMode pidmode, Double_t cut1, Double_t cut2, Double_t bound) 
  {fPIDMode = pidmode;
  if(fPIDMode == kSigma) { fNSigma1 = cut1; fNSigma2 = cut2; fNBoundP = bound;}
  if(fPIDMode == kRatio) { fNRatio1 = cut1; fNRatio2 = cut2; fNBoundP = bound;}
  }
  Bool_t IsProton (AliESDtrack *track);
  Bool_t IsAccepted (AliESDtrack *track);
  Bool_t IsPrimary (AliESDEvent *esd, const AliESDVertex *vertex, AliESDtrack* track);

Double_t Rapidity(Double_t,Double_t,Double_t,Int_t)const;

private:
  TString fOADBPath;                   // OADB path to use
  AliPIDResponse *fPIDResponse;        //! PID response Handler


  TString      fAnalysisType;                   //  ESD or AOD
  Short_t      fCollidingSystems;               //  Colliding systems 0/1 for pp/PbPb
  Bool_t       fUsePhysicsSelection;            //  Delegate event selection to AliPhysicsSelectionTask
  Float_t      fMaxPrimaryVtxPosZ;              //  Primary vertex selection in Z
  Bool_t       IsLabelUsed(TArrayI array, Int_t label);

AliPhysicsSelection* fPhysicsSelection; // event selection class

  Bool_t       fMultiplicityMode; 

  TList       *fListHist;                       //! List of histograms
  TList	      *fList;

  AliCFContainer *containerProtonsPrim;
  AliCFContainer *containerAntiProtonsPrim;

  AliCFContainer *containerProtonsPrimMulti;
  AliCFContainer *containerAntiProtonsPrimMulti;

  //Analysis containers
  AliCFManager   *fCFManagerProtonsPrim;      // CF manager protons
  AliCFManager   *fCFManagerAntiProtonsPrim;  // CF manager antiprotons

  AliCFManager   *fCFManagerProtonsPrimMulti;      // CF manager protons
  AliCFManager   *fCFManagerAntiProtonsPrimMulti;  // CF manager antiprotons

  Int_t	       fMinTPCClusters,fMinITSClusters;
  Double_t     fMaxChi2PerTPCCluster;
  Double_t     fMaxChi2PerITSCluster;
  Double_t     fMaxDCAXY, fMaxDCAZ;
  Bool_t fMaxDCAXYFlag, fMaxDCAZFlag;
  Int_t	       MINCent,MAXCent;

  PIDMode fPIDMode; 				//PID mode:dE/dx Ratio-Nsigma areas
  Double_t fNBoundP;				//Momentum bound between two pid cuts
  Double_t fNSigma1; 				//N-sigma cut in the dE/dx band
  Double_t fNSigma2; 				//N-sigma cut in the dE/dx band
  Double_t fNRatio1; 				//min value of the ratio of the measured dE/dx vs the expected
  Double_t fNRatio2; 				//min value of the ratio of the measured dE/dx vs the expected	


  Int_t mintrackrefsTPC;
  Int_t nbinsYMulti;
  Double_t fLowYMulti;
  Double_t fHighYMulti;
  Int_t nbinsYPrim;
  Double_t fLowYPrim;
  Double_t fHighYPrim;
  Int_t nbinsPtPrim;
  Double_t fLowPtPrim;
  Double_t fHighPtPrim;

  //Debug
  Bool_t fDebugMode; //Enable the debug mode

  AliAnalysisTaskProtonAbsorbtion(const AliAnalysisTaskProtonAbsorbtion&);            // not implemented
  AliAnalysisTaskProtonAbsorbtion& operator=(const AliAnalysisTaskProtonAbsorbtion&); // not implemented
  
  ClassDef(AliAnalysisTaskProtonAbsorbtion, 1);
};

#endif
