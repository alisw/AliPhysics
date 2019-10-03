#ifndef ALIANALYSISTASKPROTONQA_H
#define ALIANALYSISTASKPROTONQA_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskProtonQA class
//            This task is for QAing the protons and antiprotons from ESD only
//-----------------------------------------------------------------
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


class AliAnalysisTaskProtonQA : public AliAnalysisTaskSE {

 public:
  enum PIDMode { kRatio = 0, kSigma};

  AliAnalysisTaskProtonQA();
  AliAnalysisTaskProtonQA(const char *name);
 ~AliAnalysisTaskProtonQA();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);


  void SetPhysicsSelection(AliPhysicsSelection* physicsSelection) { fPhysicsSelection = physicsSelection; }
  AliPhysicsSelection* GetPhysicsSelection() const { return fPhysicsSelection; }
  void   SetUsePhysicsSelection(Bool_t usePhysicsSelection = 0) {fUsePhysicsSelection = usePhysicsSelection;}
  void SetDebugMode() {fDebugMode = kTRUE;}

  void   SetPtSpace(Int_t nBinsPt, Double_t gMinPt, Double_t gMaxPt){
	fNBinsPt = nBinsPt; fMinPt = gMinPt; fMaxPt = gMaxPt;
  }

  void	 SetCentralityWindow(Int_t minCent, Int_t maxCent) {MINCent = minCent, MAXCent = maxCent;}
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
  TList	      *fListQA;

  TH3F*	       gHistPrimaryProtonsDCAxyEtaPt;
  TH3F*	       gHistPrimaryAntiProtonsDCAxyEtaPt;
/*  TH3F*	       gHistFake1PrimaryProtonsDCAxyEtaPt;
  TH3F*	       gHistFake1PrimaryAntiProtonsDCAxyEtaPt;
  TH3F*	       gHistFake2PrimaryProtonsDCAxyEtaPt;
  TH3F*	       gHistFake2PrimaryAntiProtonsDCAxyEtaPt;*/


  TH3F*	       gHistSecondaryProtonsFromWeakDCAxyEtaPt;
  TH3F*	       gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt;
/*  TH3F*	       gHistFake1SecondaryProtonsFromWeakDCAxyEtaPt;
  TH3F*	       gHistFake1SecondaryAntiProtonsFromWeakDCAxyEtaPt;
  TH3F*	       gHistFake2SecondaryProtonsFromWeakDCAxyEtaPt;
  TH3F*	       gHistFake2SecondaryAntiProtonsFromWeakDCAxyEtaPt;*/


  TH3F*	       gHistSecondaryProtonsFromHadronicDCAxyEtaPt;
  TH3F*	       gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt;
/*  TH3F*	       gHistSecondaryFake1ProtonsFromHadronicDCAxyEtaPt;
  TH3F*	       gHistSecondaryFake1AntiProtonsFromHadronicDCAxyEtaPt;
  TH3F*	       gHistSecondaryFake2ProtonsFromHadronicDCAxyEtaPt;
  TH3F*	       gHistSecondaryFake2AntiProtonsFromHadronicDCAxyEtaPt;*/
  TH2F*	       gHistdEdxP;
  TH2F*	       gHistProtonsdEdxP;

  TH3F*	       gHistPrimaryProtonsDCAzEtaPt;
  TH3F*	       gHistPrimaryAntiProtonsDCAzEtaPt;
  TH3F*	       gHistSecondaryProtonsFromWeakDCAzEtaPt;
  TH3F*	       gHistSecondaryAntiProtonsFromWeakDCAzEtaPt;
  TH3F*	       gHistSecondaryProtonsFromHadronicDCAzEtaPt;
  TH3F*	       gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt;

  TH3F*	       gHistPrimaryProtonsDCAzCentPt;
  TH3F*	       gHistPrimaryAntiProtonsDCAzCentPt;
  TH3F*	       gHistSecondaryProtonsFromWeakDCAzCentPt;
  TH3F*	       gHistSecondaryAntiProtonsFromWeakDCAzCentPt;
  TH3F*	       gHistSecondaryProtonsFromHadronicDCAzCentPt;
  TH3F*	       gHistSecondaryAntiProtonsFromHadronicDCAzCentPt;
  

  Int_t	       fNBinsY;
  Int_t	       fNBinsPt;
  Double_t     fMinY,fMaxY,fMinPt,fMaxPt;
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

  //Debug
  Bool_t fDebugMode; //Enable the debug mode

  AliAnalysisTaskProtonQA(const AliAnalysisTaskProtonQA&);            // not implemented
  AliAnalysisTaskProtonQA& operator=(const AliAnalysisTaskProtonQA&); // not implemented
  
  ClassDef(AliAnalysisTaskProtonQA, 1);
};

#endif
