#ifndef ALIANALYSISTASKPROTON_H
#define ALIANALYSISTASKPROTON_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskProton class
//            This task is for raw antiproton to proton ratio from ESD only
//-----------------------------------------------------------------
class TString;
class TList;
class TH1F;
class TH2F;
class TH3F;
class TF1;

#include "AliAnalysisTaskSE.h"

#include "AliCFContainer.h"
class AliCFDataGrid;
class AliESDEvent;
class AliESDVertex;
class AliAODEvent;
class AliPhysicsSelection;
class AliTPCPIDResponse; 
class AliPIDResponse;


class AliAnalysisTaskProton : public AliAnalysisTaskSE {

 public:
  enum PIDMode { kRatio = 0, kSigma};

  AliAnalysisTaskProton();
  AliAnalysisTaskProton(const char *name);
 ~AliAnalysisTaskProton();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetPhysicsSelection(AliPhysicsSelection* physicsSelection) { fPhysicsSelection = physicsSelection; }
  AliPhysicsSelection* GetPhysicsSelection() const { return fPhysicsSelection; }
  void   SetUsePhysicsSelection(Bool_t usePhysicsSelection = 0) {fUsePhysicsSelection = usePhysicsSelection;} 
  void SetDebugMode() {fDebugMode = kTRUE;}

  void   SetCollidingSystems(Short_t collidingSystems = 0) {fCollidingSystems = collidingSystems;}
  void   SetMultiplicityMode(Bool_t flagMulti = kTRUE) {fMultiplicityMode = flagMulti;}
  void   SetPtSpace(Int_t nBinsPt, Double_t gMinPt, Double_t gMaxPt){
	nbinsPt = nBinsPt; fLowPt = gMinPt; fHighPt = gMaxPt;
  }
  
  void	 SetCentralityWindow(Int_t minCent, Int_t maxCent) {MINCent = minCent, MAXCent = maxCent;}
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

  void SetPtDependentDCAxy(Int_t nSigma, Double_t p0, 
			   Double_t p1, Double_t p2);


  Bool_t  IsUsedPtDependentDCAxy() const {return fPtDependentDcaXYFlag;}




 private:
  Int_t   EventNo;   

  TString fOADBPath;                   // OADB path to use
  AliPIDResponse *fPIDResponse;        //! PID response Handler
  Int_t   fRun;                        //! current run number
  Int_t   fOldRun;                     //! current run number
  Int_t   fRecoPass;                   //! reconstruction pass
  void SetRecoInfo();

  TString      fAnalysisType;                   //  ESD or AOD
  Short_t      fCollidingSystems;               //  Colliding systems 0/1 for pp/PbPb
  Bool_t       fUsePhysicsSelection;            //  Delegate event selection to AliPhysicsSelectionTask
  Float_t      fMaxPrimaryVtxPosZ;              //  Primary vertex selection in Z
  Bool_t       IsLabelUsed(TArrayI array, Int_t label);

  //Analysis containers
  AliCFContainer *fProtonContainer; //container for protons
  AliCFContainer *fAntiProtonContainer; //container for antiprotons

  TH1F   *fHistEventStats; //event statistics



	TList		*fListAnalysis;                       //! List of histograms
	TList		*fGlobalQAList;
	TList		*fListQA;
	TList		*fQA2DList;
//	TList		*fListSystematics;
//	TList		*fSysProtonLow;
//	TList		*fSysAntiProtonLow;
//	TList		*fSysProtonHigh;
//	TList		*fSysAntiProtonHigh;
	TH1D 		*fHistMultiplicity;//Multiplicity

	TH2F		*gHistdEdxP;
	TH2F		*gHistProtonsdEdxP;
	TH3F		*gHistProtonsDCAxyEtaPt;
	TH3F		*gHistAntiProtonsDCAxyEtaPt;
/*	TH3F		*gHistFAKEProtonsDCAxyEtaPt;
	TH3F		*gHistFAKEAntiProtonsDCAxyEtaPt;
	TH3F		*gHistFAKESPDProtonsDCAxyEtaPt;
	TH3F		*gHistFAKESPDAntiProtonsDCAxyEtaPt;*/

	TH3F		*gHistProtonsDCAzEtaPt;
	TH3F		*gHistAntiProtonsDCAzEtaPt;
	TH3F		*gHistProtonsDCAzCentPt;
	TH3F		*gHistAntiProtonsDCAzCentPt;
/*	TH3F		*gHistFAKEProtonsDCAzEtaPt;
	TH3F		*gHistFAKEAntiProtonsDCAzEtaPt;
	TH3F		*gHistFAKEProtonsDCAzCentPt;
	TH3F		*gHistFAKEAntiProtonsDCAzCentPt;
	TH3F		*gHistFAKESPDProtonsDCAzEtaPt;
	TH3F		*gHistFAKESPDAntiProtonsDCAzEtaPt;
	TH3F		*gHistFAKESPDProtonsDCAzCentPt;
	TH3F		*gHistFAKESPDAntiProtonsDCAzCentPt;*/


	TH3F		*gHistFieldProtonsEtaPt;
	TH3F		*gHistFieldAntiProtonsEtaPt;
	
	TH3F		*gHistFieldProtonsCentPt;
	TH3F		*gHistFieldAntiProtonsCentPt;

	TH3F		*gHistFieldProtonsLengthPt;
	TH3F		*gHistFieldAntiProtonsLengthPt;

	TH3F		*gHistProtonsLengthCentPt;
	TH3F		*gHistAntiProtonsLengthCentPt;

		



  AliPhysicsSelection *fPhysicsSelection; //Trigger selection: offline

  Bool_t       fMultiplicityMode;

  Int_t	       nbinsY;
  Int_t	       nbinsPt;
  Double_t     fLowY,fHighY,fLowPt,fHighPt;
  Int_t	       fMinTPCClusters,fMinITSClusters;
  Double_t     fMaxChi2PerTPCCluster;
  Double_t     fMaxChi2PerITSCluster;
  Int_t	       MINCent,MAXCent;

  Bool_t fPtDependentDcaXYFlag; //shows if this cut is used or not
  Double_t fMaxDCAXY, fMaxDCAZ, fMaxDCAXYTPC; //max DCA xy
  Bool_t fMaxDCAXYFlag, fMaxDCAZFlag, fMaxDCAXYTPCFlag; //shows if this cut is used or not
  TF1  *fPtDependentDcaXY; //pt dependence dca cut (xy)
  Int_t fNSigmaDCAXY; //n-sigma dca xy cut (pt dependent)

//  AliTPCPIDResponse fTPCpid;                   // Tool data member to manage the TPC Bethe-Bloch info
  PIDMode fPIDMode; 				//PID mode:dE/dx Ratio-Nsigma areas
  Double_t fNBoundP;				//Momentum bound between two pid cuts
  Double_t fNSigma1; 				//N-sigma cut in the dE/dx band
  Double_t fNSigma2; 				//N-sigma cut in the dE/dx band
  Double_t fNRatio1; 				//min value of the ratio of the measured dE/dx vs the expected
  Double_t fNRatio2; 				//min value of the ratio of the measured dE/dx vs the expected

  Bool_t fIsMC;
  Int_t  fUserDataRecoPass;            // forced DATA reco pass

  AliAnalysisTaskProton(const AliAnalysisTaskProton&);            // not implemented
  AliAnalysisTaskProton& operator=(const AliAnalysisTaskProton&); // not implemented

  //Debug
  Bool_t fDebugMode; //Enable the debug mode
  
  ClassDef(AliAnalysisTaskProton, 1);
};

#endif
