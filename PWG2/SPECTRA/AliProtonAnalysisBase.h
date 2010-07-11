#ifndef ALIPROTONANALYSISBASE_H
#define ALIPROTONANALYSISBASE_H

/*  See cxx source for full Copyright notice */


/* $Id: AliProtonAnalysisBase.h 31056 2009-02-16 14:31:41Z pchrist $ */

//-------------------------------------------------------------------------
//                       Class AliProtonAnalysisBase
//   This is the base class for the baryon (proton) analysis
//
//    Origin: Panos Christakoglou | Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include "TObject.h"
#include "TString.h"
class TF1;
class TCanvas;
class TList;

#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliPID.h"
class AliESDEvent;
class AliESDtrack;
class AliESDVertex;

class AliProtonAnalysisBase : public TObject {
 public:
  enum TriggerMode { kMB1 = 0, kMB2, kSPDFASTOR };
  enum AnalysisMode { kInvalid = -1, kTPC = 0, kHybrid, kFullHybrid, kGlobal };
  enum PIDMode { kBayesian = 0, kRatio, kSigma};

  AliProtonAnalysisBase();
  virtual ~AliProtonAnalysisBase();

  void SetAnalysisLevel(const char* type) {fProtonAnalysisLevel = type;}
  void SetAnalysisMode(AnalysisMode analysismode) {fProtonAnalysisMode = analysismode;}
  void SetEtaMode() {fAnalysisEtaMode = kTRUE;}
  void SetTriggerMode(TriggerMode triggermode) {
    fAnalysisMC = kTRUE; fTriggerMode = triggermode;}
  void SetPIDMode(PIDMode pidmode) {fProtonPIDMode = pidmode;}

  const char *GetAnalysisLevel() {return fProtonAnalysisLevel.Data();}
  AnalysisMode GetAnalysisMode() const {return fProtonAnalysisMode;}
  Bool_t GetEtaMode() const {return fAnalysisEtaMode;}
  TriggerMode GetTriggerMode() const {return fTriggerMode;}
  PIDMode GetPIDMode() const {return fProtonPIDMode;}
  Bool_t GetMCAnalysisMode() {return fAnalysisMC;}

  const  AliESDVertex *GetVertex(AliESDEvent *esd,
				 AnalysisMode mode,
				 Double_t gVx = 100.,
				 Double_t gVy = 100.,
				 Double_t gVz = 100.);
  void SetAcceptedVertexDiamond(Double_t gVx, Double_t gVy, Double_t gVz) {
    fVxMax = gVx; fVyMax = gVy; fVzMax = gVz;}
  Double_t GetVxMax() const {return fVxMax;}
  Double_t GetVyMax() const {return fVyMax;}
  Double_t GetVzMax() const {return fVzMax;}
  void SetMinNumOfContributors(Int_t nContributors) {
    fMinNumOfContributors = nContributors;}
  Int_t GetMinNumOfContributors() {return fMinNumOfContributors;}

  void SetPhaseSpace(Int_t nBinsX, Double_t gXmin, Double_t gXmax,
		     Int_t nBinsY, Double_t gYmin, Double_t gYmax) {
    fNBinsX = nBinsX; fMinX = gXmin; fMaxX = gXmax;
    fNBinsY = nBinsY; fMinY = gYmin; fMaxY = gYmax;
  }
  Int_t GetNBinsX() const {return fNBinsX;}
  Int_t GetNBinsY() const {return fNBinsY;}
  Double_t GetMinX() const {return fMinX;}
  Double_t GetMinY() const {return fMinY;}
  Double_t GetMaxX() const {return fMaxX;}
  Double_t GetMaxY() const {return fMaxY;}

  //Trigger
  Bool_t IsOnlineTriggerUsed() {return kUseOnlineTrigger;}
  void UseOnlineTrigger() {kUseOnlineTrigger = kTRUE;}
  Bool_t IsEventTriggered(const AliESDEvent *esd,
			  TriggerMode trigger = kMB2);
  void OfflineTriggerInit() {
    kUseOfflineTrigger = kTRUE;
    fPhysicsSelection = new AliPhysicsSelection();
    fPhysicsSelection->AddBackgroundIdentification(new AliBackgroundSelection());
    fPhysicsSelection->SetAnalyzeMC(fAnalysisMC);
  }
  Bool_t IsOfflineTriggerUsed() {return kUseOfflineTrigger;}
  AliPhysicsSelection *GetPhysicsSelectionObject() {return fPhysicsSelection;}

  Bool_t IsPrimary(AliESDEvent *esd,
		   const AliESDVertex *vertex, 
		   AliESDtrack *track);
  Bool_t IsAccepted(AliESDtrack *track);
  Bool_t IsInPhaseSpace(AliESDtrack *track);

  Float_t GetSigmaToVertex(AliESDtrack* esdTrack) const; 
  Double_t Rapidity(Double_t Px, Double_t Py, Double_t Pz) const;
  
  //Cut functions
  void    SetPointOnSPDLayers() {fPointOnSPDLayersFlag = kTRUE;}
  void    SetPointOnSDDLayers() {fPointOnSDDLayersFlag = kTRUE;}
  void    SetPointOnSSDLayers() {fPointOnSSDLayersFlag = kTRUE;}
  void    SetPointOnITSLayer1() {fPointOnITSLayer1Flag = kTRUE;}
  void    SetPointOnITSLayer2() {fPointOnITSLayer2Flag = kTRUE;}
  void    SetPointOnITSLayer3() {fPointOnITSLayer3Flag = kTRUE;}
  void    SetPointOnITSLayer4() {fPointOnITSLayer4Flag = kTRUE;}
  void    SetPointOnITSLayer5() {fPointOnITSLayer5Flag = kTRUE;}
  void    SetPointOnITSLayer6() {fPointOnITSLayer6Flag = kTRUE;}
  Bool_t  IsUsedPointOnSPDLayer() const {return fPointOnSPDLayersFlag;}
  Bool_t  IsUsedPointOnSDDLayer() const {return fPointOnSDDLayersFlag;}
  Bool_t  IsUsedPointOnSSDLayer() const {return fPointOnSSDLayersFlag;}
  Bool_t  IsUsedPointOnITSLayer1() const {return fPointOnITSLayer1Flag;}
  Bool_t  IsUsedPointOnITSLayer2() const {return fPointOnITSLayer2Flag;}
  Bool_t  IsUsedPointOnITSLayer3() const {return fPointOnITSLayer3Flag;}
  Bool_t  IsUsedPointOnITSLayer4() const {return fPointOnITSLayer4Flag;}
  Bool_t  IsUsedPointOnITSLayer5() const {return fPointOnITSLayer5Flag;}
  Bool_t  IsUsedPointOnITSLayer6() const {return fPointOnITSLayer6Flag;}
  void    SetMinITSClusters(Int_t minITSClusters) {
    fMinITSClusters = minITSClusters;
    fMinITSClustersFlag = kTRUE;
  }
  Int_t   GetMinITSClusters() const {return fMinITSClusters;}
  Bool_t  IsUsedMinITSClusters() const {return fMinITSClustersFlag;}

  void    SetMaxChi2PerITSCluster(Double_t maxChi2PerITSCluster) {
    fMaxChi2PerITSCluster = maxChi2PerITSCluster;
    fMaxChi2PerITSClusterFlag = kTRUE;
  }
  Bool_t  IsUsedMaxChi2PerITSCluster() const {return fMaxChi2PerITSClusterFlag;}
  Double_t   GetMaxChi2PerITSCluster() const {return fMaxChi2PerITSCluster;}

  void    SetMinTPCClusters(Int_t minTPCClusters) {
    fMinTPCClusters = minTPCClusters;
    fMinTPCClustersFlag = kTRUE;
  }
  Bool_t  IsUsedMinTPCClusters() const {return fMinTPCClustersFlag;}
  Int_t   GetMinTPCClusters() const {return fMinTPCClusters;}

  void    SetMaxChi2PerTPCCluster(Double_t maxChi2PerTPCCluster) {
    fMaxChi2PerTPCCluster = maxChi2PerTPCCluster;
    fMaxChi2PerTPCClusterFlag = kTRUE;
  }
  Bool_t  IsUsedMaxChi2PerTPCCluster() const {return fMaxChi2PerTPCClusterFlag;}
  Double_t   GetMaxChi2PerTPCCluster() const {return fMaxChi2PerTPCCluster;}

  void    SetMaxCov11(Double_t maxCov11) {
    fMaxCov11 = maxCov11; fMaxCov11Flag = kTRUE;}
  void    SetMaxCov22(Double_t maxCov22) {
    fMaxCov22 = maxCov22; fMaxCov22Flag = kTRUE;}
  void    SetMaxCov33(Double_t maxCov33) {
    fMaxCov33 = maxCov33; fMaxCov33Flag = kTRUE;}
  void    SetMaxCov44(Double_t maxCov44) {
    fMaxCov44 = maxCov44; fMaxCov44Flag = kTRUE;}
  void    SetMaxCov55(Double_t maxCov55) {
    fMaxCov55 = maxCov55; fMaxCov55Flag = kTRUE;}
  Bool_t  IsUsedMaxCov11() const {return fMaxCov11Flag;}
  Bool_t  IsUsedMaxCov22() const {return fMaxCov22Flag;}
  Bool_t  IsUsedMaxCov33() const {return fMaxCov33Flag;}
  Bool_t  IsUsedMaxCov44() const {return fMaxCov44Flag;}
  Bool_t  IsUsedMaxCov55() const {return fMaxCov55Flag;}
  Double_t   GetMaxCov11() const {return fMaxCov11;}
  Double_t   GetMaxCov22() const {return fMaxCov22;}
  Double_t   GetMaxCov33() const {return fMaxCov33;}
  Double_t   GetMaxCov44() const {return fMaxCov44;}
  Double_t   GetMaxCov55() const {return fMaxCov55;}

  void    SetMaxSigmaToVertex(Double_t maxSigmaToVertex) {
    fMaxSigmaToVertex = maxSigmaToVertex;
    fMaxSigmaToVertexFlag = kTRUE;
  }
  Bool_t  IsUsedMaxSigmaToVertex() const {return fMaxSigmaToVertexFlag;}
  Double_t   GetMaxSigmaToVertex() const {return fMaxSigmaToVertex;}

  void    SetMaxSigmaToVertexTPC(Double_t maxSigmaToVertex) {
    fMaxSigmaToVertexTPC = maxSigmaToVertex;
    fMaxSigmaToVertexTPCFlag = kTRUE;
  }
  Bool_t  IsUsedMaxSigmaToVertexTPC() const {return fMaxSigmaToVertexTPCFlag;}
  Double_t   GetMaxSigmaToVertexTPC() const {return fMaxSigmaToVertexTPC;}

  void    SetMaxDCAXY(Double_t maxDCAXY) {
    fMaxDCAXY = maxDCAXY;
    fMaxDCAXYFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCAXY() const {return fMaxDCAXYFlag;}
  Double_t   GetMaxDCAXY() const {return fMaxDCAXY;}

  void    SetMaxDCAXYTPC(Double_t maxDCAXY) {
    fMaxDCAXYTPC = maxDCAXY;
    fMaxDCAXYTPCFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCAXYTPC() const {return fMaxDCAXYTPCFlag;}
  Double_t   GetMaxDCAXYTPC() const {return fMaxDCAXYTPC;}

  void    SetMaxDCAZ(Double_t maxDCAZ) {
    fMaxDCAZ = maxDCAZ;
    fMaxDCAZFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCAZ() const {return fMaxDCAZFlag;}
  Double_t   GetMaxDCAZ() const {return fMaxDCAZ;}

  void    SetMaxDCAZTPC(Double_t maxDCAZ) {
    fMaxDCAZTPC = maxDCAZ;
    fMaxDCAZTPCFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCAZTPC() const {return fMaxDCAZTPCFlag;}
  Double_t   GetMaxDCAZTPC() const {return fMaxDCAZTPC;}

  void    SetMaxDCA3D(Double_t maxDCA3D) {
    fMaxDCA3D = maxDCA3D;
    fMaxDCA3DFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCA3D() const {return fMaxDCA3DFlag;}
  Double_t   GetMaxDCA3D() const {return fMaxDCA3D;}

  void SetPtDependentDCAxy(Int_t nSigma, Double_t p0, 
			   Double_t p1, Double_t p2);
  Bool_t  IsUsedPtDependentDCAxy() const {return fPtDependentDcaXYFlag;}

  void    SetMaxDCA3DTPC(Double_t maxDCA3D) {
    fMaxDCA3DTPC = maxDCA3D;
    fMaxDCA3DTPCFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCA3DTPC() const {return fMaxDCA3DTPCFlag;}
  Double_t   GetMaxDCA3DTPC() const {return fMaxDCA3DTPC;}

  void    SetMaxConstrainChi2(Double_t maxConstrainChi2) {
    fMaxConstrainChi2 = maxConstrainChi2;
    fMaxConstrainChi2Flag = kTRUE;
  }
  Bool_t  IsUsedMaxConstrainChi2() const {return fMaxConstrainChi2Flag;}
  Double_t   GetMaxConstrainChi2() const {return fMaxConstrainChi2;}

  void    SetMinTPCdEdxPoints(Int_t mindEdxpoints) {
    fMinTPCdEdxPoints = mindEdxpoints;
    fMinTPCdEdxPointsFlag = kTRUE;
  }
  Bool_t  IsUsedMinTPCdEdxPoints() const {return fMinTPCdEdxPointsFlag;}
  Int_t   GetMinTPCdEdxPoints() const {return fMinTPCdEdxPoints;}
  
  void    SetITSRefit() {fITSRefitFlag = kTRUE;}
  Bool_t  IsUsedITSRefit() const {return fITSRefitFlag;}
  void    SetTPCRefit() {fTPCRefitFlag = kTRUE;}
  Bool_t  IsUsedTPCRefit() const {return fTPCRefitFlag;}
  void    SetESDpid() {fESDpidFlag = kTRUE;}
  Bool_t  IsUsedESDpid() const {return fESDpidFlag;}
  void    SetTPCpid() {fTPCpidFlag = kTRUE;}
  Bool_t  IsUsedTPCpid() const {return fTPCpidFlag;}
  void    SetTOFpid() {fTOFpidFlag = kTRUE;}
  Bool_t  IsUsedTOFpid() const {return fTOFpidFlag;}

  TCanvas *GetListOfCuts();

  //PID related functions
  Bool_t IsProton(AliESDtrack *track);
  void SetNSigma(Int_t nsigma) {fNSigma = nsigma;}  
  Int_t GetNSigma() const {return fNSigma;}
  void SetRatio(Double_t ratio) {fNRatio = ratio;}
  Double_t GetRatio() {return fNRatio;}
  void SetPriorProbabilities(Double_t * const partFrac) {
    for(Int_t i = 0; i < AliPID::kSPECIESN; i++) fPartFrac[i] = partFrac[i];} 
  void SetPriorProbabilityFunctions(TF1 *const felectron, 
				    TF1 *const fmuon, 
				    TF1 *const fpion, 
				    TF1 *const fkaon, 
				    TF1 *const fproton) {
    fFunctionProbabilityFlag = kTRUE;
    fElectronFunction = felectron; fMuonFunction = fmuon; 
    fPionFunction = fpion; fKaonFunction = fkaon; fProtonFunction = fproton;
  }
  Bool_t IsPriorProbabilityFunctionUsed() const {return fFunctionProbabilityFlag;}
  Double_t GetParticleFraction(Int_t i, Double_t p);
  //Double_t Bethe(Double_t bg) const;

  void SetDebugMode() {fDebugMode = kTRUE;}
  Bool_t GetDebugMode() const {return fDebugMode;}

  void SetRunQA() {fRunQAAnalysis = kTRUE;}
  Bool_t IsQARun() {return fRunQAAnalysis;}
  TList *GetVertexQAList() {return fListVertexQA;}

 private:
  AliProtonAnalysisBase(const AliProtonAnalysisBase&); // Not implemented
  AliProtonAnalysisBase& operator=(const AliProtonAnalysisBase&); // Not implemented

  TString fProtonAnalysisLevel;//"ESD", "AOD" or "MC"
  Bool_t fAnalysisMC; //kTRUE if MC analysis while reading the ESDs
  TriggerMode fTriggerMode; //Trigger mode
  Bool_t kUseOnlineTrigger; //use the online trigger or not
  Bool_t kUseOfflineTrigger; //use the offline trigger or not
  AliPhysicsSelection *fPhysicsSelection; //Trigger selection: offline
  AnalysisMode fProtonAnalysisMode; //Analysis mode: TPC-Hybrid-Global
  PIDMode fProtonPIDMode; //PID mode: Bayesian-dE/dx ratio-Nsigma areas
  Bool_t fAnalysisEtaMode; //run the analysis in eta or y

  Bool_t fRunQAAnalysis; //boolnean to indicate to run the QA or not
  Double_t fVxMax, fVyMax, fVzMax; //vertex diamond constrain 
  Int_t fMinNumOfContributors;//min number of contributors

  Int_t fNBinsX; //number of bins in y or eta
  Double_t fMinX, fMaxX; //min & max value of y or eta
  Int_t fNBinsY;  //number of bins in pT
  Double_t fMinY, fMaxY; //min & max value of pT
  
  //cuts
  Int_t fMinTPCClusters, fMinITSClusters; //min TPC & ITS clusters
  Double_t fMaxChi2PerTPCCluster, fMaxChi2PerITSCluster; //max chi2 per TPC & ITS cluster
  Double_t fMaxCov11, fMaxCov22, fMaxCov33, fMaxCov44, fMaxCov55; //max values of cov. matrix
  Double_t fMaxSigmaToVertex; //max sigma to vertex cut
  Double_t fMaxSigmaToVertexTPC; //max sigma to vertex cut
  Double_t fMaxDCAXY, fMaxDCAXYTPC; //max DCA xy
  Double_t fMaxDCAZ, fMaxDCAZTPC; //max DCA z
  Double_t fMaxDCA3D, fMaxDCA3DTPC; //max DCA 3D
  Double_t fMaxConstrainChi2; //max constrain chi2 - vertex
  Int_t  fMinTPCdEdxPoints;//min number of TPC points used for the dE/dx
  Bool_t fMinTPCClustersFlag, fMinITSClustersFlag; //shows if this cut is used or not
  Bool_t fMaxChi2PerTPCClusterFlag, fMaxChi2PerITSClusterFlag; //shows if this cut is used or not
  Bool_t fMaxCov11Flag, fMaxCov22Flag, fMaxCov33Flag, fMaxCov44Flag, fMaxCov55Flag; //shows if this cut is used or not
  Bool_t fMaxSigmaToVertexFlag; //shows if this cut is used or not
  Bool_t fMaxSigmaToVertexTPCFlag; //shows if this cut is used or not
  Bool_t fMaxDCAXYFlag, fMaxDCAXYTPCFlag; //shows if this cut is used or not
  Bool_t fMaxDCAZFlag, fMaxDCAZTPCFlag; //shows if this cut is used or not
  Bool_t fMaxDCA3DFlag, fMaxDCA3DTPCFlag; //shows if this cut is used or not
  Bool_t fMaxConstrainChi2Flag; //shows if this cut is used or not
  Bool_t fITSRefitFlag, fTPCRefitFlag; //shows if this cut is used or not
  Bool_t fESDpidFlag, fTPCpidFlag, fTOFpidFlag; //shows if this cut is used or not
  Bool_t fPointOnSPDLayersFlag;//shows if this cut is used or not
  Bool_t fPointOnSDDLayersFlag;//shows if this cut is used or not
  Bool_t fPointOnSSDLayersFlag;//shows if this cut is used or not
  Bool_t fPointOnITSLayer1Flag, fPointOnITSLayer2Flag; //shows if this cut is used or not
  Bool_t fPointOnITSLayer3Flag, fPointOnITSLayer4Flag; //shows if this cut is used or not
  Bool_t fPointOnITSLayer5Flag, fPointOnITSLayer6Flag; //shows if this cut is used or not
  Bool_t fMinTPCdEdxPointsFlag; //shows if this cut is used or not
  TF1  *fPtDependentDcaXY; //pt dependence dca cut (xy)
  Bool_t fPtDependentDcaXYFlag; //shows if this cut is used or not
  Int_t fNSigmaDCAXY; //n-sigma dca xy cut (pt dependent)

  //pid
  Bool_t fFunctionProbabilityFlag; //flag: kTRUE if functions used
  Int_t fNSigma; //N-sigma cut in the dE/dx band
  Double_t fNRatio; //min value of the ratio of the measured dE/dx vs the expected
  Double_t fPartFrac[10]; //prior probabilities
  TF1  *fElectronFunction; //momentum dependence of the prior probs
  TF1  *fMuonFunction; //momentum dependence of the prior probs
  TF1  *fPionFunction; //momentum dependence of the prior probs
  TF1  *fKaonFunction; //momentum dependence of the prior probs
  TF1  *fProtonFunction; //momentum dependence of the prior probs

  //Debug
  Bool_t fDebugMode; //Enable the debug mode

  //QA list
  TList *fListVertexQA; //vertex QA list

  ClassDef(AliProtonAnalysisBase,1);
};

#endif
