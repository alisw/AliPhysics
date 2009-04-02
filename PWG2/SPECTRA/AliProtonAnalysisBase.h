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

#include "AliPID.h"
class AliESDEvent;
class AliESDtrack;
class AliESDVertex;

class AliProtonAnalysisBase : public TObject {
 public:
  enum TriggerMode { kMB1 = 0, kMB2, kSPDFASTOR };
  enum AnalysisMode { kInvalid = -1, kTPC = 0, kHybrid, kGlobal };
  enum PIDMode { kBayesian = 0, kRatio, kSigma };

  AliProtonAnalysisBase();
  virtual ~AliProtonAnalysisBase();

  void SetAnalysisLevel(const char* type) {fProtonAnalysisLevel = type;}
  void SetAnalysisMode(AnalysisMode analysismode) {fProtonAnalysisMode = analysismode;}
  void SetEtaMode() {fAnalysisEtaMode = kTRUE;}
  void SetTriggerMode(TriggerMode triggermode) {fTriggerMode = triggermode;}
  void SetPIDMode(PIDMode pidmode) {fProtonPIDMode = pidmode;}

  const char *GetAnalysisLevel() {return fProtonAnalysisLevel.Data();}
  AnalysisMode GetAnalysisMode() {return fProtonAnalysisMode;}
  Bool_t GetEtaMode() {return fAnalysisEtaMode;}
  TriggerMode GetTriggerMode() {return fTriggerMode;}
  PIDMode GetPIDMode() {return fProtonPIDMode;}

  const  AliESDVertex *GetVertex(AliESDEvent *esd,
				 AnalysisMode mode,
				 Double_t gVx = 100.,
				 Double_t gVy = 100.,
				 Double_t gVz = 100.);
  void SetAcceptedVertexDiamond(Double_t gVx, Double_t gVy, Double_t gVz) {
    fVxMax = gVx; fVyMax = gVy; fVzMax = gVz;}
  Double_t GetVxMax() {return fVxMax;}
  Double_t GetVyMax() {return fVyMax;}
  Double_t GetVzMax() {return fVzMax;}

  void SetPhaseSpace(Int_t nBinsX, Double_t gXmin, Double_t gXmax,
		     Int_t nBinsY, Double_t gYmin, Double_t gYmax) {
    fNBinsX = nBinsX; fMinX = gXmin; fMaxX = gXmax;
    fNBinsY = nBinsY; fMinY = gYmin; fMaxY = gYmax;
  }
  Int_t GetNBinsX() {return fNBinsX;}
  Int_t GetNBinsY() {return fNBinsY;}
  Double_t GetMinX() {return fMinX;}
  Double_t GetMinY() {return fMinY;}
  Double_t GetMaxX() {return fMaxX;}
  Double_t GetMaxY() {return fMaxY;}

  static Bool_t IsEventTriggered(const AliESDEvent *esd,
                                 TriggerMode trigger = kMB2);
  Bool_t IsAccepted(AliESDEvent *esd,
		    const AliESDVertex *vertex, 
		    AliESDtrack *track);
  Bool_t IsInPhaseSpace(AliESDtrack *track);

  Float_t GetSigmaToVertex(AliESDtrack* esdTrack) const; 
  Double_t Rapidity(Double_t Px, Double_t Py, Double_t Pz) const;
  
  //Cut functions
  void    SetPointOnITSLayer1() {fPointOnITSLayer1Flag = kTRUE;}
  void    SetPointOnITSLayer2() {fPointOnITSLayer2Flag = kTRUE;}
  void    SetPointOnITSLayer3() {fPointOnITSLayer3Flag = kTRUE;}
  void    SetPointOnITSLayer4() {fPointOnITSLayer4Flag = kTRUE;}
  void    SetPointOnITSLayer5() {fPointOnITSLayer5Flag = kTRUE;}
  void    SetPointOnITSLayer6() {fPointOnITSLayer6Flag = kTRUE;}
  void    SetMinITSClusters(Int_t minITSClusters) {
    fMinITSClusters = minITSClusters;
    fMinITSClustersFlag = kTRUE;
  }
  void    SetMaxChi2PerITSCluster(Double_t maxChi2PerITSCluster) {
    fMaxChi2PerITSCluster = maxChi2PerITSCluster;
    fMaxChi2PerITSClusterFlag = kTRUE;
  }
  void    SetMinTPCClusters(Int_t minTPCClusters) {
    fMinTPCClusters = minTPCClusters;
    fMinTPCClustersFlag = kTRUE;
  }
  void    SetMaxChi2PerTPCCluster(Double_t maxChi2PerTPCCluster) {
    fMaxChi2PerTPCCluster = maxChi2PerTPCCluster;
    fMaxChi2PerTPCClusterFlag = kTRUE;
  }
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
  void    SetMaxSigmaToVertex(Double_t maxSigmaToVertex) {
    fMaxSigmaToVertex = maxSigmaToVertex;
    fMaxSigmaToVertexFlag = kTRUE;
  }
  void    SetMaxSigmaToVertexTPC(Double_t maxSigmaToVertex) {
    fMaxSigmaToVertexTPC = maxSigmaToVertex;
    fMaxSigmaToVertexTPCFlag = kTRUE;
  }
  void    SetMaxDCAXY(Double_t maxDCAXY) {
    fMaxDCAXY = maxDCAXY;
    fMaxDCAXYFlag = kTRUE;
  }
  void    SetMaxDCAXYTPC(Double_t maxDCAXY) {
    fMaxDCAXYTPC = maxDCAXY;
    fMaxDCAXYTPCFlag = kTRUE;
  }
  void    SetMaxDCAZ(Double_t maxDCAZ) {
    fMaxDCAZ = maxDCAZ;
    fMaxDCAZFlag = kTRUE;
  }
  void    SetMaxDCAZTPC(Double_t maxDCAZ) {
    fMaxDCAZTPC = maxDCAZ;
    fMaxDCAZTPCFlag = kTRUE;
  }
  void    SetMaxDCA3D(Double_t maxDCA3D) {
    fMaxDCA3D = maxDCA3D;
    fMaxDCA3DFlag = kTRUE;
  }
  void    SetMaxDCA3DTPC(Double_t maxDCA3D) {
    fMaxDCA3DTPC = maxDCA3D;
    fMaxDCA3DTPCFlag = kTRUE;
  }
  void    SetMaxConstrainChi2(Double_t maxConstrainChi2) {
    fMaxConstrainChi2 = maxConstrainChi2;
    fMaxConstrainChi2Flag = kTRUE;
  }
  void    SetITSRefit() {fITSRefitFlag = kTRUE;}
  void    SetTPCRefit() {fTPCRefitFlag = kTRUE;}
  void    SetESDpid() {fESDpidFlag = kTRUE;}
  void    SetTPCpid() {fTPCpidFlag = kTRUE;}

  TCanvas *GetListOfCuts();

  //PID related functions
  Bool_t IsProton(AliESDtrack *track);
  void SetPriorProbabilities(Double_t * const partFrac) {
    for(Int_t i = 0; i < AliPID::kSPECIESN; i++) fPartFrac[i] = partFrac[i];} 
  void SetPriorProbabilityFunctions(TF1 *const felectron, 
				    TF1 *const fmuon, 
				    TF1 *const fpion, 
				    TF1 *const fkaon, 
				    TF1 *const fproton) {
    fFunctionProbabilityFlag = kTRUE;
    fElectronFunction = felectron; 
    fMuonFunction = fmuon; 
    fPionFunction = fpion;
    fKaonFunction = fkaon;
    fProtonFunction = fproton;
  }
  Bool_t IsPriorProbabilityFunctionUsed() {return fFunctionProbabilityFlag;}
  Double_t GetParticleFraction(Int_t i, Double_t p);

  void SetDebugMode() {fDebugMode = kTRUE;}
  Bool_t GetDebugMode() {return fDebugMode;}

 private:
  AliProtonAnalysisBase(const AliProtonAnalysisBase&); // Not implemented
  AliProtonAnalysisBase& operator=(const AliProtonAnalysisBase&); // Not implemented

  TString fProtonAnalysisLevel;//"ESD", "AOD" or "MC"
  TriggerMode fTriggerMode; //Trigger mode
  AnalysisMode fProtonAnalysisMode; //Analysis mode: TPC-Hybrid-Global
  PIDMode fProtonPIDMode; //PID mode: Bayesian-dE/dx ratio-Nsigma areas
  Bool_t fAnalysisEtaMode; //run the analysis in eta or y

  Double_t fVxMax, fVyMax, fVzMax; //vertex diamond constrain 

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
  Bool_t fESDpidFlag, fTPCpidFlag; //shows if this cut is used or not
  Bool_t fPointOnITSLayer1Flag, fPointOnITSLayer2Flag; //shows if this cut is used or not
  Bool_t fPointOnITSLayer3Flag, fPointOnITSLayer4Flag; //shows if this cut is used or not
  Bool_t fPointOnITSLayer5Flag, fPointOnITSLayer6Flag; //shows if this cut is used or not
  
  //pid
  Bool_t fFunctionProbabilityFlag; //flag: kTRUE if functions used
  Double_t fPartFrac[10]; //prior probabilities
  TF1  *fElectronFunction; //momentum dependence of the prior probs
  TF1  *fMuonFunction; //momentum dependence of the prior probs
  TF1  *fPionFunction; //momentum dependence of the prior probs
  TF1  *fKaonFunction; //momentum dependence of the prior probs
  TF1  *fProtonFunction; //momentum dependence of the prior probs

  //Debug
  Bool_t fDebugMode; //Enable the debug mode

  ClassDef(AliProtonAnalysisBase,0);
};

#endif
