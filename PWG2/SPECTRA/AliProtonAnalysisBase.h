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
  enum PIDMode { kBayesian = 0, kRatio, kSigma1, kSigma2 };

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
  Bool_t  IsUsedPointOnITSLayer1() {return fPointOnITSLayer1Flag;}
  Bool_t  IsUsedPointOnITSLayer2() {return fPointOnITSLayer2Flag;}
  Bool_t  IsUsedPointOnITSLayer3() {return fPointOnITSLayer3Flag;}
  Bool_t  IsUsedPointOnITSLayer4() {return fPointOnITSLayer4Flag;}
  Bool_t  IsUsedPointOnITSLayer5() {return fPointOnITSLayer5Flag;}
  Bool_t  IsUsedPointOnITSLayer6() {return fPointOnITSLayer6Flag;}
  void    SetMinITSClusters(Int_t minITSClusters) {
    fMinITSClusters = minITSClusters;
    fMinITSClustersFlag = kTRUE;
  }
  Int_t   GetMinITSClusters() {return fMinITSClusters;}
  Bool_t  IsUsedMinITSClusters() {return fMinITSClustersFlag;}

  void    SetMaxChi2PerITSCluster(Double_t maxChi2PerITSCluster) {
    fMaxChi2PerITSCluster = maxChi2PerITSCluster;
    fMaxChi2PerITSClusterFlag = kTRUE;
  }
  Bool_t  IsUsedMaxChi2PerITSCluster() {return fMaxChi2PerITSClusterFlag;}
  Double_t   GetMaxChi2PerITSCluster() {return fMaxChi2PerITSCluster;}

  void    SetMinTPCClusters(Int_t minTPCClusters) {
    fMinTPCClusters = minTPCClusters;
    fMinTPCClustersFlag = kTRUE;
  }
  Bool_t  IsUsedMinTPCClusters() {return fMinTPCClustersFlag;}
  Int_t   GetMinTPCClusters() {return fMinTPCClusters;}

  void    SetMaxChi2PerTPCCluster(Double_t maxChi2PerTPCCluster) {
    fMaxChi2PerTPCCluster = maxChi2PerTPCCluster;
    fMaxChi2PerTPCClusterFlag = kTRUE;
  }
  Bool_t  IsUsedMaxChi2PerTPCCluster() {return fMaxChi2PerTPCClusterFlag;}
  Double_t   GetMaxChi2PerTPCCluster() {return fMaxChi2PerTPCCluster;}

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
  Bool_t  IsUsedMaxCov11() {return fMaxCov11Flag;}
  Bool_t  IsUsedMaxCov22() {return fMaxCov22Flag;}
  Bool_t  IsUsedMaxCov33() {return fMaxCov33Flag;}
  Bool_t  IsUsedMaxCov44() {return fMaxCov44Flag;}
  Bool_t  IsUsedMaxCov55() {return fMaxCov55Flag;}
  Double_t   GetMaxCov11() {return fMaxCov11;}
  Double_t   GetMaxCov22() {return fMaxCov22;}
  Double_t   GetMaxCov33() {return fMaxCov33;}
  Double_t   GetMaxCov44() {return fMaxCov44;}
  Double_t   GetMaxCov55() {return fMaxCov55;}

  void    SetMaxSigmaToVertex(Double_t maxSigmaToVertex) {
    fMaxSigmaToVertex = maxSigmaToVertex;
    fMaxSigmaToVertexFlag = kTRUE;
  }
  Bool_t  IsUsedMaxSigmaToVertex() {return fMaxSigmaToVertexFlag;}
  Double_t   GetMaxSigmaToVertex() {return fMaxSigmaToVertex;}

  void    SetMaxSigmaToVertexTPC(Double_t maxSigmaToVertex) {
    fMaxSigmaToVertexTPC = maxSigmaToVertex;
    fMaxSigmaToVertexTPCFlag = kTRUE;
  }
  Bool_t  IsUsedMaxSigmaToVertexTPC() {return fMaxSigmaToVertexTPCFlag;}
  Double_t   GetMaxSigmaToVertexTPC() {return fMaxSigmaToVertexTPC;}

  void    SetMaxDCAXY(Double_t maxDCAXY) {
    fMaxDCAXY = maxDCAXY;
    fMaxDCAXYFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCAXY() {return fMaxDCAXYFlag;}
  Double_t   GetMaxDCAXY() {return fMaxDCAXY;}

  void    SetMaxDCAXYTPC(Double_t maxDCAXY) {
    fMaxDCAXYTPC = maxDCAXY;
    fMaxDCAXYTPCFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCAXYTPC() {return fMaxDCAXYTPCFlag;}
  Double_t   GetMaxDCAXYTPC() {return fMaxDCAXYTPC;}

  void    SetMaxDCAZ(Double_t maxDCAZ) {
    fMaxDCAZ = maxDCAZ;
    fMaxDCAZFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCAZ() {return fMaxDCAZFlag;}
  Double_t   GetMaxDCAZ() {return fMaxDCAZ;}

  void    SetMaxDCAZTPC(Double_t maxDCAZ) {
    fMaxDCAZTPC = maxDCAZ;
    fMaxDCAZTPCFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCAZTPC() {return fMaxDCAZTPCFlag;}
  Double_t   GetMaxDCAZTPC() {return fMaxDCAZTPC;}

  void    SetMaxDCA3D(Double_t maxDCA3D) {
    fMaxDCA3D = maxDCA3D;
    fMaxDCA3DFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCA3D() {return fMaxDCA3DFlag;}
  Double_t   GetMaxDCA3D() {return fMaxDCA3D;}

  void    SetMaxDCA3DTPC(Double_t maxDCA3D) {
    fMaxDCA3DTPC = maxDCA3D;
    fMaxDCA3DTPCFlag = kTRUE;
  }
  Bool_t  IsUsedMaxDCA3DTPC() {return fMaxDCA3DTPCFlag;}
  Double_t   GetMaxDCA3DTPC() {return fMaxDCA3DTPC;}

  void    SetMaxConstrainChi2(Double_t maxConstrainChi2) {
    fMaxConstrainChi2 = maxConstrainChi2;
    fMaxConstrainChi2Flag = kTRUE;
  }
  Bool_t  IsUsedMaxConstrainChi2() {return fMaxConstrainChi2Flag;}
  Double_t   GetMaxConstrainChi2() {return fMaxConstrainChi2;}

  void    SetMinTPCdEdxPoints(Int_t mindEdxpoints) {
    fMinTPCdEdxPoints = mindEdxpoints;
    fMinTPCdEdxPointsFlag = kTRUE;
  }
  Bool_t  IsUsedMinTPCdEdxPoints() {return fMinTPCdEdxPointsFlag;}
  Int_t   GetMinTPCdEdxPoints() {return fMinTPCdEdxPoints;}
  
  void    SetITSRefit() {fITSRefitFlag = kTRUE;}
  Bool_t  IsUsedITSRefit() {return fITSRefitFlag;}
  void    SetTPCRefit() {fTPCRefitFlag = kTRUE;}
  Bool_t  IsUsedTPCRefit() {return fTPCRefitFlag;}
  void    SetESDpid() {fESDpidFlag = kTRUE;}
  Bool_t  IsUsedESDpid() {return fESDpidFlag;}
  void    SetTPCpid() {fTPCpidFlag = kTRUE;}
  Bool_t  IsUsedTPCpid() {return fTPCpidFlag;}

  TCanvas *GetListOfCuts();

  //PID related functions
  Bool_t IsProton(AliESDtrack *track);
  void SetNSigma(Int_t nsigma) {fNSigma = nsigma;}  
  Int_t GetNSigma() {return fNSigma;}
  void SetdEdxBandInfo(const char* filename);
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
  Double_t Bethe(Double_t bg);

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
  Bool_t fESDpidFlag, fTPCpidFlag; //shows if this cut is used or not
  Bool_t fPointOnITSLayer1Flag, fPointOnITSLayer2Flag; //shows if this cut is used or not
  Bool_t fPointOnITSLayer3Flag, fPointOnITSLayer4Flag; //shows if this cut is used or not
  Bool_t fPointOnITSLayer5Flag, fPointOnITSLayer6Flag; //shows if this cut is used or not
  Bool_t fMinTPCdEdxPointsFlag; //shows if this cut is used or not

  //pid
  Bool_t fFunctionProbabilityFlag; //flag: kTRUE if functions used
  Int_t fNSigma; //N-sigma cut in the dE/dx band
  Double_t fdEdxMean[24]; //mean values of the dE/dx distributions for the proton band - P slices
  Double_t fdEdxSigma[24]; //sigma values of the dE/dx distributions for the proton band - P slices
  Double_t fPartFrac[10]; //prior probabilities
  TF1  *fElectronFunction; //momentum dependence of the prior probs
  TF1  *fMuonFunction; //momentum dependence of the prior probs
  TF1  *fPionFunction; //momentum dependence of the prior probs
  TF1  *fKaonFunction; //momentum dependence of the prior probs
  TF1  *fProtonFunction; //momentum dependence of the prior probs

  //Debug
  Bool_t fDebugMode; //Enable the debug mode

  ClassDef(AliProtonAnalysisBase,1);
};

#endif
