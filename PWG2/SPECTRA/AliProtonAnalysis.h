#ifndef ALIPROTONANALYSIS_H
#define ALIPROTONANALYSIS_H

/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                       Class AliProtonAnalysis
//   This is the class for the baryon (proton) analysis
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include "TObject.h"
#include "TH1I.h"
#include "TList.h"

#include "AliPID.h"

class TF1;
class TH2D;
class TH1F;

#include "AliCFContainer.h"
class AliCFDataGrid;
class AliAODEvent;
class AliAODtrack;
class AliESDEvent;
class AliESDtrack;
class AliExternalTrackParam;
class AliStack;
class AliESDVertex;

class AliProtonAnalysis : public TObject {
 public:
  AliProtonAnalysis();
  AliProtonAnalysis(Int_t nbinsY, Float_t fLowY, Float_t fHighY,
		    Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
  virtual ~AliProtonAnalysis();

  void UseTPCOnly() {fUseTPCOnly = kTRUE;}
  void UseHybridTPC() {fUseHybridTPC = kTRUE;}
  
  void InitAnalysisHistograms(Int_t nbinsY, Float_t fLowY, Float_t fHighY,
			      Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
  Bool_t ReadFromFile(const char* filename);
  void Analyze(AliESDEvent *fESD, 
	       const AliESDVertex *vertex);
  void Analyze(AliAODEvent *fAOD);
  void Analyze(AliStack *stack);
  
  AliCFContainer *GetProtonContainer() {return fProtonContainer;}
  AliCFContainer *GetAntiProtonContainer() {return fAntiProtonContainer;}

  TH2D *GetProtonYPtHistogram() {return fHistYPtProtons;}
  TH2D *GetAntiProtonYPtHistogram() {return fHistYPtAntiProtons;}
  TH1D *GetProtonYHistogram();
  TH1D *GetAntiProtonYHistogram();
  TH1D *GetProtonPtHistogram();
  TH1D *GetAntiProtonPtHistogram();
  TH1D *GetProtonCorrectedYHistogram();
  TH1D *GetAntiProtonCorrectedYHistogram();
  TH1D *GetProtonCorrectedPtHistogram();
  TH1D *GetAntiProtonCorrectedPtHistogram();
  
  TH1D *GetYRatioHistogram();
  TH1D *GetPtRatioHistogram();
  TH1D *GetYAsymmetryHistogram();
  TH1D *GetPtAsymmetryHistogram();

  TH1I *GetEventHistogram() {return fHistEvents;}

  Int_t   GetNumberOfAnalyzedEvents() {return (Int_t)fHistEvents->GetEntries();} 
  Bool_t  PrintMean(TH1 *hist, Double_t edge);
  Bool_t  PrintYields(TH1 *hist, Double_t edge); 

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
  void    SetMaxConstrainChi2(Double_t maxConstrainChi2) {
    fMaxConstrainChi2 = maxConstrainChi2;
    fMaxConstrainChi2Flag = kTRUE;
  }
  void    SetITSRefit() {fITSRefitFlag = kTRUE;}
  void    SetTPCRefit() {fTPCRefitFlag = kTRUE;}
  void    SetESDpid() {fESDpidFlag = kTRUE;}
  void    SetTPCpid() {fTPCpidFlag = kTRUE;}

  //Prior probabilities
  void SetPriorProbabilities(Double_t *partFrac) {
    for(Int_t i = 0; i < AliPID::kSPECIESN; i++) fPartFrac[i] = partFrac[i];} 
  void SetPriorProbabilityFunctions(TF1 *felectron, TF1 *fmuon, TF1 *fpion, TF1 *fkaon, TF1 *fproton) {
    fFunctionProbabilityFlag = kTRUE;
    fElectronFunction = felectron; 
    fMuonFunction = fmuon; 
    fPionFunction = fpion;
    fKaonFunction = fkaon;
    fProtonFunction = fproton;
  } 
  Double_t GetParticleFraction(Int_t i, Double_t p);

  //interface to the correction framework
  void Correct(Int_t step);
  Bool_t ReadCorrectionContainer(const char* filename);
  TList *GetCorrectionListProtons2D() {return fCorrectionListProtons2D;} 
  TList *GetEfficiencyListProtons1D() {return fEfficiencyListProtons1D;} 
  TList *GetCorrectionListProtons1D() {return fCorrectionListProtons1D;} 
  TList *GetCorrectionListAntiProtons2D() {return fCorrectionListAntiProtons2D;} 
  TList *GetEfficiencyListAntiProtons1D() {return fEfficiencyListAntiProtons1D;} 
  TList *GetCorrectionListAntiProtons1D() {return fCorrectionListAntiProtons1D;} 
  
  //iStep=0->MC - iStep=1->Acceptance - iStep=2->Reconstruction - iStep=3->PID
  TH1D  *GetUncorrectedProtonYHistogram(Int_t iStep) {return fProtonContainer->ShowProjection(0, iStep);}
  TH1D  *GetUncorrectedProtonPtHistogram(Int_t iStep) {return fProtonContainer->ShowProjection(1, iStep);}
  TH1D  *GetUncorrectedAntiProtonYHistogram(Int_t iStep) {return fAntiProtonContainer->ShowProjection(0, iStep);}
  TH1D  *GetUncorrectedAntiProtonPtHistogram(Int_t iStep) {return fAntiProtonContainer->ShowProjection(1, iStep);}

 private:
  AliProtonAnalysis(const AliProtonAnalysis&); // Not implemented
  AliProtonAnalysis& operator=(const AliProtonAnalysis&); // Not implemented

  Bool_t   IsAccepted(AliESDEvent *esd,
		      const AliESDVertex *vertex, 
		      AliESDtrack *track);
  Float_t  GetSigmaToVertex(AliESDtrack* esdTrack); 
  Double_t Rapidity(Double_t Px, Double_t Py, Double_t Pz);
  
  Int_t fNBinsY; //number of bins in y
  Float_t fMinY, fMaxY; //min & max value of y
  Int_t fNBinsPt;  //number of bins in pT
  Float_t fMinPt, fMaxPt; //min & max value of pT
  
  //cuts
  Int_t fMinTPCClusters, fMinITSClusters; //min TPC & ITS clusters
  Double_t fMaxChi2PerTPCCluster, fMaxChi2PerITSCluster; //max chi2 per TPC & ITS cluster
  Double_t fMaxCov11, fMaxCov22, fMaxCov33, fMaxCov44, fMaxCov55; //max values of cov. matrix
  Double_t fMaxSigmaToVertex; //max sigma to vertex cut
  Double_t fMaxSigmaToVertexTPC; //max sigma to vertex cut
  Double_t fMaxDCAXY, fMaxDCAXYTPC; //max DCA xy
  Double_t fMaxDCAZ, fMaxDCAZTPC; //max DCA z
  Double_t fMaxConstrainChi2; //max constrain chi2 - vertex
  Bool_t fMinTPCClustersFlag, fMinITSClustersFlag; //shows if this cut is used or not
  Bool_t fMaxChi2PerTPCClusterFlag, fMaxChi2PerITSClusterFlag; //shows if this cut is used or not
  Bool_t fMaxCov11Flag, fMaxCov22Flag, fMaxCov33Flag, fMaxCov44Flag, fMaxCov55Flag; //shows if this cut is used or not
  Bool_t fMaxSigmaToVertexFlag; //shows if this cut is used or not
  Bool_t fMaxSigmaToVertexTPCFlag; //shows if this cut is used or not
  Bool_t fMaxDCAXYFlag, fMaxDCAXYTPCFlag; //shows if this cut is used or not
  Bool_t fMaxDCAZFlag, fMaxDCAZTPCFlag; //shows if this cut is used or not
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

  //Detectors
  Bool_t fUseTPCOnly; //kTRUE if TPC only information is used
  Bool_t fUseHybridTPC; //kTRUE if TPC info is used for momentum - PID and ITS for vertex & points

  //Analysis containers
  AliCFContainer *fProtonContainer; //container for protons
  AliCFContainer *fAntiProtonContainer; //container for antiprotons
  TH1I *fHistEvents; //event counter
  TH2D *fHistYPtProtons; //Y-Pt of Protons
  TH2D *fHistYPtAntiProtons; // Y-Pt of Antiprotons

  //Corrections
  TList *fEffGridListProtons; //list for the efficiency grid - protons 
  TList *fCorrectionListProtons2D; //list for the 2d corrections 
  TList *fEfficiencyListProtons1D; //list for the 1d efficiencies
  TList *fCorrectionListProtons1D; //list for the 1d corrections 
  TList *fEffGridListAntiProtons; //list for the efficiency grid - antiprotons 
  TList *fCorrectionListAntiProtons2D; //list for the 2d corrections 
  TList *fEfficiencyListAntiProtons1D; //list for the 1d efficiencies
  TList *fCorrectionListAntiProtons1D; //list for the 1d corrections 
  AliCFDataGrid *fCorrectProtons; //corrected data grid for protons
  AliCFDataGrid *fCorrectAntiProtons; //corrected data grid for antiprotons

  ClassDef(AliProtonAnalysis,0);
};

#endif
