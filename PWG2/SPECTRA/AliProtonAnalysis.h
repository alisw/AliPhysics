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

#include <TObject.h>
#include "AliPID.h"

class TF1;
class TH2F;
class TH1D;
class AliAODEvent;
class AliAODtrack;
class AliESDEvent;
class AliESDtrack;
class AliExternalTrackParam;

class AliProtonAnalysis : public TObject {
 public:
  AliProtonAnalysis();
  AliProtonAnalysis(Int_t nbinsY, Float_t fLowY, Float_t fHighY,
		    Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
  virtual ~AliProtonAnalysis();

  void UseTPCOnly() {fUseTPCOnly = kTRUE;}
  
  void InitHistograms(Int_t nbinsY, Float_t fLowY, Float_t fHighY,
		      Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
  Bool_t ReadFromFile(const char* filename);
  void Analyze(AliESDEvent* fESD);
  void Analyze(AliAODEvent* fAOD);
  
  TH2F *GetProtonYPtHistogram() {return fHistYPtProtons;}
  TH2F *GetAntiProtonYPtHistogram() {return fHistYPtAntiProtons;}
  TH1D *GetProtonYHistogram();
  TH1D *GetAntiProtonYHistogram();
  TH1D *GetProtonPtHistogram();
  TH1D *GetAntiProtonPtHistogram();
  TH1D *GetYRatioHistogram();
  TH1D *GetPtRatioHistogram();
  TH1D *GetYAsymmetryHistogram();
  TH1D *GetPtAsymmetryHistogram();
  
  //Cut functions
  void    SetMinTPCClusters(Int_t minTPCClusters) {
    fMinTPCClusters = minTPCClusters;
    fMinTPCClustersFlag = kTRUE;
  }
  void    SetMinITSClusters(Int_t minITSClusters) {
    fMinITSClusters = minITSClusters;
    fMinITSClustersFlag = kTRUE;
  }
  void    SetMaxChi2PerTPCCluster(Double_t maxChi2PerTPCCluster) {
    fMaxChi2PerTPCCluster = maxChi2PerTPCCluster;
    fMaxChi2PerTPCClusterFlag = kTRUE;
  }
  void    SetMaxChi2PerITSCluster(Double_t maxChi2PerITSCluster) {
    fMaxChi2PerITSCluster = maxChi2PerITSCluster;
    fMaxChi2PerITSClusterFlag = kTRUE;
  }
  void    SetMaxCov11(Double_t maxCov11) {fMaxCov11 = maxCov11; fMaxCov11Flag = kTRUE;}
  void    SetMaxCov22(Double_t maxCov22) {fMaxCov22 = maxCov22; fMaxCov22Flag = kTRUE;}
  void    SetMaxCov33(Double_t maxCov33) {fMaxCov33 = maxCov33; fMaxCov33Flag = kTRUE;}
  void    SetMaxCov44(Double_t maxCov44) {fMaxCov44 = maxCov44; fMaxCov44Flag = kTRUE;}
  void    SetMaxCov55(Double_t maxCov55) {fMaxCov55 = maxCov55; fMaxCov55Flag = kTRUE;}
  void    SetMaxSigmaToVertex(Double_t maxSigmaToVertex) {
    fMaxSigmaToVertex = maxSigmaToVertex;
    fMaxSigmaToVertexFlag = kTRUE;
  }
  void    SetITSRefit() {fITSRefitFlag = kTRUE;}
  void    SetTPCRefit() {fTPCRefitFlag = kTRUE;}
  
  //Prior probabilities
  void    SetPriorProbabilities(Double_t *partFrac) {for(Int_t i = 0; i < AliPID::kSPECIESN; i++) fPartFrac[i] = partFrac[i];} 
  void    SetPriorProbabilityFunctions(TF1 *felectron, TF1 *fmuon, TF1 *fpion, TF1 *fkaon, TF1 *fproton) {
    fFunctionProbabilityFlag = kTRUE;
    fElectronFunction = felectron; 
    fMuonFunction = fmuon; 
    fPionFunction = fpion;
    fKaonFunction = fkaon;
    fProtonFunction = fproton;
  } 
  Double_t GetParticleFraction(Int_t i, Double_t p);

 private:
  AliProtonAnalysis(const AliProtonAnalysis&); // Not implemented
  AliProtonAnalysis& operator=(const AliProtonAnalysis&); // Not implemented

  Bool_t IsAccepted(AliESDtrack *track);
  Float_t GetSigmaToVertex(AliESDtrack* esdTrack); 
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
  Bool_t fMinTPCClustersFlag, fMinITSClustersFlag; //shows if this cut is used or not
  Bool_t fMaxChi2PerTPCClusterFlag, fMaxChi2PerITSClusterFlag; //shows if this cut is used or not
  Bool_t fMaxCov11Flag, fMaxCov22Flag, fMaxCov33Flag, fMaxCov44Flag, fMaxCov55Flag; //shows if this cut is used or not
  Bool_t fMaxSigmaToVertexFlag; //shows if this cut is used or not
  Bool_t fITSRefitFlag, fTPCRefitFlag; //shows if this cut is used or not
  
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

  TH2F *fHistYPtProtons; //Y-Pt of Protons
  TH2F *fHistYPtAntiProtons; // Y-Pt of Antiprotons
  
  ClassDef(AliProtonAnalysis,0);
};

#endif
