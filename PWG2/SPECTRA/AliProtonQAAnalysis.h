#ifndef ALIPROTONQAANALYSIS_H
#define ALIPROTONQAANALYSIS_H

/*  See cxx source for full Copyright notice */


/* $Id: AliProtonQAAnalysis.h 29114 2008-10-03 16:49:02Z pchrist $ */

//-------------------------------------------------------------------------
//                       Class AliProtonQAAnalysis
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

class AliESDEvent;
class AliESDtrack;
class AliStack;

class AliProtonQAAnalysis : public TObject {
 public:
  AliProtonQAAnalysis();
  virtual ~AliProtonQAAnalysis();

  void UseTPCOnly() {fUseTPCOnly = kTRUE;}
  
  void InitAnalysisHistograms(Int_t nbinsY, Float_t fLowY, Float_t fHighY,
			      Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);

  //Cut functions
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

  //QA histograms
  void SetQAOn();
  void SetQAYPtBins(Int_t nbinsY, Double_t minY, Double_t maxY,
		    Int_t nbinsPt, Double_t minPt, Double_t maxPt);
  void InitQA();
  void RunQA(AliStack *stack, AliESDEvent *esd);
  TList *GetGlobalQAList() {return fGlobalQAList;}

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

 private:
  AliProtonQAAnalysis(const AliProtonQAAnalysis&); // Not implemented
  AliProtonQAAnalysis& operator=(const AliProtonQAAnalysis&); // Not implemented

  Bool_t   IsAccepted(AliESDtrack *track);
  void     FillQA(AliESDtrack *track, AliStack *stack);
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
  
  //QA histograms
  //Bool_t fQAHistograms; //Boolean to activate the QA histograms
  TList *fGlobalQAList; //TList storing the directories for the QA histograms
  TList *fQA2DList; //TList storing the accepted primary/secondary (anti)protons
  TList *fQAPrimaryProtonsAcceptedList; //list of the QA histos for accepted primary protons
  TList *fQAPrimaryProtonsRejectedList; //list of the QA histos for rejected primary protons
  TList *fQASecondaryProtonsAcceptedList; //list of the QA histos for accepted secondary protons
  TList *fQASecondaryProtonsRejectedList; //list of the QA histos for rejected secondary protons
  TList *fQAPrimaryAntiProtonsAcceptedList; //list of the QA histos for accepted primary antiprotons
  TList *fQAPrimaryAntiProtonsRejectedList; //list of the QA histos for rejected primary antiprotons
  TList *fQASecondaryAntiProtonsAcceptedList; //list of the QA histos for accepted secondary antiprotons
  TList *fQASecondaryAntiProtonsRejectedList; //list of the QA histos for rejected secondary antiprotons

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

  ClassDef(AliProtonQAAnalysis,0);
};

#endif
