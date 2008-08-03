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
#include "TH1I.h"
#include "TList.h"

#include "AliPID.h"

class TF1;
class TH2F;
class TH1D;

class AliAODEvent;
class AliAODtrack;
class AliESDEvent;
class AliESDtrack;
class AliExternalTrackParam;
class AliStack;

class AliProtonAnalysis : public TObject {
 public:
  AliProtonAnalysis();
  AliProtonAnalysis(Int_t nbinsY, Float_t fLowY, Float_t fHighY,
		    Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
  virtual ~AliProtonAnalysis();

  void UseTPCOnly() {fUseTPCOnly = kTRUE;}
  
  void InitAnalysisHistograms(Int_t nbinsY, Float_t fLowY, Float_t fHighY,
			      Int_t nbinsPt, Float_t fLowPt, Float_t fHighPt);
  Bool_t ReadFromFile(const char* filename);
  void Analyze(AliESDEvent *fESD);
  void Analyze(AliAODEvent *fAOD);
  void Analyze(AliStack *stack);
  
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

  TH1I *GetEvenHtistogram() {return fHistEvents;}

  Int_t   GetNumberOfAnalyzedEvents() {return (Int_t)fHistEvents->GetEntries();} 
  Bool_t  PrintMean(TH1 *hist, Double_t edge);
  Bool_t  PrintYields(TH1 *hist, Double_t edge); 

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
  void    SetMaxCov11(Double_t maxCov11) {fMaxCov11 = maxCov11; fMaxCov11Flag = kTRUE;}
  void    SetMaxCov22(Double_t maxCov22) {fMaxCov22 = maxCov22; fMaxCov22Flag = kTRUE;}
  void    SetMaxCov33(Double_t maxCov33) {fMaxCov33 = maxCov33; fMaxCov33Flag = kTRUE;}
  void    SetMaxCov44(Double_t maxCov44) {fMaxCov44 = maxCov44; fMaxCov44Flag = kTRUE;}
  void    SetMaxCov55(Double_t maxCov55) {fMaxCov55 = maxCov55; fMaxCov55Flag = kTRUE;}
  void    SetMaxSigmaToVertex(Double_t maxSigmaToVertex) {
    fMaxSigmaToVertex = maxSigmaToVertex;
    fMaxSigmaToVertexFlag = kTRUE;
  }
  void    SetMaxSigmaToVertexTPC(Double_t maxSigmaToVertex) {
    fMaxSigmaToVertexTPC = maxSigmaToVertex;
    fMaxSigmaToVertexTPCFlag = kTRUE;
  }
  void    SetITSRefit() {fITSRefitFlag = kTRUE;}
  void    SetTPCRefit() {fTPCRefitFlag = kTRUE;}
  void    SetESDpid() {fESDpidFlag = kTRUE;}
  void    SetTPCpid() {fTPCpidFlag = kTRUE;}

  //QA histograms
  void SetQAOn() {
    fQAHistograms = kTRUE;
    fGlobalQAList = new TList();
    fQA2DList = new TList();
    fQAPrimaryProtonsAcceptedList = new TList();
    fQAPrimaryProtonsRejectedList = new TList();
    fQASecondaryProtonsAcceptedList = new TList();
    fQASecondaryProtonsRejectedList = new TList();
    fQAPrimaryAntiProtonsAcceptedList = new TList();
    fQAPrimaryAntiProtonsRejectedList = new TList();
    fQASecondaryAntiProtonsAcceptedList = new TList();
    fQASecondaryAntiProtonsRejectedList = new TList();
  }
  void SetQAYPtBins(Int_t nbinsY, Double_t minY, Double_t maxY,
		    Int_t nbinsPt, Double_t minPt, Double_t maxPt) {
      fNBinsY = nbinsY;
      fMinY = minY; fMaxY = maxY;
      fNBinsPt = nbinsPt;
      fMinPt = minPt; fMaxPt = maxPt;
    }
  void InitQA();
  void RunQA(AliStack *stack, AliESDEvent *esd);

  //Prior probabilities
  void SetPriorProbabilities(Double_t *partFrac) {for(Int_t i = 0; i < AliPID::kSPECIESN; i++) fPartFrac[i] = partFrac[i];} 
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
  Bool_t ReadCorrectionContainer(const char* filename);
  TList *GetCorrectionList2D() {return fCorrectionList2D;} 
  TList *GetEfficiencyList1D() {return fEfficiencyList1D;} 
  TList *GetCorrectionList1D() {return fCorrectionList1D;} 

 private:
  AliProtonAnalysis(const AliProtonAnalysis&); // Not implemented
  AliProtonAnalysis& operator=(const AliProtonAnalysis&); // Not implemented

  Bool_t IsAccepted(AliESDtrack *track);
  Bool_t IsAccepted(AliESDtrack *track, AliStack *stack);
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
  Double_t fMaxSigmaToVertexTPC; //max sigma to vertex cut
  Bool_t fMinTPCClustersFlag, fMinITSClustersFlag; //shows if this cut is used or not
  Bool_t fMaxChi2PerTPCClusterFlag, fMaxChi2PerITSClusterFlag; //shows if this cut is used or not
  Bool_t fMaxCov11Flag, fMaxCov22Flag, fMaxCov33Flag, fMaxCov44Flag, fMaxCov55Flag; //shows if this cut is used or not
  Bool_t fMaxSigmaToVertexFlag; //shows if this cut is used or not
  Bool_t fMaxSigmaToVertexTPCFlag; //shows if this cut is used or not
  Bool_t fITSRefitFlag, fTPCRefitFlag; //shows if this cut is used or not
  Bool_t fESDpidFlag, fTPCpidFlag; //shows if this cut is used or not
  
  //QA histograms
  Bool_t fQAHistograms; //Boolean to activate the QA histograms
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
  //primary protons
  /*TH1F *fPrimaryProtonsTPCClustersReject;        //QA histogram for the primaries rejected by the TPC cluster cut
  TH1F *fPrimaryProtonsTPCClustersPass;          //QA histogram for the primaries accepted by the TPC cluster cut
  TH1F *fPrimaryProtonsChi2PerClusterTPCReject;  //QA histogram for the primaries rejected by the chi2 per TPC cluster cut 
  TH1F *fPrimaryProtonsChi2PerClusterTPCPass;    //QA histogram for the primaries accepted by the chi2 per TPC cluster cut
  TH1F *fPrimaryProtonsExtCov11Reject;           //QA histogram for the primaries rejected by the sigma of the local Y cut
  TH1F *fPrimaryProtonsExtCov11Pass;             //QA histogram for the primaries accepted by the sigma of the local Y cut
  TH1F *fPrimaryProtonsExtCov22Reject;           //QA histogram for the primaries rejected by the sigma of the local Z cut 
  TH1F *fPrimaryProtonsExtCov22Pass;             //QA histogram for the primaries accepted by the sigma of the local Z cut 
  TH1F *fPrimaryProtonsExtCov33Reject;           //QA histogram for the primaries rejected by the sigma of the sin(phi) cut 
  TH1F *fPrimaryProtonsExtCov33Pass;             //QA histogram for the primaries accepted by the sigma of the sin(phi) cut 
  TH1F *fPrimaryProtonsExtCov44Reject;           //QA histogram for the primaries rejected by the sigma of the tan(lambda) cut 
  TH1F *fPrimaryProtonsExtCov44Pass;             //QA histogram for the primaries accepted by the sigma of the tan(lambda) cut 
  TH1F *fPrimaryProtonsExtCov55Reject;           //QA histogram for the primaries rejected by the the sigma of 1/pT cut 
  TH1F *fPrimaryProtonsExtCov55Pass;             //QA histogram for the primaries accepted by the sigma of the 1/pT cut 
  TH1F *fPrimaryProtonsSigmaToVertexReject;      //QA histogram for the primaries rejected by the sigma to vertex cut 
  TH1F *fPrimaryProtonsSigmaToVertexPass;        //QA histogram for the primaries accepted by the sigma to vertex cut 
  TH1F *fPrimaryProtonsSigmaToVertexTPCReject;   //QA histogram for the primaries rejected by the sigma to vertex (TPC) cut 
  TH1F *fPrimaryProtonsSigmaToVertexTPCPass;     //QA histogram for the primaries accepted by the sigma to vertex (TPC) cut 
  TH1F *fPrimaryProtonsITSRefitReject;           //QA histogram for the primaries rejected by the ITS refit cut 
  TH1F *fPrimaryProtonsITSRefitPass;             //QA histogram for the primaries accepted by the ITS refit cut
  TH1F *fPrimaryProtonsTPCRefitReject;           //QA histogram for the primaries rejected by the TPC refit cut
  TH1F *fPrimaryProtonsTPCRefitPass;             //QA histogram for the primaries accepted by the TPC refit cut
  TH1F *fPrimaryProtonsESDpidReject;             //QA histogram for the primaries rejected by the ESD pid cut
  TH1F *fPrimaryProtonsESDpidPass;               //QA histogram for the primaries accepted by the ESD pid cut
  TH1F *fPrimaryProtonsTPCpidReject;             //QA histogram for the primaries rejected by the TPC pid cut
  TH1F *fPrimaryProtonsTPCpidPass;               //QA histogram for the primaries accepted by the TPC pid cut
  //secondary protons
  TH1F *fSecondaryProtonsTPCClustersReject;        //QA histogram for the secondaries rejected by the TPC cluster cut
  TH1F *fSecondaryProtonsTPCClustersPass;          //QA histogram for the secondaries accepted by the TPC cluster cut
  TH1F *fSecondaryProtonsChi2PerClusterTPCReject;  //QA histogram for the secondaries rejected by the chi2 per TPC cluster cut 
  TH1F *fSecondaryProtonsChi2PerClusterTPCPass;    //QA histogram for the secondaries accepted by the chi2 per TPC cluster cut
  TH1F *fSecondaryProtonsExtCov11Reject;           //QA histogram for the secondaries rejected by the sigma of the local Y cut
  TH1F *fSecondaryProtonsExtCov11Pass;             //QA histogram for the secondaries accepted by the sigma of the local Y cut
  TH1F *fSecondaryProtonsExtCov22Reject;           //QA histogram for the secondaries rejected by the sigma of the local Z cut 
  TH1F *fSecondaryProtonsExtCov22Pass;             //QA histogram for the secondaries accepted by the sigma of the local Z cut 
  TH1F *fSecondaryProtonsExtCov33Reject;           //QA histogram for the secondaries rejected by the sigma of the sin(phi) cut 
  TH1F *fSecondaryProtonsExtCov33Pass;             //QA histogram for the secondaries accepted by the sigma of the sin(phi) cut 
  TH1F *fSecondaryProtonsExtCov44Reject;           //QA histogram for the secondaries rejected by the sigma of the tan(lambda) cut 
  TH1F *fSecondaryProtonsExtCov44Pass;             //QA histogram for the secondaries accepted by the sigma of the tan(lambda) cut 
  TH1F *fSecondaryProtonsExtCov55Reject;           //QA histogram for the secondaries rejected by the the sigma of 1/pT cut 
  TH1F *fSecondaryProtonsExtCov55Pass;             //QA histogram for the secondaries accepted by the sigma of the 1/pT cut 
  TH1F *fSecondaryProtonsSigmaToVertexReject;      //QA histogram for the secondaries rejected by the sigma to vertex cut 
  TH1F *fSecondaryProtonsSigmaToVertexPass;        //QA histogram for the secondaries accepted by the sigma to vertex cut 
  TH1F *fSecondaryProtonsSigmaToVertexTPCReject;   //QA histogram for the secondaries rejected by the sigma to vertex (TPC) cut 
  TH1F *fSecondaryProtonsSigmaToVertexTPCPass;     //QA histogram for the secondaries accepted by the sigma to vertex (TPC) cut 
  TH1F *fSecondaryProtonsITSRefitReject;           //QA histogram for the secondaries rejected by the ITS refit cut 
  TH1F *fSecondaryProtonsITSRefitPass;             //QA histogram for the secondaries accepted by the ITS refit cut
  TH1F *fSecondaryProtonsTPCRefitReject;           //QA histogram for the secondaries rejected by the TPC refit cut
  TH1F *fSecondaryProtonsTPCRefitPass;             //QA histogram for the secondaries accepted by the TPC refit cut
  TH1F *fSecondaryProtonsESDpidReject;             //QA histogram for the secondaries rejected by the ESD pid cut
  TH1F *fSecondaryProtonsESDpidPass;               //QA histogram for the secondaries accepted by the ESD pid cut
  TH1F *fSecondaryProtonsTPCpidReject;             //QA histogram for the secondaries rejected by the TPC pid cut
  TH1F *fSecondaryProtonsTPCpidPass;               //QA histogram for the secondaries accepted by the TPC pid cut
  //primary antiprotons
  TH1F *fPrimaryAntiProtonsTPCClustersReject;        //QA histogram for the primaries rejected by the TPC cluster cut
  TH1F *fPrimaryAntiProtonsTPCClustersPass;          //QA histogram for the primaries accepted by the TPC cluster cut
  TH1F *fPrimaryAntiProtonsChi2PerClusterTPCReject;  //QA histogram for the primaries rejected by the chi2 per TPC cluster cut 
  TH1F *fPrimaryAntiProtonsChi2PerClusterTPCPass;    //QA histogram for the primaries accepted by the chi2 per TPC cluster cut
  TH1F *fPrimaryAntiProtonsExtCov11Reject;           //QA histogram for the primaries rejected by the sigma of the local Y cut
  TH1F *fPrimaryAntiProtonsExtCov11Pass;             //QA histogram for the primaries accepted by the sigma of the local Y cut
  TH1F *fPrimaryAntiProtonsExtCov22Reject;           //QA histogram for the primaries rejected by the sigma of the local Z cut 
  TH1F *fPrimaryAntiProtonsExtCov22Pass;             //QA histogram for the primaries accepted by the sigma of the local Z cut 
  TH1F *fPrimaryAntiProtonsExtCov33Reject;           //QA histogram for the primaries rejected by the sigma of the sin(phi) cut 
  TH1F *fPrimaryAntiProtonsExtCov33Pass;             //QA histogram for the primaries accepted by the sigma of the sin(phi) cut 
  TH1F *fPrimaryAntiProtonsExtCov44Reject;           //QA histogram for the primaries rejected by the sigma of the tan(lambda) cut 
  TH1F *fPrimaryAntiProtonsExtCov44Pass;             //QA histogram for the primaries accepted by the sigma of the tan(lambda) cut 
  TH1F *fPrimaryAntiProtonsExtCov55Reject;           //QA histogram for the primaries rejected by the the sigma of 1/pT cut 
  TH1F *fPrimaryAntiProtonsExtCov55Pass;             //QA histogram for the primaries accepted by the sigma of the 1/pT cut 
  TH1F *fPrimaryAntiProtonsSigmaToVertexReject;      //QA histogram for the primaries rejected by the sigma to vertex cut 
  TH1F *fPrimaryAntiProtonsSigmaToVertexPass;        //QA histogram for the primaries accepted by the sigma to vertex cut 
  TH1F *fPrimaryAntiProtonsSigmaToVertexTPCReject;   //QA histogram for the primaries rejected by the sigma to vertex (TPC) cut 
  TH1F *fPrimaryAntiProtonsSigmaToVertexTPCPass;     //QA histogram for the primaries accepted by the sigma to vertex (TPC) cut 
  TH1F *fPrimaryAntiProtonsITSRefitReject;           //QA histogram for the primaries rejected by the ITS refit cut 
  TH1F *fPrimaryAntiProtonsITSRefitPass;             //QA histogram for the primaries accepted by the ITS refit cut
  TH1F *fPrimaryAntiProtonsTPCRefitReject;           //QA histogram for the primaries rejected by the TPC refit cut
  TH1F *fPrimaryAntiProtonsTPCRefitPass;             //QA histogram for the primaries accepted by the TPC refit cut
  TH1F *fPrimaryAntiProtonsESDpidReject;             //QA histogram for the primaries rejected by the ESD pid cut
  TH1F *fPrimaryAntiProtonsESDpidPass;               //QA histogram for the primaries accepted by the ESD pid cut
  TH1F *fPrimaryAntiProtonsTPCpidReject;             //QA histogram for the primaries rejected by the TPC pid cut
  TH1F *fPrimaryAntiProtonsTPCpidPass;               //QA histogram for the primaries accepted by the TPC pid cut
  //secondary antiprotons
  TH1F *fSecondaryAntiProtonsTPCClustersReject;        //QA histogram for the secondaries rejected by the TPC cluster cut
  TH1F *fSecondaryAntiProtonsTPCClustersPass;          //QA histogram for the secondaries accepted by the TPC cluster cut
  TH1F *fSecondaryAntiProtonsChi2PerClusterTPCReject;  //QA histogram for the secondaries rejected by the chi2 per TPC cluster cut 
  TH1F *fSecondaryAntiProtonsChi2PerClusterTPCPass;    //QA histogram for the secondaries accepted by the chi2 per TPC cluster cut
  TH1F *fSecondaryAntiProtonsExtCov11Reject;           //QA histogram for the secondaries rejected by the sigma of the local Y cut
  TH1F *fSecondaryAntiProtonsExtCov11Pass;             //QA histogram for the secondaries accepted by the sigma of the local Y cut
  TH1F *fSecondaryAntiProtonsExtCov22Reject;           //QA histogram for the secondaries rejected by the sigma of the local Z cut 
  TH1F *fSecondaryAntiProtonsExtCov22Pass;             //QA histogram for the secondaries accepted by the sigma of the local Z cut 
  TH1F *fSecondaryAntiProtonsExtCov33Reject;           //QA histogram for the secondaries rejected by the sigma of the sin(phi) cut 
  TH1F *fSecondaryAntiProtonsExtCov33Pass;             //QA histogram for the secondaries accepted by the sigma of the sin(phi) cut 
  TH1F *fSecondaryAntiProtonsExtCov44Reject;           //QA histogram for the secondaries rejected by the sigma of the tan(lambda) cut 
  TH1F *fSecondaryAntiProtonsExtCov44Pass;             //QA histogram for the secondaries accepted by the sigma of the tan(lambda) cut 
  TH1F *fSecondaryAntiProtonsExtCov55Reject;           //QA histogram for the secondaries rejected by the the sigma of 1/pT cut 
  TH1F *fSecondaryAntiProtonsExtCov55Pass;             //QA histogram for the secondaries accepted by the sigma of the 1/pT cut 
  TH1F *fSecondaryAntiProtonsSigmaToVertexReject;      //QA histogram for the secondaries rejected by the sigma to vertex cut 
  TH1F *fSecondaryAntiProtonsSigmaToVertexPass;        //QA histogram for the secondaries accepted by the sigma to vertex cut 
  TH1F *fSecondaryAntiProtonsSigmaToVertexTPCReject;   //QA histogram for the secondaries rejected by the sigma to vertex (TPC) cut 
  TH1F *fSecondaryAntiProtonsSigmaToVertexTPCPass;     //QA histogram for the secondaries accepted by the sigma to vertex (TPC) cut 
  TH1F *fSecondaryAntiProtonsITSRefitReject;           //QA histogram for the secondaries rejected by the ITS refit cut 
  TH1F *fSecondaryAntiProtonsITSRefitPass;             //QA histogram for the secondaries accepted by the ITS refit cut
  TH1F *fSecondaryAntiProtonsTPCRefitReject;           //QA histogram for the secondaries rejected by the TPC refit cut
  TH1F *fSecondaryAntiProtonsTPCRefitPass;             //QA histogram for the secondaries accepted by the TPC refit cut
  TH1F *fSecondaryAntiProtonsESDpidReject;             //QA histogram for the secondaries rejected by the ESD pid cut
  TH1F *fSecondaryAntiProtonsESDpidPass;               //QA histogram for the secondaries accepted by the ESD pid cut
  TH1F *fSecondaryAntiProtonsTPCpidReject;             //QA histogram for the secondaries rejected by the TPC pid cut
  TH1F *fSecondaryAntiProtonsTPCpidPass;*/               //QA histogram for the secondaries accepted by the TPC pid cut

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

  TH1I *fHistEvents; //event counter
  TH2F *fHistYPtProtons; //Y-Pt of Protons
  TH2F *fHistYPtAntiProtons; // Y-Pt of Antiprotons

  //Corrections
  TList *fCorrectionList2D; //list for the 2d corrections 
  TList *fEfficiencyList1D; //list for the 1d efficiencies
  TList *fCorrectionList1D; //list for the 1d corrections 
  
  ClassDef(AliProtonAnalysis,0);
};

#endif
