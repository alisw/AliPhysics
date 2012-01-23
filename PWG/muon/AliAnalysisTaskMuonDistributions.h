#ifndef ALIANALYSISTASKMUONDISTRIBUTIONS_H
#define ALIANALYSISTASKMUONDISTRIBUTIONS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

/* Analysis Task for Muon/Dimuon distributions */

#include "AliAnalysisTaskSE.h"

class TH1D;
class TParticle ;
class TLorentzVector ;
class TFile ;
class AliStack ;
class AliESDtrack;
class AliVParticle;

class AliAnalysisTaskMuonDistributions : public AliAnalysisTaskSE {
  public:

  AliAnalysisTaskMuonDistributions();
  AliAnalysisTaskMuonDistributions(const Char_t* name);
  AliAnalysisTaskMuonDistributions& operator= (const AliAnalysisTaskMuonDistributions& c);
  AliAnalysisTaskMuonDistributions(const AliAnalysisTaskMuonDistributions& c);
  virtual ~AliAnalysisTaskMuonDistributions();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  void UserCreateOutputObjects();
  
  void SetBeamEnergy(Double_t en) {fBeamEnergy=en;}
  void SetAnalysisType(const char* type) {fkAnalysisType=type;}
  void SetInvMassFitLimits(Double_t xmin, Double_t xmax) {fInvMassFitLimitMin=xmin; fInvMassFitLimitMax=xmax;}
  void SetPsiFitLimits(Double_t xmin, Double_t xmax) {fPsiFitLimitMin=xmin; fPsiFitLimitMax=xmax;}
  void SetPsiPFitLimits(Double_t xmin, Double_t xmax) {fPsiPFitLimitMin=xmin; fPsiPFitLimitMax=xmax;}
  void SetBckFitLimits(Double_t xmin, Double_t xmax) {fBckFitLimitMin=xmin; fBckFitLimitMax=xmax;}
  void FitInvariantMassSpectrum(Bool_t massfit=kFALSE) {fInvariantMassFit=massfit;}
 
 protected:
  
  Double_t fBeamEnergy;   // Energy of the beam (required for the CS angle)
  Double_t fInvMassFitLimitMin;  // invariant mass spectrum fit lower limit 
  Double_t fInvMassFitLimitMax;  // invariant mass spectrum fit upper limit 
  Double_t fPsiFitLimitMin;  // psi fit lower limits 
  Double_t fPsiFitLimitMax;  // psi fit upper limits 
  Double_t fPsiPFitLimitMin;  // psi(2S) fit lower limits 
  Double_t fPsiPFitLimitMax;  // psi(2S) fit upper limits 
  Double_t fBckFitLimitMin;  // bck fit lower limits 
  Double_t fBckFitLimitMax;  // bck fit upper limits
  Bool_t fInvariantMassFit; // flag to perform or not inv. mass fit
    
  const char* fkAnalysisType; //ESD or AOD based analysis
  TList *fOutput;  // output file
      
  Float_t InvMass (Float_t e1, Float_t px1, Float_t py1, Float_t pz1, Float_t e2, Float_t px2, Float_t py2, Float_t pz2) const;
  Float_t Rapidity (Float_t e, Float_t pz) const;
  Double_t CostCS (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  Double_t CostHE (Double_t px1, Double_t py1, Double_t pz1, Double_t e1, Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2);
  void FitInvMass(TH1D* histo);
  
  ClassDef(AliAnalysisTaskMuonDistributions,1);
};

#endif
