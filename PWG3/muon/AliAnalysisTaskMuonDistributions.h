#ifndef ALIANALYSISTASKMUONDISTRIBUTIONS_H
#define ALIANALYSISTASKMUONDISTRIBUTIONS_H

#include "AliAnalysisTaskSE.h"
#include "TMath.h"

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
  void SetAnalysisType(const char* type) {fAnalysisType=type;}
  void SetInvMassFitLimits(Double_t xmin, Double_t xmax) {fInvMassFitLimitMin=xmin; fInvMassFitLimitMax=xmax;}
  void SetPsiFitLimits(Double_t xmin, Double_t xmax) {fPsiFitLimitMin=xmin; fPsiFitLimitMax=xmax;}
  void SetBckFitLimits(Double_t xmin, Double_t xmax) {fBckFitLimitMin=xmin; fBckFitLimitMax=xmax;}
  void FitInvariantMassSpectrum(Bool_t massfit=kFALSE) {fInvariantMassFit=massfit;}
 
 protected:
  
  Double_t fBeamEnergy;   // Energy of the beam (required for the CS angle)
  Double_t fInvMassFitLimitMin;  // invariant mass spectrum fit limits 
  Double_t fInvMassFitLimitMax;
  Double_t fPsiFitLimitMin;  // psi fit limits 
  Double_t fPsiFitLimitMax;
  Double_t fBckFitLimitMin;  // bck fit limits 
  Double_t fBckFitLimitMax;
  Bool_t fInvariantMassFit; // flag to perform or not inv. mass fit
    
  const char* fAnalysisType; //ESD or AOD based analysis
  TList *fOutput;
      
  Float_t InvMass (Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
  Float_t Rapidity (Float_t, Float_t);
  Double_t CostCS (Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
  Double_t CostHE (Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
  void FitInvMass(TH1D* );
  
  ClassDef(AliAnalysisTaskMuonDistributions,1);
};

#endif
