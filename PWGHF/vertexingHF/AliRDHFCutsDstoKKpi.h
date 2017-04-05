#ifndef ALIRDHFCUTSDSTOKKPI_H
#define ALIRDHFCUTSDSTOKKPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
/// \class Class AliRDHFCutsDstoKKpi
/// \brief class for cuts on AOD reconstructed Ds->KKpi
/// \author Author: A.Dainese, andrea.dainese@pd.infn.it
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsDstoKKpi : public AliRDHFCuts 
{
 public:

  AliRDHFCutsDstoKKpi(const char* name="CutsDstoKKpi");
  
  virtual ~AliRDHFCutsDstoKKpi(){}

  AliRDHFCutsDstoKKpi(const AliRDHFCutsDstoKKpi& source);
  AliRDHFCutsDstoKKpi& operator=(const AliRDHFCutsDstoKKpi& source); 
 
  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters){
    return GetCutVarsForOpt(d,vars,nvars,pdgdaughters,0x0);
  }
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters,AliAODEvent *aod);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel){
    return IsSelected(obj,selectionLevel,0x0);
  }
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent* aod);


  virtual Int_t IsSelectedPID(AliAODRecoDecayHF *rd);
  virtual Int_t IsSelectedPIDBayes(AliAODRecoDecayHF *rd);
  virtual void SetStandardCutsPP2010();
   
  virtual Bool_t IsInFiducialAcceptance(Double_t pt,Double_t y) const;
  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(11,iPtBin)] : 1.e6);}
  UInt_t GetPIDTrackTPCTOFBitMap(AliAODTrack *track) const;

  void Setd0MeasMinusExpCut(Int_t nPtBins, Float_t *cutval);
  void Setd0Cut(Int_t nPtBins, Float_t *cutval);
    
  enum TrackPIDBit{kTPCPionLess1,kTPCPionMore1Less2,kTPCPionMore2Less3,kTPCPionMore3,
                   kTPCKaonLess1,kTPCKaonMore1Less2,kTPCKaonMore2Less3,kTPCKaonMore3,
                   kTPCProtonLess1,kTPCProtonMore1Less2,kTPCProtonMore2Less3,kTPCProtonMore3,
                   kTOFPionLess1,kTOFPionMore1Less2,kTOFPionMore2Less3,kTOFPionMore3,
                   kTOFKaonLess1,kTOFKaonMore1Less2,kTOFKaonMore2Less3,kTOFKaonMore3,
                   kTOFProtonLess1,kTOFProtonMore1Less2,kTOFProtonMore2Less3,kTOFProtonMore3};
                
                 
  enum EDsPid {kConservative, kStrong, kStrongPDep, kBayesianMaxProb, kBayesianThreshold, kBayesianWeights};
  void SetPidOption(Int_t opt){
    fPidOption=opt;
  }
  void SetUseBayesianPIDWithMaxProb(Double_t dist=0.01){
    fPidOption=kBayesianMaxProb;
    fDistToMaxProb=dist;
    fPidHF->SetUseCombined(kTRUE);
  }

  void SetUseBayesianPIDWithThresholds(Double_t thr=0.05){
    fPidOption=kBayesianThreshold;
    fBayesThreshold=thr;
    fPidHF->SetUseCombined(kTRUE);
  }
  void SetUseBayesianPIDWithWeights(){
    fPidOption=kBayesianWeights;
    fPidHF->SetUseCombined(kTRUE);
  }
  void SetMaxPtStrongPid(Float_t spid){fMaxPtStrongPid=spid;}
  void SetMaxPStrongPidK(Float_t spid){fMaxPStrongPidK=spid;}
  void SetMaxPStrongPidpi(Float_t spid){fMaxPStrongPidpi=spid;}
  
  void SetUseReferencePhiMass(Double_t value){fUseRefPhiMass = kTRUE; fPhiMassRef = value;}
    
  Int_t GetPidOption() const {return fPidOption;}
  Double_t GetWeightForKKpi() const {return fWeightKKpi;}
  Double_t GetWeightForpiKK() const {return fWeightpiKK;}

  Bool_t IsCutOnResonancesApplied() const {return fCutOnResonances;}
  void  ApplyCutOnResonances(Bool_t opt=kTRUE){
    fCutOnResonances=opt;
  }


 protected:
 
  Bool_t fCutOnResonances;  /// switch for the cuts on phi and K0* inv. mass
  Int_t fPidOption;         /// pid option
  Float_t fMaxPtStrongPid; /// Maximum pt of candidate to apply strong Pid p dependent
  Float_t fMaxPStrongPidK; /// Maximum P of track to apply strong Pid on K
  Float_t fMaxPStrongPidpi; /// Maximum P of track to apply strong Pid on pi
  Double_t fDistToMaxProb; /// Difference between max probability
  Double_t fBayesThreshold;/// Threshold for Bayesian PID probability
  Double_t fWeightKKpi; /// weight for KKpi for kBayesianWeights
  Double_t fWeightpiKK; /// weight for piKK for kBayesianWeights
  Double_t fPhiMassRef; /// Reference Phi mass to be used for the cut on delta phi mass (instead of PDG value)
  Bool_t fUseRefPhiMass; ///swicth to the usage of Reference Phi mass (instead of PDG value)
  Bool_t fUsed0MeasMinusExpCut; /// switch for cut on d0meas-d0exp
  Float_t* fMaxd0MeasMinusExp;  //[fnPtBins] cut values on d0meas-d0exp
  Bool_t fUsed0Cut; /// switch for cut on d0
  Float_t* fMaxd0;  //[fnPtBins] cut values on d0

  /// \cond CLASSIMP     
  ClassDef(AliRDHFCutsDstoKKpi,6);  /// class for cuts on AOD reconstructed Ds->KKpi
  /// \endcond
};

#endif
