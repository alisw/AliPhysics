#ifndef ALIRDHFCUTSXICTOPKPI_H
#define ALIRDHFCUTSXICTOPKPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
/// \class Class AliRDHFCutsXictopKpi
/// \brief class for cuts on AOD reconstructed Xic->pKpi
/// \First implementation by copying the AliRDHFCutsLctopKpi class
/// \author Author: M. Faggin (mfaggin@cern.ch)
//***********************************************************

#include "AliRDHFCuts.h"
#include "AliAODPidHF.h"
#include "AliAODRecoDecayHF3Prong.h"
//#include "AliVertexerTracks.h"

class AliRDHFCutsXictopKpi : public AliRDHFCuts 
{
 public:

 enum EPIDStrategy {
  kNSigma,
  kCombined,
  kCombinedSoft,
  kNSigmaStrong,
  kCombinedpPb,
  kCombinedpPb2,
  kNSigmaPbPb,
  kNSigmaMin,
  kCombinedProb
 };
 enum ECutsStrategy {
  kStandard,
  kKF
 };

  AliRDHFCutsXictopKpi(const char* name="CutsLctopKpi");
  
  virtual ~AliRDHFCutsXictopKpi();

  AliRDHFCutsXictopKpi(const AliRDHFCutsXictopKpi& source);
  AliRDHFCutsXictopKpi& operator=(const AliRDHFCutsXictopKpi& source); 
 
  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters){
    return GetCutVarsForOpt(d,vars,nvars,pdgdaughters,0x0);
  }
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters,AliAODEvent *aod);

  void SetPidpion(AliAODPidHF* pidPion) { 
      if(fPidObjpion) delete fPidObjpion;
      fPidObjpion=new AliAODPidHF(*pidPion);
      }
  void SetPidprot(AliAODPidHF* pidProt) {
      if(fPidObjprot) delete fPidObjprot;
      fPidObjprot=new AliAODPidHF(*pidProt);
  }

  virtual void SetStandardCutsPP2010();
  virtual void SetStandardCutsPbPb2010();
  virtual void SetStandardCutsPbPb2011();
  virtual void SetStandardCutsPPb2013(); 


  AliAODPidHF* GetPidpion() const {return fPidObjpion;}
  AliAODPidHF* GetPidprot() const {return fPidObjprot;}
  void SetPIDStrategy(EPIDStrategy pidStrategy) {
   fPIDStrategy=pidStrategy;
  }
  EPIDStrategy GetPIDStrategy() const {
   return fPIDStrategy;
  }
  void SetCutsStrategy(ECutsStrategy cutsStrategy) {
   fCutsStrategy=cutsStrategy;
  }
  ECutsStrategy GetCutsStrategy() const {
   return fCutsStrategy;
  }
  void SetPIDThreshold(AliPID::EParticleType species,Double_t threshold) {
   fPIDThreshold[static_cast<Int_t>(species)]=threshold;
  }
 Double_t GetPIDThreshold(AliPID::EParticleType species) const {
  return fPIDThreshold[static_cast<Int_t>(species)];
 }
 Bool_t GetUseSpecialCut(){return fUseSpecialCut;}
 void SetUseSpecialCut(Bool_t useSpecialCut=kTRUE){fUseSpecialCut=useSpecialCut;}

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel)
                           {return IsSelected(obj,selectionLevel,0);}
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent *aod);
  using AliRDHFCuts::IsSelectedPID;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPID(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPIDSoft(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPIDpPb(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPIDpPb2(AliAODRecoDecayHF* obj);
  Int_t IsSelectedPIDStrong(AliAODRecoDecayHF* obj);
  Int_t IsSelectedNSigmaPbPb(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPIDProb(AliAODRecoDecayHF* obj);
  Int_t CombinePIDCuts (Int_t returnvalue, Int_t returnvaluePID) const;

  virtual Bool_t IsInFiducialAcceptance(Double_t pt,Double_t y) const;
  
  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(11,iPtBin)] : 1.e6);}

  void SetUseImpParProdCorrCut(Bool_t use){
    fUseImpParProdCorrCut=use;
  }
  Bool_t GetUseImpParProdCorrCut() const {
    return fUseImpParProdCorrCut;
  }

  AliKFParticle* ReconstructKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field,Bool_t constraint) const;

  // function for PID cuts exploration
  void ExplorePID(AliPIDResponse* pid_resp, AliAODRecoDecayHF3Prong* cand, UInt_t PIDcase, Bool_t &is_pKpi_passed, Bool_t &is_piKp_passed);

 protected:
  AliAODPidHF *fPidObjprot;
  AliAODPidHF *fPidObjpion;
  Bool_t fUseImpParProdCorrCut; /// switch for cut on d0p*d0K vs. d0K*d0pi

  // --- functions for PID cuts exploration ---
  // Home-made functions for PID
  //
  // NB: functions tuned on TH2 PID plots from pp @ 5 TeV!
  //
  Float_t func_TPCprot_up(Float_t x)
  {
      Float_t value;
      if(x<4)     value = 3.9-1.938*x+0.585*x*x-0.0764*x*x*x+0.003*x*x*x*x;
      else        value = 3.9-1.938*4+0.585*4*4-0.0764*4*4*4+0.003*4*4*4*4;
      return TMath::Min(value,(Float_t) 3.);
  }
  Float_t func_TPCprot_down(Float_t x)
  {
      Float_t value;
      if(x<3)     value = 5.6-57.89*x+123.78*x*x-116.83*x*x*x+55.44*x*x*x*x-13*x*x*x*x*x+1.2*x*x*x*x*x*x;
      else        value = -2.375;
      return TMath::Max(value,(Float_t) -3.);
  }
  Float_t func_TPCpion_up(Float_t x)
  {
      Float_t value;
      value = 4.23-3.85*x+2.6*x*x-0.59*x*x*x+0.044*x*x*x*x;
      return TMath::Min(value,(Float_t) 3.);
  }
  Float_t func_TPCpion_down(Float_t x)
  {
      Float_t value;
      if(x<1.4)   value = -2-16.33*x+20.93*x*x-6.63*x*x*x;
      else        value = -1.63636;
      return TMath::Max(value,(Float_t) -3.);
  }
  Float_t func_TPCkaon_up(Float_t x)
  {
      Float_t value;
      if(x<2.1)   value = 2.24+3.05*x-4.47*x*x+1.37*x*x*x+9.16e-8;
      else        value = 1.625;
      return TMath::Min(value,(Float_t) 3.);
  }
  Float_t func_TPCkaon_down(Float_t x)
  {
      Float_t value;
      if(x<2)   value = -9.04 + 31.41*x - 43.81*x*x+23.83*x*x*x-4.44*x*x*x*x;
      else      value = -1.7875;
      return TMath::Max(value,(Float_t) -3.);
  }
  Float_t func_TOFprot_up(Float_t x)
  {
      Float_t value;
      if(x<6)     value = 1.45+2.11*x-1.19*x*x+0.24*x*x*x-0.016*x*x*x*x;
      else        value = 1.45+2.11*6-1.19*6*6+0.24*6*6*6-0.016*6*6*6*6;
      return TMath::Min(value,(Float_t) 3.);
  }
  Float_t func_TOFprot_down(Float_t x)
  {
      Float_t value;
      if(x<5)     value = 0.019-6.93*x+5.15*x*x-1.61*x*x*x+0.24*x*x*x*x-0.014*x*x*x*x*x;
      else        value = 0.019-6.93*5+5.15*5*5-1.61*5*5*5+0.24*5*5*5*5-0.014*5*5*5*5*5;
      return TMath::Max(value,(Float_t) -3.);
  }
  Float_t func_TOFpion_up(Float_t x)
  {
      Float_t value;
      if(x<1.7)                   value = 2.4;
      else if(1.7<x && x<3.42)    value = 7.1848 -4.15*x+0.78*x*x;
      else                        value = 2.18;
      return TMath::Min(value,(Float_t) 3.);
  }
  Float_t func_TOFpion_down(Float_t x)
  {
      Float_t value = -3.;
      return TMath::Max(value,(Float_t) -3.);
  }
  Float_t func_TOFkaon_up(Float_t x)
  {
      Float_t value;
      if(x<6)     value = 1.68+1.19*x-0.644*x*x+0.115*x*x*x-0.0065*x*x*x*x;
      else        value = 1.68+1.19*6-0.644*6*6+0.115*6*6*6-0.0065*6*6*6*6;
      return TMath::Min(value,(Float_t) 3.);
  }
  Float_t func_TOFkaon_down(Float_t x)
  {
      Float_t value;
      if(x<3)     value = -0.64-5*x+3.78*x*x-0.94*x*x*x+0.072*x*x*x*x;
      else        value = -0.64-5*3+3.78*3*3-0.94*3*3*3+0.072*3*3*3*3;
      return TMath::Max(value,(Float_t) -3.);
  }
  //-------------------------


  EPIDStrategy fPIDStrategy;                /// PIS strategy (nsigma, combined)
  Double_t fPIDThreshold[AliPID::kSPECIES]; /// PID threshold for each species
  ECutsStrategy fCutsStrategy;              /// cut strategy (standard or KF)
  Bool_t fUseSpecialCut;

private:

  /// \cond CLASSIMP    
  ClassDef(AliRDHFCutsXictopKpi,1);  /// class for cuts on AOD reconstructed Xic->pKpi
  /// \endcond
};

#endif
