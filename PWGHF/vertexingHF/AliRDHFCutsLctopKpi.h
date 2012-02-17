#ifndef ALIRDHFCUTSLCTOPKPI_H
#define ALIRDHFCUTSLCTOPKPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliRDHFCutsLctopKpi
// class for cuts on AOD reconstructed Lc->pKpi
// Author: A.Dainese, andrea.dainese@pd.infn.it
//***********************************************************

#include "AliRDHFCuts.h"
#include "AliAODPidHF.h"
#include "AliAODRecoDecayHF3Prong.h"

class AliRDHFCutsLctopKpi : public AliRDHFCuts 
{
 public:

 enum EPIDStrategy {
  kNSigma,
  kCombined
 };
 enum ECutsStrategy {
  kStandard,
  kKF
 };


  AliRDHFCutsLctopKpi(const char* name="CutsLctopKpi");
  
  virtual ~AliRDHFCutsLctopKpi();

  AliRDHFCutsLctopKpi(const AliRDHFCutsLctopKpi& source);
  AliRDHFCutsLctopKpi& operator=(const AliRDHFCutsLctopKpi& source); 
 
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


  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel)
                           {return IsSelected(obj,selectionLevel,0);}
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent *aod);
  using AliRDHFCuts::IsSelectedPID;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPID(AliAODRecoDecayHF* obj);
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
 protected:
  AliAODPidHF *fPidObjprot;
  AliAODPidHF *fPidObjpion;
  Bool_t fUseImpParProdCorrCut; //switch for cut on d0p*d0K vs. d0K*d0pi 

private:
  EPIDStrategy fPIDStrategy;                // PIS strategy (nsigma, combined)
  Double_t fPIDThreshold[AliPID::kSPECIES]; // PID threshold for each species
  ECutsStrategy fCutsStrategy;              // cut strategy (standard or KF)

  ClassDef(AliRDHFCutsLctopKpi,5);  // class for cuts on AOD reconstructed Lc->pKpi
};

#endif

