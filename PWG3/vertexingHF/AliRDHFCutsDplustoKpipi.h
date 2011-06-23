#ifndef ALIRDHFCUTSDPLUSTOKPIPI_H
#define ALIRDHFCUTSDPLUSTOKPIPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 


/////////////////////////////////////////////////////////////
//
// Class for cuts on AOD reconstructed D+->Kpipi
//
// Author: R. Bala bala@to.infn.it
//         G. Ortona ortona@to.infn.it
/////////////////////////////////////////////////////////////


#include "AliRDHFCuts.h"
#include "AliAODPidHF.h"
class AliRDHFCutsDplustoKpipi : public AliRDHFCuts 
{
 public:

  AliRDHFCutsDplustoKpipi(const char* name="CutsDplustoKpipi");
  
  virtual ~AliRDHFCutsDplustoKpipi(){};
  AliRDHFCutsDplustoKpipi(const AliRDHFCutsDplustoKpipi& source);
  AliRDHFCutsDplustoKpipi& operator=(const AliRDHFCutsDplustoKpipi& source); 

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

  virtual Bool_t IsInFiducialAcceptance(Double_t pt,Double_t y) const;
  virtual void SetStandardCutsPP2010();
  virtual void SetStandardCutsPbPb2010();

  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(1,iPtBin)] : 1.e6);}
  void SetUseStrongPid(Int_t spid){fUseStrongPid=spid;}
  void SetMaxPtStrongPid(Float_t spid){fMaxPtStrongPid=spid;}
  Int_t GetStrongPid() const {return fUseStrongPid;}
  Float_t GetMaxPtStrongPid(){return fMaxPtStrongPid;}
  void SetUseImpParProdCorrCut(Bool_t use){
    fUseImpParProdCorrCut=use;
  }
  Bool_t GetUseImpParProdCorrCut() const {
    return fUseImpParProdCorrCut;
  }
  
 protected:

 private:
  Int_t fUseStrongPid; //use strong pid 0 no,1 only for K,2 both pi and K
  Float_t fMaxPtStrongPid;//Maximum pt to apply strong Pid
  Bool_t fUseImpParProdCorrCut; //switch for d0K*d0pi1 vs. d0K*d0pi2 cut

  ClassDef(AliRDHFCutsDplustoKpipi,4);  // class for cuts on AOD reconstructed 
                                   // D+->Kpipi
};

#endif
