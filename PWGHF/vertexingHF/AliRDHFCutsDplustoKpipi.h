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
class AliAODPidHF;

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
  virtual void SetStandardCutsPbPb2011();

  Int_t GetPIDBitMask(AliAODRecoDecayHF *rd);
  UInt_t GetPIDTrackTPCTOFBitMap(AliAODTrack *track) const;
  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(1,iPtBin)] : 1.e6);}
  void SetUseStrongPid(Int_t spid){fUseStrongPid=spid;}
  void SetMaxPtStrongPid(Float_t spid){fMaxPtStrongPid=spid;}
  void SetMaxPStrongPidK(Float_t spid){fMaxPStrongPidK=spid;}
  void SetMaxPStrongPidpi(Float_t spid){fMaxPStrongPidpi=spid;}
  Int_t GetStrongPid() const {return fUseStrongPid;}
  Float_t GetMaxPtStrongPid() const {return fMaxPtStrongPid;}
  Float_t GetMaxPtStrongPidK() const {return fMaxPStrongPidK;}
  Float_t GetMaxPtStrongPidpi() const {return fMaxPStrongPidpi;}
  void SetUseImpParProdCorrCut(Bool_t use){
    fUseImpParProdCorrCut=use;
  }
  Bool_t GetUseImpParProdCorrCut() const {
    return fUseImpParProdCorrCut;
  }
  
  enum TrackPIDBit{kTPCPionLess1,kTPCPionMore1Less2,kTPCPionMore2Less3,kTPCPionMore3,
                   kTPCKaonLess1,kTPCKaonMore1Less2,kTPCKaonMore2Less3,kTPCKaonMore3,
                   kTPCProtonLess1,kTPCProtonMore1Less2,kTPCProtonMore2Less3,kTPCProtonMore3,
                   kTOFPionLess1,kTOFPionMore1Less2,kTOFPionMore2Less3,kTOFPionMore3,
                   kTOFKaonLess1,kTOFKaonMore1Less2,kTOFKaonMore2Less3,kTOFKaonMore3,
                   kTOFProtonLess1,kTOFProtonMore1Less2,kTOFProtonMore2Less3,kTOFProtonMore3};
 protected:

 private:
  Int_t fUseStrongPid; //use strong pid 0 no,1 only for K,2 pi 3 both
  Float_t fMaxPtStrongPid;//Maximum pt of candidate to apply strong Pid
  Float_t fMaxPStrongPidK;//Maximum P of track to apply strong Pid on K
  Float_t fMaxPStrongPidpi;//Maximum P of track to apply strong Pid on pi
  Bool_t fUseImpParProdCorrCut; //switch for d0K*d0pi1 vs. d0K*d0pi2 cut

  ClassDef(AliRDHFCutsDplustoKpipi,6);  // class for cuts on AOD reconstructed 
                                   // D+->Kpipi
};

#endif

