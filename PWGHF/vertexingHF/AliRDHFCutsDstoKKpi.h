#ifndef ALIRDHFCUTSDSTOKKPI_H
#define ALIRDHFCUTSDSTOKKPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliRDHFCutsDstoKKpi
// class for cuts on AOD reconstructed Ds->KKpi
// Author: A.Dainese, andrea.dainese@pd.infn.it
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
  virtual void SetStandardCutsPP2010();
   
  virtual Bool_t IsInFiducialAcceptance(Double_t pt,Double_t y) const;
  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(11,iPtBin)] : 1.e6);}
  UInt_t GetPIDTrackTPCTOFBitMap(AliAODTrack *track) const;

  
  enum TrackPIDBit{kTPCPionLess1,kTPCPionMore1Less2,kTPCPionMore2Less3,kTPCPionMore3,
                   kTPCKaonLess1,kTPCKaonMore1Less2,kTPCKaonMore2Less3,kTPCKaonMore3,
                   kTPCProtonLess1,kTPCProtonMore1Less2,kTPCProtonMore2Less3,kTPCProtonMore3,
                   kTOFPionLess1,kTOFPionMore1Less2,kTOFPionMore2Less3,kTOFPionMore3,
                   kTOFKaonLess1,kTOFKaonMore1Less2,kTOFKaonMore2Less3,kTOFKaonMore3,
                   kTOFProtonLess1,kTOFProtonMore1Less2,kTOFProtonMore2Less3,kTOFProtonMore3};
                
                 
  enum EDsPid {kConservative, kStrong};
  void SetPidOption(Int_t opt){
    fPidOption=opt;
  }
  Int_t GetPidOption() const {return fPidOption;}
  
 protected:
  Int_t fPidOption;         //pid option


  ClassDef(AliRDHFCutsDstoKKpi,2);  // class for cuts on AOD reconstructed Ds->KKpi
};

#endif
