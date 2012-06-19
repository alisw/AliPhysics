#ifndef ALIRDHFCUTSLCTOV0_H
#define ALIRDHFCUTSLCTOV0_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliRDHFCutsLctoV0
// class for cuts on AOD reconstructed Lc-> V0 + bachelor
//***********************************************************

#include "AliAODPidHF.h"
#include "AliRDHFCuts.h"

class AliRDHFCutsLctoV0 : public AliRDHFCuts 
{
 public:

  enum {
    kLcToK0Spr=0x0001,
    kLcToLBarpi=0x0002,
    kLcToLpi=0x0004
  };

  AliRDHFCutsLctoV0(const char* name="CutsLctoV0", Short_t v0channel=0);
  
  virtual ~AliRDHFCutsLctoV0();

  AliRDHFCutsLctoV0(const AliRDHFCutsLctoV0& source);
  AliRDHFCutsLctoV0& operator=(const AliRDHFCutsLctoV0& source); 
 
  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);
  
  using AliRDHFCuts::IsSelectedPID;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);
  Int_t IsSelected(TObject* obj, Int_t selectionLevel, Int_t cutIndex);

  Int_t CombinePIDCuts (Int_t returnvalue, Int_t returnvaluePID) const;

  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(7,iPtBin)] : 1.e6);}

  void SetPidV0pos(AliAODPidHF* pidV0pos) {
    if (fPidHFV0pos) delete fPidHFV0pos;
    fPidHFV0pos = new AliAODPidHF(*pidV0pos);
  }
  void SetPidV0neg(AliAODPidHF* pidV0neg) {
    if (fPidHFV0neg) delete fPidHFV0neg;
    fPidHFV0neg = new AliAODPidHF(*pidV0neg);
  }

  AliAODPidHF * GetPidV0pos() { return fPidHFV0pos; }
  AliAODPidHF * GetPidV0neg() { return fPidHFV0neg; }

 protected:

 private:

  AliAODPidHF *fPidHFV0pos;
  AliAODPidHF *fPidHFV0neg;
  //UShort_t fV0channel;

  ClassDef(AliRDHFCutsLctoV0,2);  // class for cuts on AOD reconstructed Lc->V0+bachelor
};

#endif
