#ifndef ALIRDHFCUTSLCTOV0_H
#define ALIRDHFCUTSLCTOV0_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliRDHFCutsLctoV0
// class for cuts on AOD reconstructed Lc-> V0 + bachelor
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsLctoV0 : public AliRDHFCuts 
{
 public:

  AliRDHFCutsLctoV0(const char* name="CutsLctoV0");
  
  virtual ~AliRDHFCutsLctoV0(){;}

  AliRDHFCutsLctoV0(const AliRDHFCutsLctoV0& source);
  AliRDHFCutsLctoV0& operator=(const AliRDHFCutsLctoV0& source); 
 
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);
  
  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(7,iPtBin)] : 1.e6);}

 protected:


  ClassDef(AliRDHFCutsLctoV0,1);  // class for cuts on AOD reconstructed Lc->V0+bachelor
};

#endif
