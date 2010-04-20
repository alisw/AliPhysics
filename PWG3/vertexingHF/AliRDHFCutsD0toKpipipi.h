#ifndef ALIRDHFCUTSD0TOKPIPIPI_H
#define ALIRDHFCUTSD0TOKPIPIPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliRDHFCutsD0toKpipipi
// class for cuts on AOD reconstructed D0->Kpipipi
// Author: A.Dainese, andrea.dainese@pd.infn.it
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsD0toKpipipi : public AliRDHFCuts 
{
 public:

  AliRDHFCutsD0toKpipipi(const char* name="CutsD0toKpipipi");
  
  virtual ~AliRDHFCutsD0toKpipipi(){}

  AliRDHFCutsD0toKpipipi(const AliRDHFCutsD0toKpipipi& source);
  AliRDHFCutsD0toKpipipi& operator=(const AliRDHFCutsD0toKpipipi& source); 
 
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);
  
  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(1,iPtBin)] : 1.e6);}

 protected:


  ClassDef(AliRDHFCutsD0toKpipipi,1);  // class for cuts on AOD reconstructed D0->Kpipipi
};

#endif
