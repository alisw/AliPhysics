#ifndef ALIRDHFCUTSLCTOPKPI_H
#define ALIRDHFCUTSLCTOPKPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliRDHFCutsLctopKpi
// class for cuts on AOD reconstructed Lc->pKpi
// Author: A.Dainese, andrea.dainese@pd.infn.it
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsLctopKpi : public AliRDHFCuts 
{
 public:

  AliRDHFCutsLctopKpi(const char* name="CutsLctopKpi");
  
  virtual ~AliRDHFCutsLctopKpi(){}

  AliRDHFCutsLctopKpi(const AliRDHFCutsLctopKpi& source);
  AliRDHFCutsLctopKpi& operator=(const AliRDHFCutsLctopKpi& source); 
 
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);
  
  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(11,iPtBin)] : 1.e6);}

 protected:


  ClassDef(AliRDHFCutsLctopKpi,1);  // class for cuts on AOD reconstructed Lc->pKpi
};

#endif
