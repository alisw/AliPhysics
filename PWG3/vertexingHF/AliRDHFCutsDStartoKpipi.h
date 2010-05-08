#ifndef ALIRDHFCUTSDSTARTOKPIPI_H
#define ALIRDHFCUTSDSTARTOKPIPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliRDHFCutsDStartoKpipi
// class for cuts on AOD reconstructed DStar->Kpipi
// Author: A.Grelli, alessandro.grelli@uu.nl
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsDStartoKpipi : public AliRDHFCuts 
{
 public:

  AliRDHFCutsDStartoKpipi();
  
  virtual ~AliRDHFCutsDStartoKpipi(){}

  AliRDHFCutsDStartoKpipi(const AliRDHFCutsDStartoKpipi& source);
  AliRDHFCutsDStartoKpipi& operator=(const AliRDHFCutsDStartoKpipi& source); 
 
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgD0daughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);
  virtual Int_t IsD0FromDStarSelected(Double_t pt, TObject* obj,Int_t selectionLevel); 
  
  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(1,iPtBin)] : 1.e6);}

 protected:


  ClassDef(AliRDHFCutsDStartoKpipi,1);  // class for cuts on AOD reconstructed D0->Kpipi
};

#endif
