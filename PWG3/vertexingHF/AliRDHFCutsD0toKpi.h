#ifndef ALIRDHFCUTSD0TOKPI_H
#define ALIRDHFCUTSD0TOKPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliRDHFCutsD0toKpi
// class for cuts on AOD reconstructed D0->Kpi
// Author: A.Dainese, andrea.dainese@pd.infn.it
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsD0toKpi : public AliRDHFCuts 
{
 public:

  AliRDHFCutsD0toKpi();
  
  virtual ~AliRDHFCutsD0toKpi(){}

  AliRDHFCutsD0toKpi(const AliRDHFCutsD0toKpi& source);
  AliRDHFCutsD0toKpi& operator=(const AliRDHFCutsD0toKpi& source); 
 
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);

 protected:


  ClassDef(AliRDHFCutsD0toKpi,1);  // class for cuts on AOD reconstructed 
                                   // D0->Kpi
};

#endif
