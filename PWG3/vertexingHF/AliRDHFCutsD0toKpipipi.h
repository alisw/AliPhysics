#ifndef ALIRDHFCUTSD0TOKPIPIPI_H
#define ALIRDHFCUTSD0TOKPIPIPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliRDHFCutsD0toKpipipi
// class for cuts on AOD reconstructed D0->Kpipipi
// Author: A.Dainese, andrea.dainese@pd.infn.it
//	   F.Colamaria, fabio.colamaria@ba.infn.it
//***********************************************************

#include "AliRDHFCuts.h"
#include "AliAODRecoDecayHF4Prong.h"

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
  virtual Int_t IsSelectedFromPID(AliAODRecoDecayHF4Prong *d, Int_t *hyp1, Int_t *hyp2, Int_t *hyp3, Int_t *hyp4);
  virtual Int_t D01Selected(TObject* obj,Int_t selectionLevel);
  virtual Int_t D02Selected(TObject* obj,Int_t selectionLevel);
  virtual Int_t D0bar1Selected(TObject* obj,Int_t selectionLevel);
  virtual Int_t D0bar2Selected(TObject* obj,Int_t selectionLevel);
  
  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(1,iPtBin)] : 1.e6);}
  Bool_t GetUsePID(Int_t iPtBin=0) const { return (GetCuts() ? (Bool_t)(fCutsRD[GetGlobalIndex(8,iPtBin)]) : kFALSE);}

  virtual Bool_t IsInFiducialAcceptance(Double_t pt,Double_t y) const;

 protected:


  ClassDef(AliRDHFCutsD0toKpipipi,1);  // class for cuts on AOD reconstructed D0->Kpipipi
};

#endif
