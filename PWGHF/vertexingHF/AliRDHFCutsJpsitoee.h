#ifndef ALIRDHFCUTSJPSITOEE_H
#define ALIRDHFCUTSJPSITOEE_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
/// \class Class AliRDHFCutsJpsitoee
/// \brief class for cuts on AOD reconstructed J/psi->ee (from B)
/// \author Author: A.Dainese, andrea.dainese@pd.infn.it
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsJpsitoee : public AliRDHFCuts 
{
 public:

  AliRDHFCutsJpsitoee(const char* name="CutsJpsitoee");
  
  virtual ~AliRDHFCutsJpsitoee(){}

  AliRDHFCutsJpsitoee(const AliRDHFCutsJpsitoee& source);
  AliRDHFCutsJpsitoee& operator=(const AliRDHFCutsJpsitoee& source); 
 
  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel) 
                         {return IsSelected(obj,selectionLevel,0);}
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent* aod);
  
  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(1,iPtBin)] : 1.e6);}

 protected:


  /// \cond CLASSIMP
  ClassDef(AliRDHFCutsJpsitoee,1);  /// class for cuts on AOD reconstructed J/psi->ee
  /// \endcond
};

#endif
