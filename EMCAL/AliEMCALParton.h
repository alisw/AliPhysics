#ifndef ALIEMCALPARTON_H
#define ALIEMCALPARTON_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

#include <TObject.h>
//*-- Author: Renan Cabrera (Creighton U.)


class AliEMCALParton : public TObject {
public:
  AliEMCALParton();
  AliEMCALParton(Float_t energy, Float_t phi, Float_t eta);
  virtual ~AliEMCALParton();
  void SetEnergy(Float_t val) {fEnergy = val;}
  void SetPhi(Float_t val)    {fPhi    = val;}
  void SetEta(Float_t val)    {fEta    = val;}
  Float_t Energy()  {return fEnergy;}
  Float_t Phi()     {return fPhi;}
  Float_t Eta()     {return fEta;}
  
protected:
  Float_t  fEnergy;   // Jet Energy
  Float_t  fEta;      // Jet Phi
  Float_t  fPhi;      // Jet Eta
  ClassDef(AliEMCALParton,3) // Jet for EMCAL
    
} ;

#endif // ALIEMCALParton_H
