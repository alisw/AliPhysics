#ifndef ALIPMDHIT_H
#define ALIPMDHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//              hits classes for set:PMD      //
////////////////////////////////////////////////
 
#include "AliHit.h"
#include <cstdio>

class TClonesArray;

class AliPMDhit : public AliHit {

 public:
  AliPMDhit();
  AliPMDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  AliPMDhit(AliPMDhit* oldhit);
  virtual ~AliPMDhit() {}
  Int_t   GetVolume(Int_t i) const {return fVolume[i];}
  Float_t GetEnergy() const {return fEnergy;}
  Float_t GetTime() const {return fTime;}
  Int_t operator == (AliPMDhit &cell) const;
  AliPMDhit operator + (AliPMDhit &cell) {
    fEnergy+=cell.GetEnergy();
    return *this;
  }
  void Print(Option_t *) const;
 protected:
  Int_t   fVolume[6];    //array of volumes
  Float_t fEnergy;       //Total energy deposited in eV
  Float_t fTime;         //time information for the event (pile-up cal)

  ClassDef(AliPMDhit,6)  //Hits object for set:PMD
};
#endif
