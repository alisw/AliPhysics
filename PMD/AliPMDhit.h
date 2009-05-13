#ifndef ALIPMDHIT_H
#define ALIPMDHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//              hits classes for set:PMD      //
////////////////////////////////////////////////
 
#include "AliHit.h"
#include "Riostream.h"

class TClonesArray;

class AliPMDhit : public AliHit {

 public:
  AliPMDhit();
  AliPMDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits);
  AliPMDhit(AliPMDhit* oldhit);
  virtual ~AliPMDhit() {}
  Int_t GetVolume(Int_t i) const {return fVolume[i];}
  Float_t GetEnergy() const {return fEnergy;}
  int operator == (AliPMDhit &cell) const;
  AliPMDhit operator + (AliPMDhit &cell) {
    fEnergy+=cell.GetEnergy();
    return *this;
  }
  void Print(Option_t *) const {
    printf("PMD Cell %d %d %d %d %d %d\n   Primary %d -   Energy %f\n",
	   fVolume[0],fVolume[1],fVolume[2],fVolume[3],
	   fVolume[4],fVolume[5],fTrack,fEnergy);
  }
  
 protected:
  Int_t      fVolume[6];  //array of volumes
  Float_t    fEnergy;      //Total energy deposited in eV
  
  ClassDef(AliPMDhit,5)  //Hits object for set:PMD
};
#endif
