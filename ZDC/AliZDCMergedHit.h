#ifndef ALIZDCMERGEDHIT_H
#define ALIZDCMERGEDHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  Hits classes for set ZDC                  //
////////////////////////////////////////////////

#include "TObject.h"

class AliZDCMergedHit : public TObject {

public:
  AliZDCMergedHit() {}
  AliZDCMergedHit(Int_t *sector, Float_t *mhits);
  AliZDCMergedHit(AliZDCMergedHit* oldhit) {*this=*oldhit;}
  virtual ~AliZDCMergedHit() {}

  // Getters 
  virtual Int_t   GetSector(Int_t i) {return fSector[i];}
  virtual Float_t GetPrimKinEn()     {return fPrimKinEn;}
  virtual Float_t GetXImpact()       {return fXImpact;}
  virtual Float_t GetYImpact()       {return fYImpact;}
  virtual Float_t GetSFlag()         {return fSFlag;}
  virtual Float_t GetLightPMQ()      {return fLightPMQ;}
  virtual Float_t GetLightPMC()      {return fLightPMC;}
  virtual Float_t GetEnergy()        {return fEnergy;}


  // Data members
  Int_t      fSector[2];    //Array of volumes
  Float_t    fPrimKinEn;    //Primary particle energy
  Float_t    fXImpact;      //x-coord. of the impact point over the ZDC
  Float_t    fYImpact;      //y-coord. of the impact point over the ZDC
  Float_t    fSFlag;        //Secondary flag
  Float_t    fLightPMQ;     //Cerenkov light produced in each quadrant
  Float_t    fLightPMC;     //Cerenkov light seen by the common PM
  Float_t    fEnergy;       //Total energy deposited in eV
 

  // Operators
  Int_t operator == (AliZDCMergedHit &quad) {
     Int_t i;
     for(i=0; i<2; i++) if(fSector[i]!=quad.GetSector(i)) return 0;
     return 1;
  }
  
  virtual AliZDCMergedHit& operator + (AliZDCMergedHit &quad) {
     fLightPMQ+=quad.GetLightPMQ();
     fLightPMC+=quad.GetLightPMC();
     fEnergy+=quad.GetEnergy();
     return *this;
  }

  // Print method
  virtual void Print(Option_t *) const {
     printf(" ### MergedHit: sector[0] =  %d sector[1] =  %d "
	    "  LightPMQ = %f, LightPMC = %f,  Deposited E = %f\n ", 
	    fSector[0],fSector[1],fLightPMQ,fLightPMC,fEnergy);
  }

  ClassDef(AliZDCMergedHit,1)  // MergedHits for the Zero Degree Calorimeters
};
 
#endif
