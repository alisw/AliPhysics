#ifndef ALITRDTRACKINGSECTOR_H
#define ALITRDTRACKINGSECTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackingSector.h,v */

#include <TObject.h>

class AliTRDtimeBin;
class AliTRDgeometry;

class AliTRDtrackingSector : public TObject {

// Provides tools to address clusters which lay within one sector

protected:
  Int_t fN;
  AliTRDtimeBin           *fTimeBin;    // Pointer to array of AliTRDtimeBin
  AliTRDgeometry          *fGeom;       // Pointer to TRD geometry
  Float_t                  fTimeBinSize;  // Time bin size in cm  

public:
  AliTRDtrackingSector() {fN=0; fTimeBin=0; fGeom=0; fTimeBinSize=0;}
  virtual ~AliTRDtrackingSector();
  virtual void SetUp();
 
  AliTRDtimeBin& operator[](Int_t i);
  Int_t GetNtimeBins() const { return fN; }
  Double_t GetX(Int_t l) const;
  Double_t GetMaxY(Int_t l) const; 
  //Double_t GetAlpha() const { return 2*TMath::Pi()/kNsect; } 
  Int_t GetTimeBinNumber(Double_t x) const;
  Int_t GetTimeBin(Int_t det, Int_t local_tb) const;
  Float_t GetPitch() const {return fTimeBinSize;}   

  ClassDef(AliTRDtrackingSector,1)  // Provides tools to address clusters which lay within one sector

}; 

#endif 
