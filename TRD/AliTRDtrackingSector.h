#ifndef ALITRDTRACKINGSECTOR_H
#define ALITRDTRACKINGSECTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackingSector.h,v */

/////////////////////////////////////////////////////////////////////////
//                                                                     //
//  Tracking sector                                                    //
//                                                                     //
/////////////////////////////////////////////////////////////////////////                       

#include <TObject.h>

class AliTRDtimeBin;
class AliTRDgeometry;
class AliTRDparameter;

class AliTRDtrackingSector : public TObject {

// Provides tools to address clusters which lay within one sector

public:

  AliTRDtrackingSector() {fN=0; fTimeBin=0; fGeom=0; fTimeBinSize=0;}
  virtual ~AliTRDtrackingSector();
  virtual void SetUp();
 
  AliTRDtimeBin& operator[](Int_t i);
  Int_t GetNtimeBins() const { return fN; }
  Double_t GetX(Int_t tb) const;
  Int_t   GetTimeBinNumber(Double_t x) const;
  Int_t   GetTimeBin(Int_t det, Int_t local_tb) const;
  Bool_t  TECframe(Int_t tb, Double_t y, Double_t z) const;

  virtual void     SetParameter(AliTRDparameter *par)      { fPar           = par; };
  AliTRDparameter *GetParameter()                    const { return fPar;          };

protected:

  Int_t fN;                                // ???????
  AliTRDgeometry          *fGeom;          // Pointer to TRD geometry
  AliTRDtimeBin           *fTimeBin;       // Pointer to array of AliTRDtimeBin
  Float_t                  fTimeBinSize;   // Time bin size in cm  
  AliTRDparameter         *fPar;           // TRD parameter
					      
  ClassDef(AliTRDtrackingSector,2)  // Provides tools to address clusters which lay within one sector

}; 


#endif 
