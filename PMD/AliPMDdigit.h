#ifndef ALIPMDDIGIT_H
#define ALIPMDDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store digits  for PMD                              //
//                                                     //
//-----------------------------------------------------//

#include "TObject.h"
class TClonesArray;

class AliPMDdigit : public TObject
{
 public:
  AliPMDdigit();
  AliPMDdigit(Int_t trnumber, Int_t det, Int_t smnumber,
	      Int_t cellnumber, Float_t adc);
  AliPMDdigit(AliPMDdigit *pmddigit) {*this = *pmddigit;}
  AliPMDdigit (const AliPMDdigit &pmddigit);  // copy constructor
  AliPMDdigit &operator=(const AliPMDdigit &pmddigit); // assignment op

  virtual ~AliPMDdigit();

  Int_t   GetTrackNumber() const;
  Int_t   GetDetector() const;
  Int_t   GetSMNumber() const;
  Int_t   GetCellNumber() const;
  Float_t GetADC() const;

 protected:
  Int_t   fTrNumber;    // Parent Track number
  Int_t   fDet;         // Detecor Number (0:PRE, 1:CPV)
  Int_t   fSMNumber;    // Serial Module Number
  Int_t   fCellNumber;  // Cell Number (row(0-47)*96 + col(0-95))
  Float_t fADC;         // Energy deposition(ADC) in a hexagonal cell
  
  ClassDef(AliPMDdigit,2) // Digits object for Detector set:PMD
};

#endif
