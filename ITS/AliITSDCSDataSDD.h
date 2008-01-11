#ifndef ALIITSDCSDATASDD_H
#define ALIITSDCSDATASDD_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to define object containing SDD DCS data                //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include<TObject.h>
#include<TArrayF.h>
#include<TArrayI.h>

class AliITSDCSDataSDD : public TObject { 

 public:
  AliITSDCSDataSDD();
  AliITSDCSDataSDD(Int_t npts);
  ~AliITSDCSDataSDD(){};
  void SetNPoints(Int_t npts);
  void SetValues(Int_t time, Float_t field, Float_t templ, Float_t tempr);
  void Compress();

  Int_t GetNumberOfValues() const {return fSetPoints;}
  Int_t GetTimeStamp(Int_t i) const {return fTimeStamp.At(i);}
  Float_t GetDriftField(Int_t i) const {return fDriftField.At(i);}
  Float_t GetLeftTemperature(Int_t i) const {return fTemperatureLeft.At(i);}
  Float_t GetRightTemperature(Int_t i) const {return fTemperatureRight.At(i);}
  void PrintValues() const;

 private:


  Int_t fNpts;   // number of values for DCS data points
  Int_t fSetPoints; // number of set values
  TArrayI fTimeStamp; // DCS time stamp
  TArrayF fDriftField; // drift field (calculated from HV and MV)
  TArrayF fTemperatureLeft; // temperature from sensor on left hybrid
  TArrayF fTemperatureRight; // temperature from sensor on right hybrid

  ClassDef(AliITSDCSDataSDD, 1);
};

#endif
