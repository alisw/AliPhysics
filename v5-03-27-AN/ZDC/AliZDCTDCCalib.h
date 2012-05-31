#ifndef ALIZDCTDCCALIB_H
#define ALIZDCTDCCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////
//  class for ZDC calibration -> TDC mean values  //
////////////////////////////////////////////////////

#include "TNamed.h"
#include "AliCDBEntry.h"

class AliZDC;

class AliZDCTDCCalib: public TNamed {

 public:
  AliZDCTDCCalib();
  AliZDCTDCCalib(const char* name);
  AliZDCTDCCalib(const AliZDCTDCCalib &calibda);
  AliZDCTDCCalib& operator= (const AliZDCTDCCalib &calibda);
  virtual ~AliZDCTDCCalib();
  void Reset();
  virtual void  Print(Option_t *) const; 

  Float_t GetMeanTDC(Int_t ch) const 
  	{if(ch<6) return fMeanTDC[ch];
	 else return 0.;};
  Float_t GetWidthTDC(Int_t ch) const 
  	{if(ch<6) return fWidthTDC[ch];
	 else return 0.;};
	 
  void SetMeanTDC(Int_t ch, Float_t val) {fMeanTDC[ch]=val;}
  void SetWidthTDC(Int_t ch, Float_t val) {fWidthTDC[ch]=val;}
  void SetMeanTDC(Float_t* mean);
  void SetWidthTDC(Float_t* width);
  
 protected:
  // --- Pedestals
  Float_t  fMeanTDC[6];	 	// Mean TDC values 
  Float_t  fWidthTDC[6];	// TDC widths 
  //
  ClassDef(AliZDCTDCCalib,1)    // ZDC TDC calibration data
};

#endif
