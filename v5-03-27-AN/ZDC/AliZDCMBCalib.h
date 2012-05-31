#ifndef ALIZDCMBCALIB_H
#define ALIZDCMBCALIB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////
//  class for ZDC calibration from MB standalone runs  //
/////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TH2F.h"
#include "AliCDBEntry.h"

class AliZDC;

class AliZDCMBCalib: public TNamed {

 public:
  AliZDCMBCalib();
  AliZDCMBCalib(const char* name);
  AliZDCMBCalib(const char* name, 
  	        TH2F *hzdcvszem, TH2F *hzdccvszem, TH2F *hzdcavszem);
  AliZDCMBCalib(const AliZDCMBCalib &calibda);
  AliZDCMBCalib& operator= (const AliZDCMBCalib &calibda);
  virtual ~AliZDCMBCalib();
  void Reset();

  TH2F* GethZDCvsZEM()   const {return fhZDCvsZEM;}  
  TH2F* GethZDCCvsZEM()  const {return fhZDCCvsZEM;}  
  TH2F* GethZDCAvsZEM()  const {return fhZDCAvsZEM;}  
  //
  void SetZDCvsZEM(TH2F *hCorr)  {fhZDCvsZEM  = hCorr;}    
  void SetZDCCvsZEM(TH2F *hCorr) {fhZDCCvsZEM = hCorr;}    
  void SetZDCAvsZEM(TH2F *hCorr) {fhZDCAvsZEM = hCorr;}   

 protected:
  TH2F *  fhZDCvsZEM;	// E_ZDC (total) vs. E_ZEM 
  TH2F *  fhZDCCvsZEM;  // E_ZDC vs. E_ZEM sideC
  TH2F *  fhZDCAvsZEM;  // E_ZDC vs. E_ZEM sideA
  //
  ClassDef(AliZDCMBCalib,1)    // ZDC calibration calibration data
};

#endif
