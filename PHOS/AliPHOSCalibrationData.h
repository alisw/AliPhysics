#ifndef ALIPHOSCALIBRATIONDATA_H
#define ALIPHOSCALIBRATIONDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Base Class for PHOS     
//  for Calibration data                 
//*-- Author:D.Peressounko


// --- ROOT system ---
#include "TArrayF.h"
#include "TString.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "TObject.h"

class AliPHOSCalibrationData : public TObject {

public:
  AliPHOSCalibrationData() ;          // Default ctor (empty)
  AliPHOSCalibrationData(const char* category, const char * version="v1" ,Int_t nchanels = 17920) ; 
  AliPHOSCalibrationData(const AliPHOSCalibrationData & cd) ;
  virtual ~AliPHOSCalibrationData() ;

  virtual const char* GetSubsystem(void)const{return "PHOS" ;}
  virtual const char* GetVersion(void)  const{return fVersion ;}
  virtual const char* GetCategory(void) const {return fCategory ;} 
  virtual void  GetValidityRange(Int_t &begin,Int_t &end) const {begin=fBegin;end=fEnd ;}
  
  Float_t Data(Int_t channel)const ;
  Float_t DataCheck(Int_t channel) const ;
  Int_t   NChannels(void){if(fData) return fData->GetSize() ;
                          else return 0 ;}

  void SetData(Int_t channel,Float_t data){fData->AddAt(data,channel); }
  void SetData(TArrayF &array){if(fData) delete fData; fData=new TArrayF(array) ;} ;
  void SetDataCheck(Int_t channel,Float_t check){fDataCheck->AddAt(check,channel) ;}
  void SetDataCheck(TArrayF &array){if(fData) delete fDataCheck; fDataCheck=new TArrayF(array) ;} ;
  void SetValidityRange(Int_t begin,Int_t end){fBegin=begin;fEnd=end;}

  AliPHOSCalibrationData & operator = (const AliPHOSCalibrationData & rvalue) ;

private:
  TString fCategory; //e.g. Gains, Pedestals,...
  TString fVersion ; //Version (title)
  TArrayF * fData ;      //Data themself
  TArrayF * fDataCheck ; //Parameter to check Data validity (e.g. width of pedestal peak)
  Int_t   fBegin ; // validity period
  Int_t   fEnd ;   // validity period

  ClassDef(AliPHOSCalibrationData,1)  // description 

};

#endif // AliPHOSCALIBRATIONDATA_H
