/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


//_________________________________________________________________________
// Calibration data and their quality  
//
//*-- Author : D.Peressounko
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSCalibrationData.h"


ClassImp(AliPHOSCalibrationData)


//____________________________________________________________________________ 
  AliPHOSCalibrationData::AliPHOSCalibrationData():TObject() {
  fBegin=0;
  fEnd=0;
  fData = 0 ;
  fDataCheck = 0 ;
  fCategory=""; 
  fVersion="" ; 
}

//____________________________________________________________________________ 
  AliPHOSCalibrationData::AliPHOSCalibrationData(const char * category, const char * version, Int_t nchannels){
  fData      = new TArrayF(nchannels) ;
  fDataCheck = new TArrayF(nchannels) ;
  fCategory=category; 
  fVersion=version ; 
}
//____________________________________________________________________________ 
AliPHOSCalibrationData::AliPHOSCalibrationData(const AliPHOSCalibrationData & cd){
  fData = new TArrayF(*cd.fData) ;
  fDataCheck = new TArrayF(*cd.fDataCheck) ;
  fCategory=cd.fCategory; 
  fVersion=cd.fVersion ; 
}
//____________________________________________________________________________ 
  AliPHOSCalibrationData::~AliPHOSCalibrationData()
{
  if(fData){
    delete fData ;
    fData=0 ;
  }
  if(fDataCheck){
    delete fDataCheck ;
    fDataCheck=0;
  }
}
//____________________________________________________________________________ 
Float_t AliPHOSCalibrationData::Data(Int_t channel)const {
  return fData->At(channel) ;
}
//____________________________________________________________________________ 
Float_t AliPHOSCalibrationData::DataCheck(Int_t channel)const {
  return fDataCheck->At(channel) ;
}
//____________________________________________________________________________ 
AliPHOSCalibrationData & AliPHOSCalibrationData::operator = (const AliPHOSCalibrationData & rvalue){
  if(fData)
    delete fData; 
  fData      = new TArrayF(*rvalue.fData) ;
  if(fDataCheck)
    delete fDataCheck ;
  fDataCheck = new TArrayF(*rvalue.fDataCheck) ;
  fCategory=rvalue.fCategory; 
  fVersion=rvalue.fVersion ; 
  fBegin=rvalue.fBegin ;
  fEnd = rvalue.fEnd ;
  return *this ;
}
