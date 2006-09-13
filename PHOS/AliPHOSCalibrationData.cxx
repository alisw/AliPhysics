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
// Calibration data class contains two arrays: 
// data themselves and their quality (or possibility to check 
// goodness of the data. 
// There are two kinds of data in PHOS "Gains" and "Pedestals"
// each of them have as well title and validity range to distinguish
// different calibrations.
//
//*-- Author : D.Peressounko
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSCalibrationData.h"


ClassImp(AliPHOSCalibrationData)


//____________________________________________________________________________ 
AliPHOSCalibrationData::AliPHOSCalibrationData():
  fCategory(), 
  fVersion(),
  fData(0),
  fDataCheck(0),
  fBegin(0),
  fEnd(0)
{
  // default ctor
}

//____________________________________________________________________________ 
AliPHOSCalibrationData::AliPHOSCalibrationData(const char * category, const char * version, Int_t nchannels):
  fCategory(category), 
  fVersion(version),
  fData(new TArrayF(nchannels)),
  fDataCheck(new TArrayF(nchannels)),
  fBegin(0),
  fEnd(0)
{
  // ctor: sets up the calibration IO
}

//____________________________________________________________________________ 
AliPHOSCalibrationData::AliPHOSCalibrationData(const AliPHOSCalibrationData & cd):
  TObject(cd),
  fCategory(cd.fCategory),
  fVersion(cd.fVersion),
  fData(new TArrayF(*cd.fData)),
  fDataCheck(new TArrayF(*cd.fDataCheck)),
  fBegin(cd.fBegin),
  fEnd(cd.fEnd)
{
  //copy ctor
}

//____________________________________________________________________________ 
  AliPHOSCalibrationData::~AliPHOSCalibrationData()
{
  // dtor: deletes the arrays
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
Float_t AliPHOSCalibrationData::Data(Int_t channel)const 
{
  // returns calibration data for a given channel   
  return fData->At(channel) ;
}
//____________________________________________________________________________ 
Float_t AliPHOSCalibrationData::DataCheck(Int_t channel)const 
{
  // returns calibration data check for a given channel   
  return fDataCheck->At(channel) ;
}
//____________________________________________________________________________ 
AliPHOSCalibrationData & AliPHOSCalibrationData::operator = (const AliPHOSCalibrationData & rvalue)
{ 
  // overload of =
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
