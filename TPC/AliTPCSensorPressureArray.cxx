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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TPC calibration class for pressure sensors                               //
//  Authors: Marian Ivanov and Haavard Helstrup                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCSensorPressureArray.h"

ClassImp(AliTPCSensorPressureArray)

const char kFname[] = "PressureSensor.txt";
const char kAmandaString[] = "system:gcs-ALITPC-ALITPC-AnalogInput-0%d.ProcessInput.PosSt";

//_____________________________________________________________________________
AliTPCSensorPressureArray::AliTPCSensorPressureArray():AliDCSSensorArray(),
 fAmandaString(kAmandaString)
{
  //
  // AliTPCSensorPressureArray default constructor
  //
 
}
//_____________________________________________________________________________
AliTPCSensorPressureArray::AliTPCSensorPressureArray(Int_t prevRun) : 
                AliDCSSensorArray(prevRun,"TPC/Calib/Pressure"),
 fAmandaString(kAmandaString)
{
}
//_____________________________________________________________________________
AliTPCSensorPressureArray::AliTPCSensorPressureArray(UInt_t startTime, UInt_t endTime,
                       const char *filepath)
             :AliDCSSensorArray(),
     fAmandaString(kAmandaString)
{
  //
  // AliTPCSensorPressureArray default constructor
  //
  char *expPath = gSystem->ExpandPathName(filepath);
  TString filename(expPath);
  filename.Append('/');
  filename.Append(kFname);
  fSensors =  AliTPCSensorPressure::ReadList(filename.Data());
  fStartTime = TTimeStamp(startTime);
  fEndTime   = TTimeStamp(endTime);
  delete expPath;
}

//_____________________________________________________________________________
AliTPCSensorPressureArray::AliTPCSensorPressureArray(const char *fname) : 
                                                  AliDCSSensorArray(),
 fAmandaString(kAmandaString)
{
  //
  // AliTPCSensorPressureArray constructor
  //
  fSensors = AliTPCSensorPressure::ReadList(fname);
  fSensors->BypassStreamer(kFALSE);
}


//_____________________________________________________________________________
AliTPCSensorPressureArray::AliTPCSensorPressureArray(const AliTPCSensorPressureArray &c):
  AliDCSSensorArray(c),
  fAmandaString(c.fAmandaString)
{
  //
  // AliTPCSensorPressureArray copy constructor
  //

}

///_____________________________________________________________________________
AliTPCSensorPressureArray::~AliTPCSensorPressureArray()
{
  //
  // AliTPCSensorPressureArray destructor
  //
}

//_____________________________________________________________________________
AliTPCSensorPressureArray &AliTPCSensorPressureArray::operator=(const AliTPCSensorPressureArray &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTPCSensorPressureArray &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTPCSensorPressureArray::Copy(TObject &c) const
{
  //
  // Copy function
  //

  TObject::Copy(c);
}
//_____________________________________________________________________________
void AliTPCSensorPressureArray::ReadSensors(const char *fname) 
{
  //
  // Read list of temperature sensors from text file
  //
  fSensors = AliTPCSensorPressure::ReadList(fname);
}  
//_____________________________________________________________________________
void AliTPCSensorPressureArray::SetGraph(TMap *map) 
{
  // 
  // Read graphs from DCS maps 
  //
  AliDCSSensorArray::SetGraph(map,fAmandaString.Data());
}  
//_____________________________________________________________________________
void AliTPCSensorPressureArray::MakeSplineFit(TMap *map) 
{
  // 
  // Make spline fits from DCS maps 
  //
  AliDCSSensorArray::MakeSplineFit(map,fAmandaString.Data());
}  


//_____________________________________________________________________________
TMap* AliTPCSensorPressureArray::ExtractDCS(TMap *dcsMap) 
{
 //
 // Extract temperature graphs from DCS maps
 //

 TMap *values = AliDCSSensorArray::ExtractDCS(dcsMap,fAmandaString.Data());
 return values;
}

//_____________________________________________________________________________
AliTPCSensorPressure* AliTPCSensorPressureArray::GetSensor(Int_t type, Int_t side, Int_t sector, Int_t num) 
{
 //
 //  Return sensor information for sensor specified by type, side, sector and num
 //
 Int_t nsensors = fSensors->GetEntries();
 for (Int_t isensor=0; isensor<nsensors; isensor++) {
   AliTPCSensorPressure *entry = (AliTPCSensorPressure*)fSensors->At(isensor);
   if (entry->GetSide() == side &&
       entry->GetType() == type &&
       entry->GetSector() == sector &&
       entry->GetNum() == num ) return entry;
 }
 return 0;
}
//_____________________________________________________________________________
AliTPCSensorPressure* AliTPCSensorPressureArray::GetSensor(Int_t IdDCS){
  return dynamic_cast<AliTPCSensorPressure*>(AliDCSSensorArray::GetSensor(IdDCS));
}
//_____________________________________________________________________________
AliTPCSensorPressure* AliTPCSensorPressureArray::GetSensor(Double_t x, Double_t y, Double_t z){
  return dynamic_cast<AliTPCSensorPressure*>(AliDCSSensorArray::GetSensor(x,y,z));
}
