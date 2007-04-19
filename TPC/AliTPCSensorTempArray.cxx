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
//  TPC calibration class for parameters which saved per pad                 //
//  Authors: Marian Ivanov and Haavard Helstrup                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCSensorTempArray.h"

ClassImp(AliTPCSensorTempArray)

const char kFname[] = "TempSensor.txt";
const char kAmandaString[] = "tpc_temp:PT_%d.Temperature";
const Int_t  kMinGraph = 10;       // minimum #points of graph to be fitted
const Int_t  kMinPoints = 10;      // minimum number of points per knot in fit
const Int_t  kIter = 10;           // number of iterations for spline fit
const Double_t  kMaxDelta = 0.00;  // precision parameter for spline fit
const Int_t  kFitReq = 2;          // fit requirement, 2 = continuous 2nd derivative

//_____________________________________________________________________________
AliTPCSensorTempArray::AliTPCSensorTempArray():AliDCSSensorArray(),
 fAmandaString(kAmandaString)
{
  //
  // AliTPCSensorTempArray default constructor
  //
 
}
//_____________________________________________________________________________
AliTPCSensorTempArray::AliTPCSensorTempArray(Int_t prevRun) : 
                AliDCSSensorArray(prevRun,"TPC/Calib/Temperature"),
 fAmandaString(kAmandaString)
{
}
//_____________________________________________________________________________
AliTPCSensorTempArray::AliTPCSensorTempArray(UInt_t startTime, UInt_t endTime,
                       const char *filepath)
             :AliDCSSensorArray(),
     fAmandaString(kAmandaString)
{
  //
  // AliTPCSensorTempArray default constructor
  //
  char *expPath = gSystem->ExpandPathName(filepath);
  TString filename(expPath);
  filename.Append('/');
  filename.Append(kFname);
  fSensors =  AliTPCSensorTemp::ReadListInd(filename.Data(),fFirstSensor,fLastSensor);
  fStartTime = TTimeStamp(startTime);
  fEndTime   = TTimeStamp(endTime);
  delete expPath;
}

//_____________________________________________________________________________
AliTPCSensorTempArray::AliTPCSensorTempArray(const char *fname) : 
                                                  AliDCSSensorArray(),
 fAmandaString(kAmandaString)
{
  //
  // AliTPCSensorTempArray constructor
  //
  fSensors = AliTPCSensorTemp::ReadListInd(fname,fFirstSensor,fLastSensor);
  fSensors->BypassStreamer(kFALSE);
}


//_____________________________________________________________________________
AliTPCSensorTempArray::AliTPCSensorTempArray(const AliTPCSensorTempArray &c):
  AliDCSSensorArray(c),
  fAmandaString(c.fAmandaString)
{
  //
  // AliTPCSensorTempArray copy constructor
  //

}

///_____________________________________________________________________________
AliTPCSensorTempArray::~AliTPCSensorTempArray()
{
  //
  // AliTPCSensorTempArray destructor
  //
}

//_____________________________________________________________________________
AliTPCSensorTempArray &AliTPCSensorTempArray::operator=(const AliTPCSensorTempArray &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTPCSensorTempArray &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTPCSensorTempArray::Copy(TObject &c) const
{
  //
  // Copy function
  //

  TObject::Copy(c);
}
//_____________________________________________________________________________
void AliTPCSensorTempArray::ReadSensors(const char *fname) 
{
  //
  // Read list of temperature sensors from text file
  //
  fSensors = AliTPCSensorTemp::ReadListInd(fname,fFirstSensor,fLastSensor);
}  
//_____________________________________________________________________________
void AliTPCSensorTempArray::SetGraph(TMap *map) 
{
  // 
  // Read graphs from DCS maps 
  //
  AliDCSSensorArray::SetGraph(map,fAmandaString.Data());
}  
//_____________________________________________________________________________
void AliTPCSensorTempArray::MakeSplineFit(TMap *map) 
{
  // 
  // Make spline fits from DCS maps 
  //
  AliDCSSensorArray::MakeSplineFit(map,fAmandaString.Data());
}  


//_____________________________________________________________________________
TMap* AliTPCSensorTempArray::ExtractDCS(TMap *dcsMap) 
{
 //
 // Extract temperature graphs from DCS maps
 //

 TMap *values = AliDCSSensorArray::ExtractDCS(dcsMap,fAmandaString.Data());
 return values;
}

