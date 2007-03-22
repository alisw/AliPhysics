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
const Int_t  kMinGraph = 10;       // minimum #points of graph to be fitted
const Int_t  kMinPoints = 10;      // minimum number of points per knot in fit
const Int_t  kIter = 10;           // number of iterations for spline fit
const Double_t  kMaxDelta = 0.00;  // precision parameter for spline fit
const Int_t  kFitReq = 2;          // fit requirement, 2 = continuous 2nd derivative

//_____________________________________________________________________________
AliTPCSensorTempArray::AliTPCSensorTempArray():AliDCSSensorArray()
{
  //
  // AliTPCSensorTempArray default constructor
  //
  fSensors = 0;
  fFirstSensor = 0;
  fLastSensor = 0;
  TTimeStamp defTime(2000,1,1,0,0,0);
  fStartTime = defTime;
  fEndTime = defTime;

}
//_____________________________________________________________________________
AliTPCSensorTempArray::AliTPCSensorTempArray(Int_t prevRun) : 
                AliDCSSensorArray(prevRun,"TPC/Calib/Temperature")
{
}
//_____________________________________________________________________________
AliTPCSensorTempArray::AliTPCSensorTempArray(UInt_t startTime, UInt_t endTime)
             :AliDCSSensorArray()
{
  //
  // AliTPCSensorTempArray default constructor
  //
  fSensors = AliTPCSensorTemp::ReadListInd(kFname,fFirstSensor,fLastSensor);
  fStartTime = TTimeStamp(startTime);
  fEndTime   = TTimeStamp(endTime);

}

//_____________________________________________________________________________
AliTPCSensorTempArray::AliTPCSensorTempArray(const char *fname) : 
                                                  AliDCSSensorArray()
{
  //
  // AliTPCSensorTempArray constructor
  //
  fSensors = AliTPCSensorTemp::ReadListInd(fname,fFirstSensor,fLastSensor);
  fSensors->BypassStreamer(kFALSE);
  TTimeStamp defTime(2000,1,1,0,0,0);
  fStartTime = defTime;
  fEndTime = defTime;

}


//_____________________________________________________________________________
AliTPCSensorTempArray::AliTPCSensorTempArray(const AliTPCSensorTempArray &c):
  AliDCSSensorArray(c)
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
  fSensors->Delete();
  delete fSensors;

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
  AliDCSSensorArray::SetGraph(map,"tpc_temp:PT_%d.Temperature");
}  
//_____________________________________________________________________________
void AliTPCSensorTempArray::MakeSplineFit(TMap *map) 
{
  // 
  // Make spline fits from DCS maps 
  //
  AliDCSSensorArray::MakeSplineFit(map,"tpc_temp:PT_%d.Temperature");
}  


//_____________________________________________________________________________
TMap* AliTPCSensorTempArray::ExtractDCS(TMap *dcsMap) 
{
 //
 // Extract temperature graphs from DCS maps
 //

 TMap *values = AliDCSSensorArray::ExtractDCS(dcsMap,"tpc_temp:PT_%d.Temperature");
 return values;
}

