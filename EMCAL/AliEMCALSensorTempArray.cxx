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
//  EMCAL calibration class for saved temperature sensor parameters          //
//  Authors: David Silvermyr, copied from TPC (Ivanov, Helstrup)             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliEMCALSensorTempArray.h"
#include "TLinearFitter.h"
#include "TVectorD.h"
#include "AliLog.h"

ClassImp(AliEMCALSensorTempArray)

//_____________________________________________________________________________
AliEMCALSensorTempArray::AliEMCALSensorTempArray():AliDCSSensorArray()
{
  //
  // AliEMCALSensorTempArray default constructor
  //
}
//_____________________________________________________________________________
AliEMCALSensorTempArray::AliEMCALSensorTempArray(Int_t run) : AliDCSSensorArray() 
{
  //
  // Read configuration from OCDB
  //
  
  AliCDBEntry *entry =
    AliCDBManager::Instance()->Get("EMCAL/Config/Temperature",run); 
  if(entry){
    TTree *tree = (TTree*) entry->GetObject();
    fSensors = AliEMCALSensorTemp::ReadTree(tree);
    fSensors->BypassStreamer(kFALSE);
  }
  else AliFatal("CDB entry null!");
}
//_____________________________________________________________________________
AliEMCALSensorTempArray::AliEMCALSensorTempArray(UInt_t startTime, UInt_t endTime,
						 TTree* confTree, const TString& amandaString)
  :AliDCSSensorArray()
{
  //
  // AliEMCALSensorTempArray constructor for Shuttle preprocessor 
  //  (confTree read from OCDB)
  //
  fSensors = AliEMCALSensorTemp::ReadTree(confTree,amandaString);
  fSensors->BypassStreamer(kFALSE);
  fStartTime = TTimeStamp((time_t)startTime,0);
  fEndTime   = TTimeStamp((time_t)endTime,0);
}

//_____________________________________________________________________________
AliEMCALSensorTempArray::AliEMCALSensorTempArray(const char *fname,
                                          const TString& amandaString) :
                                                  AliDCSSensorArray()
{
  //
  // AliEMCALSensorTempArray constructor
  //
  fSensors = AliEMCALSensorTemp::ReadList(fname,amandaString);
  fSensors->BypassStreamer(kFALSE);
}

//_____________________________________________________________________________
AliEMCALSensorTempArray::AliEMCALSensorTempArray(const AliEMCALSensorTempArray &c):
  AliDCSSensorArray(c)
{
  //
  // AliEMCALSensorTempArray copy constructor
  //
}

//_____________________________________________________________________________
AliEMCALSensorTempArray::~AliEMCALSensorTempArray()
{
  //
  // AliEMCALSensorTempArray destructor
  //
}

//_____________________________________________________________________________
AliEMCALSensorTempArray &AliEMCALSensorTempArray::operator=(const AliEMCALSensorTempArray &c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    fSensors->Delete();
    new (this) AliEMCALSensorTempArray(c);
    fSensors = (TClonesArray*)c.fSensors->Clone();
  }
  return *this;
}

//_____________________________________________________________________________
void AliEMCALSensorTempArray::ReadSensors(const char *dbEntry)
{
  //
  // Read list of temperature sensors from text file
  //
  AliCDBEntry *entry = AliCDBManager::Instance()->Get(dbEntry);
  if(entry){
    TTree *tree = (TTree*) entry->GetObject();
    fSensors = AliEMCALSensorTemp::ReadTree(tree);
  }
  else AliFatal("NULL CDB entry!");
}  

//_____________________________________________________________________________
AliEMCALSensorTemp* AliEMCALSensorTempArray::GetSensor(Int_t side, Int_t sector, Int_t num) 
{
  //
  //  Return sensor information for sensor specified by side, sector and num
  //
  Int_t nsensors = fSensors->GetEntries();
  for (Int_t isensor=0; isensor<nsensors; isensor++) {
    AliEMCALSensorTemp *entry = (AliEMCALSensorTemp*)fSensors->At(isensor);
    if (entry->GetSide() == side &&
	entry->GetSector() == sector &&
	entry->GetNum() == num ) return entry;
  }
  return 0;
}

//_____________________________________________________________________________

AliEMCALSensorTemp* AliEMCALSensorTempArray::GetSensor(Int_t IdDCS){
  return dynamic_cast<AliEMCALSensorTemp*>(AliDCSSensorArray::GetSensor(IdDCS));
}

//_____________________________________________________________________________
AliEMCALSensorTemp* AliEMCALSensorTempArray::GetSensor(Double_t x, Double_t y, Double_t z){
  return dynamic_cast<AliEMCALSensorTemp*>(AliDCSSensorArray::GetSensor(x,y,z));
}

//_____________________________________________________________________________
Double_t AliEMCALSensorTempArray::GetTempGradientY(UInt_t timeSec, Int_t side){
  //
  // Extract Linear Vertical Temperature Gradient [K/cm] within the EMCAL on
  // Shaft Side(A): 0
  // Muon  Side(C): 1
  // Values based on TemperatureSensors within the EMCAL 
  //
  // FIXME: Also return residual-distribution, covariance Matrix
  //        or simply chi2 for validity check?
  //
  
  TLinearFitter fitter(3,"x0++x1++x2");
  TVectorD param(3);
  Int_t i = 0;
  
  Int_t nsensors = fSensors->GetEntries();
  for (Int_t isensor=0; isensor<nsensors; isensor++) { // loop over all sensors
    AliEMCALSensorTemp *entry = (AliEMCALSensorTemp*)fSensors->At(isensor);
    
    if (entry->GetSide()==side) { // only the selected side
      Double_t x[3];
      x[0]=1;
      x[1]=entry->GetX();
      x[2]=entry->GetY();
      Double_t y = entry->GetValue(timeSec); // get temperature value
      fitter.AddPoint(x,y,1); // add values to LinearFitter
      i++;
    }
    
  }
  fitter.Eval();
  fitter.GetParameters(param);
  
  return param[2]; // return vertical (Y) tempGradient
  
 }
