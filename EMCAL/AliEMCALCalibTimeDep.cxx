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

/* $Id: AliEMCALCalibTimeDep.cxx $ */

//_________________________________________________________________________
///*-- Author: 
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for EMCAL time-dep calibration                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <TGraphSmooth.h>
#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliEMCALSensorTempArray.h"
#include "AliEMCALCalibTimeDep.h"

/* first a bunch of constants.. */
const double fkSecToHour = 1.0/3600.0; // conversion factor from seconds to hours

// some global variables for APD handling; assume the same for all, at least to start with
const double fkTempSlope = 0.017; // 1.7% per deg. C, seems about right for all APDs, from studies at Catania, and by Rachid. 
// Note that this is only valid at gain M~30; i.e. voltages at V30. 
// At V50, it's more like 2.25%/C, and at V20 about 1.35%/C (also info from Catania)
const double fkNormTemp = 20; // let's say that 20 degrees C is normal as default

const double fkErrorCode = -999; // to indicate that something went wrong

using namespace std;

ClassImp(AliEMCALCalibTimeDep)

//________________________________________________________________
AliEMCALCalibTimeDep::AliEMCALCalibTimeDep() :
  fRun(0),
  fStartTime(0),
  fEndTime(0),
  fMinTemp(0),
  fMaxTemp(0),
  fMinTime(0),
  fMaxTime(0),
  fRefTemp(fkNormTemp), 
  fTempArray(NULL)
{
  // Constructor
}

//________________________________________________________________
AliEMCALCalibTimeDep::AliEMCALCalibTimeDep(const AliEMCALCalibTimeDep& calibt) :
  TObject(calibt),
  fRun(calibt.GetRunNumber()),
  fStartTime(calibt.GetStartTime()),
  fEndTime(calibt.GetEndTime()),
  fMinTemp(calibt.GetMinTemp()),
  fMaxTemp(calibt.GetMaxTemp()),
  fMinTime(calibt.GetMinTime()),
  fMaxTime(calibt.GetMaxTime()),
  fRefTemp(calibt.GetRefTemp()),
  fTempArray(calibt.GetTempArray())
{
  // copy constructor
}


//________________________________________________________________
AliEMCALCalibTimeDep &AliEMCALCalibTimeDep::operator =(const AliEMCALCalibTimeDep& calibt)
{
  // assignment operator; use copy ctor
  if (&calibt == this) return *this;

  new (this) AliEMCALCalibTimeDep(calibt);
  return *this;
}

//________________________________________________________________
AliEMCALCalibTimeDep::~AliEMCALCalibTimeDep()
{
  // Destructor
}

//________________________________________________________________
void  AliEMCALCalibTimeDep::Reset() 
{
  // clear variables to default
  fRun = 0;
  fStartTime = 0;
  fEndTime = 0;
  fMinTemp = 0;
  fMaxTemp = 0;
  fMinTime = 0;
  fMaxTime = 0;
  fRefTemp = fkNormTemp;
  fTempArray = NULL;
  return;
}

//________________________________________________________________
void  AliEMCALCalibTimeDep::PrintInfo() const
{
  // print some info
  cout << endl << " AliEMCALCalibTimeDep::Print() " << endl;
  // basic variables, all 'publicly available' also
  cout << " VARIABLE DUMP: " << endl
       << " GetStartTime() " << GetStartTime() << endl
       << " GetEndTime() " << GetEndTime() << endl
       << " GetMinTemp() " << GetMinTemp() << endl
       << " GetMaxTemp() " << GetMaxTemp() << endl
       << " GetRefTemp() " << GetRefTemp() << endl;
  // run ranges
  cout << " RUN INFO: " << endl
       << " length (in hours) " << GetLengthOfRunInHours() << endl
       << " range of temperature measurements (in hours) " << GetRangeOfTempMeasureInHours()
       << " (in deg. C) " << GetRangeOfTempMeasureInDegrees()
       << endl;
  // range in correction values
  double corrAtMinTemp = GetCorrection( fMinTemp );
  double corrAtMaxTemp = GetCorrection( fMaxTemp );
  double corrMaxMinDiff = 100*(corrAtMinTemp - corrAtMaxTemp);
  cout << " CORRECTION INFO : " << endl
       << " corrAtMinTemp " << corrAtMinTemp << endl
       << " corrAtMaxTemp " << corrAtMaxTemp << endl
       << " corrMaxMinDiff (~%) [=(corrAtMin - corrAtMax)*100] " 
       << corrMaxMinDiff << endl;

  return;
}
//________________________________________________________________ 
double AliEMCALCalibTimeDep::GetLengthOfRunInHours() const
{
  return (fEndTime - fStartTime)*fkSecToHour;
}
//________________________________________________________________ 
double AliEMCALCalibTimeDep::GetRangeOfTempMeasureInHours() const
{
  return (fMaxTime - fMinTime)*fkSecToHour;
}
//________________________________________________________________ 
double AliEMCALCalibTimeDep::GetRangeOfTempMeasureInDegrees() const
{
  return (fMaxTemp - fMinTemp);
}

//________________________________________________________________
void AliEMCALCalibTimeDep::Initialize(Int_t run, 
				      UInt_t startTime, UInt_t endTime)
{
  Reset(); // start fresh

  fRun = run;
  fStartTime = startTime;
  fEndTime = endTime;
  
  // collect the needed information
  GetTemperatureInfo(); // temperature readings during the run

  return;
}

//________________________________________________________________
double AliEMCALCalibTimeDep::GetTemperature(UInt_t timeStamp) const
{// return estimate for all SuperModules and sensors, that had data 

  // first convert from seconds to hours..
  double timeHour = (timeStamp - fStartTime) * fkSecToHour;

  double average = 0;
  int n = 0;

  for (int i=0; i<fTempArray->NumSensors(); i++) {
    
    AliEMCALSensorTemp *st = fTempArray->GetSensor(i);

    // check if we had valid data for the time that is being asked for
    if ( timeStamp>=st->GetStartTime() && timeStamp<=st->GetEndTime() ) {
      AliSplineFit *f = st->GetFit();
      if (f) { // ok, looks like we have valid data/info
	// let's check what the expected value at the time appears to be
	double val = f->Eval(timeHour);
	average += val;
	n++;
      }
    } // time
  } // loop over fTempArray
  
  if (n>0) { // some valid data was found
    average /= n;
    return average;
  }
  else { // no good data
    return fkErrorCode;
  }

}

//________________________________________________________________
double AliEMCALCalibTimeDep::GetTemperatureSM(int imod, UInt_t timeStamp) const
{// return estimate for this one SuperModule, if it had data 

  // first convert from seconds to hours..
  double timeHour = (timeStamp - fStartTime) * fkSecToHour;

  double average = 0;
  int n = 0;

  for (int i=0; i<fTempArray->NumSensors(); i++) {
    
    AliEMCALSensorTemp *st = fTempArray->GetSensor(i);
    int module = st->GetSector()*2 + st->GetSide();
    if ( module == imod ) { // right module
      // check if we had valid data for the time that is being asked for
      if ( timeStamp>=st->GetStartTime() && timeStamp<=st->GetEndTime() ) {
	AliSplineFit *f = st->GetFit();
	if (f) { // ok, looks like we have valid data/info
	  // let's check what the expected value at the time appears to be
	  double val = f->Eval(timeHour);
	  cout << " i " << i << " val " << val << endl;
	  average += val;
	  n++;
	}
      } // time
    }
    
  } // loop over fTempArray
  
  if (n>0) { // some valid data was found
    average /= n;
    return average;
  }
  else { // no good data
    return fkErrorCode;
  }

}

//________________________________________________________________
double AliEMCALCalibTimeDep::GetTemperatureSMSensor(int imod, int isens, UInt_t timeStamp) const
{// return estimate for this one SuperModule and sensor, if it had data 

  // first convert from seconds to hours..
  double timeHour = (timeStamp - fStartTime) * fkSecToHour;

  for (int i=0; i<fTempArray->NumSensors(); i++) {
    
    AliEMCALSensorTemp *st = fTempArray->GetSensor(i);
    int module = st->GetSector()*2 + st->GetSide();
    if ( module == imod && st->GetNum()==isens ) { // right module, and sensor
      // check if we had valid data for the time that is being asked for
      if ( timeStamp>=st->GetStartTime() && timeStamp<=st->GetEndTime() ) {
	AliSplineFit *f = st->GetFit();
	if (f) { // ok, looks like we have valid data/info
	  // let's check what the expected value at the time appears to be
	  double val = f->Eval(timeHour);

	  return val; // no point to move further in for loop, we have found the sensor we were looking for
	}
      } // time
    }
    
  } // loop over fTempArray
  
  // if we made it all here, it means that we didn't find the sensor we were looking for
  // i.e. no good data
  return fkErrorCode;

}

//________________________________________________________________
double AliEMCALCalibTimeDep::GetCorrection(double temperature) const
{ // gain decreases with increasing temperature
  double gain = (1.0 - (temperature-fRefTemp)*fkTempSlope);  
  // i.e. multiplicative correction to ADC increases with increasing temperature
  return 1.0/gain;
}

/* Next comes the method that does the work in picking up all the needed info..*/
//________________________________________________________________
void AliEMCALCalibTimeDep::GetTemperatureInfo() 
{
  // pick up Preprocessor output, based on fRun (most recent version)
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("EMCAL/Calib/Temperature", fRun);
  if (entry) {
    fTempArray = (AliEMCALSensorTempArray *) entry->GetObject();
  }

  if (fTempArray) { 
    AliInfo( Form("NumSensors %d - IdDCS: first %d last %d",
		  fTempArray->NumSensors(),
		  fTempArray->GetFirstIdDCS(), fTempArray->GetLastIdDCS() ) );
  }
  else {
    AliWarning( Form("AliEMCALSensorTempArray not found!") );
  }
  
  return;
}
