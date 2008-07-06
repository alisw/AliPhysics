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
const double fkTempSlope = 0.017; // 1.7% per deg. C, seems about right for all APDs, from studies at Catania, and by Rachid
const double fkNormTemp = 20; // let's say that 20 degrees C is normal as default

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
void  AliEMCALCalibTimeDep::Print() const
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
double AliEMCALCalibTimeDep::GetTemperature(int secSinceRunStart) const
{// return estimate for all SuperModules and sensors, that had data 

  // first convert from seconds to hours..
  double timeHour = secSinceRunStart * fkSecToHour;

  //  return fTempSpline->Eval(timeHour);
  return timeHour; // DStmp - FIXME - just return time for now
}

//________________________________________________________________
double AliEMCALCalibTimeDep::GetTemperatureSM(int imod, int secSinceRunStart) const
{// return estimate for this one SuperModule, if it had data 

  // first convert from seconds to hours..
  double timeHour = secSinceRunStart * fkSecToHour;

  //  return fTempSpline->Eval(timeHour);
  return timeHour; // DStmp - FIXME - just return time for now
}

//________________________________________________________________
double AliEMCALCalibTimeDep::GetTemperatureSMSensor(int imod, int isens, int secSinceRunStart) const
{// return estimate for this one SuperModule and sensor, if it had data 
  // first convert from seconds to hours..
  double timeHour = secSinceRunStart * fkSecToHour;
  
  //  return fTempSpline->Eval(timeHour);
  return timeHour; // DStmp - FIXME - just return time for now
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

  if (fTempArray) { // DStmp - FIXME - not implemented what should be done here yet
    AliInfo( Form("NumSensors %d - IdDCS: first %d last %d",
		  fTempArray->NumSensors(),
		  fTempArray->GetFirstIdDCS(), fTempArray->GetLastIdDCS() ) );

    /* // below is just examples on what could be done to access the data..
 
  AliEMCALSensorTemp *o = arr->GetSensor(0);
  o->Print();
  cout << " side " << o->GetSide()
       << " sector " << o->GetSector()
       << " num " << o->GetNum()
       << " startTime " << o->GetStartTime()
       << " endTime " << o->GetEndTime()
       << endl;

  AliSplineFit *f = o->GetFit();
  int np = 0;

  if (f) {
    f->SplineFit(0);

    np = f->GetKnots();
    cout << " np " << np << endl;
    Double_t *x = f->GetX();
    Double_t *y0 = f->GetY0();
    Double_t *y1 = f->GetY1();
    for (int i=0; i<np; i++) {
      cout << " i " << i
           << " x " << x[i]
           << " y0 " << y0[i]
           << " y1 " << y1[i]
           << endl;
      //      g->SetPoint(i, x[i], y0[i]);
    }
    double start = x[0];
    double stop = x[np-1];

    TGraph *g = f->MakeGraph(start, stop, np);
    f->MakeSmooth(g, ratio, "normal");
    np = f->GetKnots();
    cout << " np " << np << endl;
    TGraph *gSmooth = f->MakeGraph(start, stop, np);

  // compare the raw and smoothed versions and select one..
    // plot in 2 canvases/pads - CONTINUE HERE

    g->Draw("ALP");
    gSmooth->Draw("ALP");
    f->Draw();
    }

     */

  }
  
  return;
}
