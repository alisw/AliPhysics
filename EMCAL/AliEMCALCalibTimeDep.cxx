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
// class for EMCAL time-dep calibration   
// - supposed to run in preprocessor
// we use input from the following sources:
// AliEMCALBiasAPD (bias values),  AliEMCALCalibMapAPD (APD calibration and location info),
// AliCaloCalibSignal (LED DA), AliEMCALSensorTempArray (ELMB DCS)
// AliEMCALCalibReference: LED amplitude and temperature info at reference time
//
// output/result is in AliEMCALCalibTimeDepCorrection 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <TGraphSmooth.h>
#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliEMCALSensorTempArray.h"
#include "AliCaloCalibSignal.h"
#include "AliEMCALBiasAPD.h"
#include "AliEMCALCalibMapAPD.h"
#include "AliEMCALCalibReference.h"
#include "AliEMCALCalibTimeDepCorrection.h" 
#include "AliEMCALCalibTimeDep.h"

/* first a bunch of constants.. */
const double kSecToHour = 1.0/3600.0; // conversion factor from seconds to hours

// some global variables for APD handling; values from Catania studies, best fit
// TempCoeff = p0+p1*M (M=gain), where p0 and and p1 are functions of the dark current
const double kTempCoeffP0Const = -0.903; // 
const double kTempCoeffP0Factor = -1.381e7; // 
const double kTempCoeffP1Const = -0.023; // 
const double kTempCoeffP1Factor = -4.966e5; //
 
const double kErrorCode = -999; // to indicate that something went wrong

using namespace std;

ClassImp(AliEMCALCalibTimeDep)

//________________________________________________________________
AliEMCALCalibTimeDep::AliEMCALCalibTimeDep() :
  fRun(0),
  fStartTime(0),
  fEndTime(0),
  fMinTemp(0),
  fMaxTemp(0),
  fMinTempVariation(0),
  fMaxTempVariation(0),
  fMinTime(0),
  fMaxTime(0),
  fTemperatureResolution(0.1), // 0.1 deg C is default
  fTimeBinsPerHour(2), // 2 30-min bins per hour is default
  fHighLowGainFactor(16), // factor ~16 between High gain and low gain
  fTempArray(NULL),
  fCalibSignal(NULL),
  fBiasAPD(NULL),
  fCalibMapAPD(NULL),
  fCalibReference(NULL),
  fCalibTimeDepCorrection(NULL),
  fVerbosity(0)
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
  fMinTempVariation(calibt.GetMinTempVariation()),
  fMaxTempVariation(calibt.GetMaxTempVariation()),
  fMinTime(calibt.GetMinTime()),
  fMaxTime(calibt.GetMaxTime()),
  fTemperatureResolution(calibt.GetTemperatureResolution()),
  fTimeBinsPerHour(calibt.GetTimeBinsPerHour()),
  fHighLowGainFactor(calibt.GetHighLowGainFactor()),
  fTempArray(calibt.GetTempArray()),
  fCalibSignal(calibt.GetCalibSignal()),
  fBiasAPD(calibt.GetBiasAPD()),
  fCalibMapAPD(calibt.GetCalibMapAPD()),
  fCalibReference(calibt.GetCalibReference()),
  fCalibTimeDepCorrection(calibt.GetCalibTimeDepCorrection()),
  fVerbosity(calibt.GetVerbosity())
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
  fMinTempVariation = 0;
  fMaxTempVariation = 0;
  fMinTime = 0;
  fMaxTime = 0;
  fTemperatureResolution = 0.1; // 0.1 deg C is default
  fTimeBinsPerHour = 2; // 2 30-min bins per hour is default
  fTempArray = NULL;
  fCalibSignal = NULL;
  fBiasAPD = NULL;
  fCalibMapAPD = NULL;
  fCalibReference = NULL;
  fCalibTimeDepCorrection = NULL;
  fVerbosity = 0;
  return;
}

//________________________________________________________________
void  AliEMCALCalibTimeDep::PrintInfo() const
{
  // print some info
  cout << endl << " AliEMCALCalibTimeDep::PrintInfo() " << endl;
  // basic variables, all 'publicly available' also
  cout << " VARIABLE DUMP: " << endl
       << " GetStartTime() " << GetStartTime() << endl
       << " GetEndTime() " << GetEndTime() << endl
       << " GetMinTime() " << GetMinTime() << endl
       << " GetMaxTime() " << GetMaxTime() << endl
       << " GetMinTemp() " << GetMinTemp() << endl
       << " GetMaxTemp() " << GetMaxTemp() << endl
       << " GetMinTempVariation() " << GetMinTempVariation() << endl
       << " GetMaxTempVariation() " << GetMaxTempVariation() << endl;
  // run ranges
  cout << " RUN INFO: " << endl
       << " runnumber " << GetRunNumber() << endl
       << " length (in hours) " << GetLengthOfRunInHours() << endl
       << " range of temperature measurements (in hours) " << GetRangeOfTempMeasureInHours()
       << " (in deg. C) " << GetRangeOfTempMeasureInDegrees()
       << endl;

  return;
}

//________________________________________________________________ 
Double_t AliEMCALCalibTimeDep::GetLengthOfRunInHours() const
{
  return (fEndTime - fStartTime)*kSecToHour;
}

//________________________________________________________________ 
Double_t AliEMCALCalibTimeDep::GetLengthOfRunInBins() const
{
  return (fEndTime - fStartTime)*kSecToHour*fTimeBinsPerHour;
}

//________________________________________________________________ 
Double_t AliEMCALCalibTimeDep::GetRangeOfTempMeasureInHours() const
{
  return (fMaxTime - fMinTime)*kSecToHour;
}

//________________________________________________________________ 
Double_t AliEMCALCalibTimeDep::GetRangeOfTempMeasureInDegrees() const
{
  return (fMaxTemp - fMinTemp);
}

//________________________________________________________________
void AliEMCALCalibTimeDep::Initialize(Int_t run, 
				      UInt_t startTime, UInt_t endTime)
{ // setup, and get temperature info
  Reset(); // start fresh

  fRun = run;
  fStartTime = startTime;
  fEndTime = endTime;
  
  // collect the needed information
  GetTemperatureInfo(); // temperature readings during the run
  ScanTemperatureInfo(); // see what the boundaries are (Min/Max Time/Temp)

  return;
}

//________________________________________________________________
Double_t AliEMCALCalibTimeDep::GetTemperature(UInt_t timeStamp) const
{// return estimate for all SuperModules and sensors, that had data 

  // first convert from seconds to hours..
  Double_t timeHour = (timeStamp - fStartTime) * kSecToHour;

  Double_t average = 0;
  int n = 0;

  for (int i=0; i<fTempArray->NumSensors(); i++) {
    
    AliEMCALSensorTemp *st = fTempArray->GetSensor(i);

    // check if we had valid data for the time that is being asked for
    if ( timeStamp>=st->GetStartTime() && timeStamp<=st->GetEndTime() ) {
      AliSplineFit *f = st->GetFit();
      if (f) { // ok, looks like we have valid data/info
	// let's check what the expected value at the time appears to be
	Double_t val = f->Eval(timeHour);
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
    return kErrorCode;
  }

}

//________________________________________________________________
Double_t AliEMCALCalibTimeDep::GetTemperatureSM(int imod, UInt_t timeStamp) const
{// return estimate for this one SuperModule, if it had data 

  // first convert from seconds to hours..
  Double_t timeHour = (timeStamp - fStartTime) * kSecToHour;

  Double_t average = 0;
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
	  Double_t val = f->Eval(timeHour);
	  if ( fVerbosity > 0 ) {
	    cout << " sensor i " << i << " val " << val << endl;
	  }
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
    return kErrorCode;
  }

}

//________________________________________________________________
Double_t AliEMCALCalibTimeDep::GetTemperatureSMSensor(int imod, int isens, UInt_t timeStamp) const
{// return estimate for this one SuperModule and sensor, if it had data 

  // first convert from seconds to hours..
  Double_t timeHour = (timeStamp - fStartTime) * kSecToHour;

  for (int i=0; i<fTempArray->NumSensors(); i++) {
    
    AliEMCALSensorTemp *st = fTempArray->GetSensor(i);
    int module = st->GetSector()*2 + st->GetSide();
    if ( module == imod && st->GetNum()==isens ) { // right module, and sensor
      // check if we had valid data for the time that is being asked for
      if ( timeStamp>=st->GetStartTime() && timeStamp<=st->GetEndTime() ) {
	AliSplineFit *f = st->GetFit();
	if (f) { // ok, looks like we have valid data/info
	  // let's check what the expected value at the time appears to be
	  Double_t val = f->Eval(timeHour);

	  return val; // no point to move further in for loop, we have found the sensor we were looking for
	}
      } // time
    }
    
  } // loop over fTempArray
  
  // if we made it all here, it means that we didn't find the sensor we were looking for
  // i.e. no good data
  return kErrorCode;

}

//________________________________________________________________
Int_t AliEMCALCalibTimeDep::CalcCorrection()
{ // OK, this is where the real action takes place - the heart of this class..
  /* The philosophy is as follows:
     0. Init corrections to 1.0 values, and see how many correction bins we need
     1. Check how large temperature variations we have through the run - do we really need all the correction bias (otherwise adjust to single bin)
     2. try to use temperature info + APD bias and calibration info, to estimate correction.
     For now (from Dec 2009), we do not use LED info, since we are not taking LED triggers during the run.
   */

  // 0: Init
  // how many SuperModules do we have?
  Int_t nSM = fCalibReference->GetNSuperModule();
  // how many time-bins should we have for this run?
  Int_t nBins = (Int_t) GetLengthOfRunInBins(); // round-down (from double to int)
  Int_t binSize = (Int_t) (3600/fTimeBinsPerHour); // in seconds

  // 1: get info on how much individual sensors might have changed during
  // the run (compare max-min for each sensor separately)
  if (fMaxTempVariation < fTemperatureResolution) {
    nBins = 1; // just one bin needed..
    binSize = fEndTime - fStartTime;
  }

  // set up a reasonable default (correction = 1.0)
  fCalibTimeDepCorrection = new AliEMCALCalibTimeDepCorrection(nSM);
  fCalibTimeDepCorrection->InitCorrection(nSM, nBins, 1.0);
  fCalibTimeDepCorrection->SetStartTime(fStartTime);
  fCalibTimeDepCorrection->SetNTimeBins(nBins);
  fCalibTimeDepCorrection->SetTimeBinSize(binSize);

  // 2: try with Temperature correction 
  Int_t nRemaining = CalcTemperatureCorrection(nSM, nBins, binSize);

  return nRemaining;
}


//________________________________________________________________
Double_t AliEMCALCalibTimeDep::GetTempCoeff(Double_t IDark, Double_t M) const
{ // estimate the Temperature Coefficient, based on the dark current (IDark)
  // and the gain (M), based on Catania parameterizations

  Double_t dP0 = kTempCoeffP0Const + kTempCoeffP0Factor * IDark;
  Double_t dP1 = kTempCoeffP1Const + kTempCoeffP1Factor * IDark;
 
  Double_t dTC = dP0 + dP1*M;

  return dTC;
}

/* Next come the methods that do the work in picking up all the needed info..*/
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

//________________________________________________________________
Int_t AliEMCALCalibTimeDep::ScanTemperatureInfo() 
{// assign max/min time and temperature values

  fMinTemp = 999; // init to some large value (999 deg C)
  fMaxTemp = 0;
  fMinTempVariation = 999; // init to some large value (999 deg C)
  fMaxTempVariation = 0;
  fMinTime = 2147483647; // init to a large value in the far future (0x7fffffff), year 2038 times..
  fMaxTime = 0;

  Int_t n = 0; // number of valid readings

  for (int i=0; i<fTempArray->NumSensors(); i++) {
    
    AliEMCALSensorTemp *st = fTempArray->GetSensor(i);
    if ( st->GetStartTime() == 0 ) { // no valid data
      continue;
    }

    // check time ranges
    if (fMinTime > st->GetStartTime()) { fMinTime = st->GetStartTime(); }
    if (fMaxTime < st->GetEndTime()) { fMaxTime = st->GetEndTime(); }
 
    // check temperature ranges
    AliSplineFit *f = st->GetFit();

    if (f) { // ok, looks like we have valid data/info
      int np = f->GetKnots();
      Double_t *y0 = f->GetY0();
      // min and max values within the single sensor
      Double_t min = 999;
      Double_t max = 0;
      for (int ip=0; ip<np; ip++) { 
	if (min > y0[ip]) { min = y0[ip]; }
	if (max < y0[ip]) { max = y0[ip]; }
      }
      if (fMinTemp > min) { fMinTemp = min; }
      if (fMaxTemp < max) { fMaxTemp = max; }
      Double_t variation = max - min;
      if (fMinTempVariation > variation) { fMinTempVariation = variation; }
      if (fMaxTempVariation < variation) { fMaxTempVariation = variation; }

      n++;
    }
  } // loop over fTempArray
  
  if (n>0) { // some valid data was found
    return n;
  }
  else { // no good data
    return (Int_t) kErrorCode;
  }

}

//________________________________________________________________
void AliEMCALCalibTimeDep::GetCalibSignalInfo() 
{
  // pick up Preprocessor output, based on fRun (most recent version)
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("EMCAL/Calib/LED", fRun);
  if (entry) {
    fCalibSignal = (AliCaloCalibSignal *) entry->GetObject();
  }

  if (fCalibSignal) { 
    AliInfo( Form("CalibSignal: NEvents %d NAcceptedEvents %d Entries %lld AvgEntries LEDRefEntries %lld LEDRefEntries %lld, LEDRefAvgEntries %lld",
		  fCalibSignal->GetNEvents(), fCalibSignal->GetNAcceptedEvents(),
		  fCalibSignal->GetTreeAmpVsTime()->GetEntries(),
		  fCalibSignal->GetTreeAvgAmpVsTime()->GetEntries(),
		  fCalibSignal->GetTreeLEDAmpVsTime()->GetEntries(),
		  fCalibSignal->GetTreeLEDAvgAmpVsTime()->GetEntries() ) );		  
  }
  else {
    AliWarning( Form("AliCaloCalibSignal not found!") );
  }
  
  return;
}

//________________________________________________________________
void AliEMCALCalibTimeDep::GetBiasAPDInfo() 
{
  // pick up Preprocessor output, based on fRun (most recent version)
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("EMCAL/Calib/BiasAPD", fRun);
  if (entry) {
    fBiasAPD = (AliEMCALBiasAPD *) entry->GetObject();
  }

  if (fBiasAPD) { 
    AliInfo( Form("BiasAPD: NSuperModule %d ", fBiasAPD->GetNSuperModule() ) );
  }
  else {
    AliWarning( Form("AliEMCALBiasAPD not found!") );
  }
  
  return;
}

//________________________________________________________________
void AliEMCALCalibTimeDep::GetCalibMapAPDInfo() 
{
  // pick up Preprocessor output, based on fRun (most recent version)
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("EMCAL/Calib/MapAPD", fRun);
  // stored object should be a TTree; read the info
  if (entry) {
    fCalibMapAPD = (AliEMCALCalibMapAPD *) entry->GetObject();
  }

  if (fCalibMapAPD) { 
    AliInfo( Form("CalibMapAPD: NSuperModule %d ", fCalibMapAPD->GetNSuperModule() ) );
  }
  else {
    AliWarning( Form("AliEMCALCalibMapAPD not found!") );
  }
  
  return;
}

//________________________________________________________________
void AliEMCALCalibTimeDep::GetCalibReferenceInfo() 
{
  // pick up Preprocessor output, based on fRun (most recent version)
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("EMCAL/Calib/Reference", fRun);
  if (entry) {
    fCalibReference = (AliEMCALCalibReference *) entry->GetObject();
  }

  if (fCalibReference) { 
    AliInfo( Form("CalibReference: NSuperModule %d ", fCalibReference->GetNSuperModule() ) );
  }
  else {
    AliWarning( Form("AliEMCALCalibReference not found!") );
  }
  
  return;
}

//________________________________________________________________
Int_t AliEMCALCalibTimeDep::CalcLEDCorrection(Int_t nSM, Int_t nBins) 
{// Construct normalized ratios R(t)=LED(t)/LEDRef(t), for current time T and calibration time t0 
  // The correction factor we keep is c(T) = R(t0)/R(T)
  // T info from fCalibSignal, t0 info from fCalibReference

  // NOTE: for now we don't use the RMS info either from fCalibSignal or fCalibReference
  // but one could upgrade this in the future
  Int_t nRemaining = 0; // we count the towers for which we could not get valid data

  // sanity check; same SuperModule indices for corrections as for regular calibrations
  for (int i = 0; i < nSM; i++) {
    AliEMCALSuperModuleCalibReference * dataCalibReference = fCalibReference->GetSuperModuleCalibReferenceNum(i);
    AliEMCALSuperModuleCalibTimeDepCorrection * dataCalibTimeDepCorrection = fCalibTimeDepCorrection->GetSuperModuleCalibTimeDepCorrectionNum(i);

    int iSMRef = dataCalibReference->GetSuperModuleNum();
    int iSMCorr = dataCalibTimeDepCorrection->GetSuperModuleNum();
    if (iSMRef != iSMCorr) {
      AliWarning( Form("AliEMCALCalibTimeDep - SuperModule index mismatch: %d != %d", iSMRef, iSMCorr) );
      nRemaining = nSM * AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows * nBins;
      return nRemaining;
    }
  }
 
  int iSM = 0;
  int iCol = 0;
  int iRow = 0;
  int iStrip = 0;
  int iGain = 0;

  // The fCalibSignal info is stored in TTrees
  // Note that the time-bins for the TTree's may not exactly match with our correction time bins 
  int timeDiff = fCalibSignal->GetStartTime() - fStartTime; // in seconds
  // fCalibSignal time info in seconds: Hour/kSecToHour
  // corrected for startTime difference: Hour/kSecToHour + timeDiff
  // converted into a time-bin we use: (Hour + timeDiff*kSecToHour) * fTimeBinsPerHour

  // values for R(T), size of TArray = nBins
  // the [2] dimension below is for the low or high gain 
  TArrayF ampT[AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows][2]; 
  TArrayF nT[AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows][2]; 
  TArrayF ampLEDRefT[AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALLEDRefs][2]; 
  TArrayF nLEDRefT[AliEMCALGeoParams::fgkEMCALModules][AliEMCALGeoParams::fgkEMCALLEDRefs][2]; 

  // set up TArray's first
  for (iSM = 0; iSM < AliEMCALGeoParams::fgkEMCALModules; iSM++) {
    for (iCol = 0; iCol < AliEMCALGeoParams::fgkEMCALCols; iCol++) {
      for (iRow = 0; iRow < AliEMCALGeoParams::fgkEMCALRows; iRow++) {
	for (iGain = 0; iGain < 2; iGain++) {
	  // length of arrays
	  ampT[iSM][iCol][iRow][iGain].Set(nBins);
	  nT[iSM][iCol][iRow][iGain].Set(nBins);
	  // content of arrys
	  for (int j = 0; j < nBins; j++) {
	    ampT[iSM][iCol][iRow][iGain].AddAt(0, j);
	    nT[iSM][iCol][iRow][iGain].AddAt(0, j);
	  }
	}
      }
    }//iCol
    for (iStrip = 0; iStrip < AliEMCALGeoParams::fgkEMCALLEDRefs; iStrip++) {
      for (iGain = 0; iGain < 2; iGain++) {
	// length of arrays
	ampLEDRefT[iSM][iStrip][iGain].Set(nBins);
	nLEDRefT[iSM][iStrip][iGain].Set(nBins);
	// content of arrys
	for (int j = 0; j < nBins; j++) {
	  ampLEDRefT[iSM][iStrip][iGain].AddAt(0, j);
	  nLEDRefT[iSM][iStrip][iGain].AddAt(0, j);
	}
      }
    }//iStrip
  }

  // OK, now loop over the TTrees and fill the arrays needed for R(T)
  TTree *treeAvg = fCalibSignal->GetTreeAvgAmpVsTime();
  TTree *treeLEDRefAvg = fCalibSignal->GetTreeAvgAmpVsTime();

  int iChannelNum; // for regular towers
  int iRefNum; // for LED
  double dHour;
  double dAvgAmp;

  treeAvg->SetBranchAddress("fChannelNum", &iChannelNum);
  treeAvg->SetBranchAddress("fHour", &dHour);
  treeAvg->SetBranchAddress("fAvgAmp",&dAvgAmp);

  int iBin = 0;
  // counters for how many values were seen per SuperModule
  int nCount[AliEMCALGeoParams::fgkEMCALModules] = {0};
  int nCountLEDRef[AliEMCALGeoParams::fgkEMCALModules] = {0};

  for (int ient=0; ient<treeAvg->GetEntries(); ient++) {
    treeAvg->GetEntry(ient);
    // figure out where this info comes from
    fCalibSignal->DecodeChannelNum(iChannelNum, &iSM, &iCol, &iRow, &iGain);
    iBin = (int) ( (dHour + timeDiff*kSecToHour) * fTimeBinsPerHour  ); // CHECK!!!
    // add value in the arrays
    ampT[iSM][iCol][iRow][iGain].AddAt( ampT[iSM][iCol][iRow][iGain].At(iBin)+dAvgAmp, iBin );
    nT[iSM][iCol][iRow][iGain].AddAt( nT[iSM][iCol][iRow][iGain].At(iBin)+1, iBin );
    nCount[iSM]++;
  }

  treeLEDRefAvg->SetBranchAddress("fRefNum", &iRefNum);
  treeLEDRefAvg->SetBranchAddress("fHour", &dHour);
  treeLEDRefAvg->SetBranchAddress("fAvgAmp",&dAvgAmp);

  for (int ient=0; ient<treeLEDRefAvg->GetEntries(); ient++) {
    treeLEDRefAvg->GetEntry(ient);
    // figure out where this info comes from
    fCalibSignal->DecodeRefNum(iRefNum, &iSM, &iStrip, &iGain);
    iBin = (int) ( (dHour + timeDiff*kSecToHour) * fTimeBinsPerHour  ); // CHECK!!!
    // add value in the arrays
    ampLEDRefT[iSM][iStrip][iGain].AddAt( ampLEDRefT[iSM][iStrip][iGain].At(iBin)+dAvgAmp, iBin );
    nLEDRefT[iSM][iStrip][iGain].AddAt( nLEDRefT[iSM][iStrip][iGain].At(iBin)+1, iBin );
    nCountLEDRef[iSM]++;
  }

  // Normalize TArray values, and calculate average also
  Float_t norm = 0; // extra var, for readability

  for (iSM = 0; iSM < AliEMCALGeoParams::fgkEMCALModules; iSM++) {
    if (nCount[iSM]>0 && nCountLEDRef[iSM]>0) { // avoid SuperModules with no data..
      for (iCol = 0; iCol < AliEMCALGeoParams::fgkEMCALCols; iCol++) {
	//	iStrip = AliEMCALGeoParams::GetStripModule(iSM, iCol);
	iStrip = (iSM%2==0) ? iCol/2 : AliEMCALGeoParams::fgkEMCALLEDRefs - 1 - iCol/2; //TMP, FIXME
	for (iRow = 0; iRow < AliEMCALGeoParams::fgkEMCALRows; iRow++) {
	  for (iGain = 0; iGain < 2; iGain++) {
	    // content of arrys
	    for (int j = 0; j < nBins; j++) {
	      if (nT[iSM][iCol][iRow][iGain].At(j) > 0) { 
		norm = ampT[iSM][iCol][iRow][iGain].At(j) / nT[iSM][iCol][iRow][iGain].At(j); 
		ampT[iSM][iCol][iRow][iGain].AddAt(norm, j); // AddAt = SetAt
	      }	
	    }
	  }
	}
      }//iCol
      for (iStrip = 0; iStrip < AliEMCALGeoParams::fgkEMCALLEDRefs; iStrip++) {
	for (iGain = 0; iGain < 2; iGain++) {
	  for (int j = 0; j < nBins; j++) {
	    if (nLEDRefT[iSM][iStrip][iGain].At(j) > 0) {
	      norm = ampLEDRefT[iSM][iStrip][iGain].At(j) / nLEDRefT[iSM][iStrip][iGain].At(j); 
	      ampLEDRefT[iSM][iStrip][iGain].AddAt(norm, j); // AddAt = SetAt
	    }
	  }
	}
      }//iStrip
    }
  } // iSM


  // Calculate correction values, and store them
  // set kErrorCode values for those that could not be set

  Float_t ratiot0 = 0;
  Float_t ratioT = 0;
  Float_t correction = 0; // c(T) = R(t0)/R(T)
  Int_t refGain = 0; // typically use low gain for LED reference amplitude (high gain typically well beyond saturation)

  for (int i = 0; i < nSM; i++) {
    AliEMCALSuperModuleCalibReference * dataCalibReference = fCalibReference->GetSuperModuleCalibReferenceNum(i);
    AliEMCALSuperModuleCalibTimeDepCorrection * dataCalibTimeDepCorrection = fCalibTimeDepCorrection->GetSuperModuleCalibTimeDepCorrectionNum(i);
    iSM = dataCalibReference->GetSuperModuleNum();

    for (iCol = 0; iCol < AliEMCALGeoParams::fgkEMCALCols; iCol++) {
      //      iStrip = AliEMCALGeoParams::GetStripModule(iSM, iCol);
      iStrip = (iSM%2==0) ? iCol/2 : AliEMCALGeoParams::fgkEMCALLEDRefs - 1 - iCol/2; //TMP, FIXME
	refGain = dataCalibReference->GetLEDRefHighLow(iStrip); // LED reference gain value used for reference calibration	

      for (iRow = 0; iRow < AliEMCALGeoParams::fgkEMCALRows; iRow++) {

	// Calc. R(t0):
	AliEMCALCalibReferenceVal * refVal = dataCalibReference->GetAPDVal(iCol, iRow);
	iGain = refVal->GetHighLow(); // gain value used for reference calibration	
	// valid amplitude values should be larger than 0
	if (refVal->GetLEDAmp()>0 && dataCalibReference->GetLEDRefAmp(iStrip)>0) {
	  ratiot0 =  refVal->GetLEDAmp() / dataCalibReference->GetLEDRefAmp(iStrip);
	}
	else {
	  ratiot0 = kErrorCode;
	}

	// Calc. R(T)
	for (int j = 0; j < nBins; j++) {

	  // calculate R(T) also; first try with individual tower:
	  // same gain as for reference calibration is the default
	  if (ampT[iSM][iCol][iRow][iGain].At(j)>0 && ampLEDRefT[iSM][iStrip][refGain].At(j)>0) {
	    // looks like valid data with the right gain combination
	    ratioT = ampT[iSM][iCol][iRow][iGain].At(j) / ampLEDRefT[iSM][iStrip][refGain].At(j); 

	    // if data appears to be saturated, and we are in high gain, then try with low gain instead
	    int newGain = iGain;
	    int newRefGain = refGain;
	    if ( ampT[iSM][iCol][iRow][iGain].At(j)>AliEMCALGeoParams::fgkOverflowCut && iGain==1 ) {
	      newGain = 0;
	    }
	    if ( ampLEDRefT[iSM][iStrip][refGain].At(j)>AliEMCALGeoParams::fgkOverflowCut && refGain==1 ) { 
	      newRefGain = 0;
	    }

	    if (newGain!=iGain || newRefGain!=refGain) {
	      // compensate for using different gain than in the reference calibration
	      // we may need to have a custom H/L ratio value for each tower
	      // later, but for now just use a common value, as the rest of the code does.. 
	      ratioT = ampT[iSM][iCol][iRow][newGain].At(j) / ampLEDRefT[iSM][iStrip][newRefGain].At(j); 

	      if (newGain<iGain) {
		ratioT *= fHighLowGainFactor;
	      }
	      else if (newRefGain<refGain) {
		ratioT /= fHighLowGainFactor;
	      }
	    }
	  }
	  else {
	    ratioT = kErrorCode;
	  }

	  // Calc. correction factor
	  if (ratiot0>0 && ratioT>0) { 
	    correction = ratiot0/ratioT;
	  }
	  else {
	    correction = kErrorCode;
	    nRemaining++;
	  }

	  // Store the value
	  dataCalibTimeDepCorrection->GetCorrection(iCol,iRow)->AddAt(correction, j);
	  /* Check that
	  fTimeDepCorrection->SetCorrection(i, iCol, iRow, j, correction);
	  also works OK */
	} // nBins
      }
    }
  }

  nRemaining = CalcLEDCorrectionStripBasis(nSM, nBins);
  return nRemaining;
}

//________________________________________________________________
Int_t AliEMCALCalibTimeDep::CalcLEDCorrectionStripBasis(Int_t nSM, Int_t nBins) 
{ // use averages for each strip if no good values exist for some single tower 

  // go over fTimeDepCorrection info
  Int_t nRemaining = 0; // we count the towers for which we could not get valid data

  // for calculating StripAverage info
  int nValidTower = 0;
  Float_t stripAverage = 0;
  Float_t val = 0;

  int iSM = 0;
  int iCol = 0;
  int iRow = 0;
  int iStrip = 0;
  int firstCol = 0;
  int lastCol = 0;

  for (int i = 0; i < nSM; i++) {
    AliEMCALSuperModuleCalibTimeDepCorrection * dataCalibTimeDepCorrection = fCalibTimeDepCorrection->GetSuperModuleCalibTimeDepCorrectionNum(i);
    iSM = dataCalibTimeDepCorrection->GetSuperModuleNum();

    for (int j = 0; j < nBins; j++) {

      nValidTower = 0;
      stripAverage = 0;

      for (iStrip = 0; iStrip < AliEMCALGeoParams::fgkEMCALLEDRefs; iStrip++) {
	firstCol = iStrip*2;
	if ((iSM%2)==1) { // C side
	  firstCol = (AliEMCALGeoParams::fgkEMCALLEDRefs-1 - iStrip)*2;
	}
	lastCol = firstCol+1;

	for (iCol = firstCol; iCol <= lastCol; iCol++) {
	  for (iRow = 0; iRow < AliEMCALGeoParams::fgkEMCALRows; iRow++) {
	    val = dataCalibTimeDepCorrection->GetCorrection(iCol,iRow)->At(j);
	    if (val>0) { // valid value; error code is negative
	      stripAverage += val;
	      nValidTower++;
	    } 
	  }
	}

	// calc average over strip
	if (nValidTower>0) {
	  stripAverage /= nValidTower;
	  for (iCol = firstCol; iCol <= lastCol; iCol++) {
	    for (iRow = 0; iRow < AliEMCALGeoParams::fgkEMCALRows; iRow++) {
	      val = dataCalibTimeDepCorrection->GetCorrection(iCol,iRow)->At(j);
	      if (val<0) { // invalid value; error code is negative
		dataCalibTimeDepCorrection->GetCorrection(iCol,iRow)->AddAt(val, j);
	      }
	    }
	  }
	}
	else { // could not fill in unvalid entries
	  nRemaining += 2*AliEMCALGeoParams::fgkEMCALRows;
	}

      } // iStrip
    } // j, bins
  } // iSM

  return nRemaining;
}

//________________________________________________________________
Int_t AliEMCALCalibTimeDep::CalcTemperatureCorrection(Int_t nSM, Int_t nBins, Int_t binSize) 
{ // OK, so we didn't have valid LED data that allowed us to do the correction only 
  // with that info.
  // So, instead we'll rely on the temperature info and try to do the correction
  // based on that instead.
  // For this, we'll need the info from 3 classes (+temperature array), and output the values in a 4th class 
  Int_t nRemaining = 0;

  int iSM = 0;
  int iCol = 0;
  int iRow = 0;

  Double_t dTempCoeff[AliEMCALGeoParams::fgkEMCALCols][AliEMCALGeoParams::fgkEMCALRows];
  memset(dTempCoeff, 0, sizeof(dTempCoeff));
  Float_t gainM = 0; 
  Double_t correction = 0;
  Double_t secondsPerBin = (Double_t) binSize;

  for (int i = 0; i < nSM; i++) {
    AliEMCALSuperModuleCalibTimeDepCorrection * dataCalibTimeDepCorrection = fCalibTimeDepCorrection->GetSuperModuleCalibTimeDepCorrectionNum(i);
    iSM = dataCalibTimeDepCorrection->GetSuperModuleNum();

    AliEMCALSuperModuleCalibReference * dataCalibReference = fCalibReference->GetSuperModuleCalibReferenceNum(iSM);
    AliEMCALSuperModuleCalibMapAPD * dataCalibMapAPD = fCalibMapAPD->GetSuperModuleCalibMapAPDNum(iSM);
    AliEMCALSuperModuleBiasAPD * dataBiasAPD = fBiasAPD->GetSuperModuleBiasAPDNum(iSM);    

    // first calculate the M=Gain values, and TemperatureCoeff, for all towers in this SuperModule, from BiasAPD and CalibMapAPD info
    for (iCol = 0; iCol < AliEMCALGeoParams::fgkEMCALCols; iCol++) {
      for (iRow = 0; iRow < AliEMCALGeoParams::fgkEMCALRows; iRow++) {
	AliEMCALCalibMapAPDVal * mapAPD = dataCalibMapAPD->GetAPDVal(iCol, iRow);
	gainM =  fCalibMapAPD->GetGain(mapAPD->GetPar(0), mapAPD->GetPar(1), mapAPD->GetPar(2), 
				       dataBiasAPD->GetVoltage(iCol, iRow));
	dTempCoeff[iCol][iRow] = GetTempCoeff(mapAPD->GetDarkCurrent(), gainM);
	if (fVerbosity > 1) {
	  cout << " iSM " << iSM << " iCol " << iCol << " iRow " << iRow 
	       << " par0 " << mapAPD->GetPar(0) 
	       << " par1 " << mapAPD->GetPar(1) 
	       << " par2 " << mapAPD->GetPar(2)
	       << " bias " << dataBiasAPD->GetVoltage(iCol, iRow) 
	       << " gainM " << gainM << " dTempCoeff " << dTempCoeff[iCol][iRow] << endl; 
	}
      }
    }
    
    // figure out what the reference temperature is, from fCalibReference
    Double_t referenceTemperature = 0;
    int nVal = 0;
    for (int iSensor = 0; iSensor<AliEMCALGeoParams::fgkEMCALTempSensors ; iSensor++) {
      if (dataCalibReference->GetTemperature(iSensor)>0) { // hopefully OK value
	referenceTemperature += dataCalibReference->GetTemperature(iSensor);
	nVal++;
      }
    }
    
    if (nVal>0) {
      referenceTemperature /= nVal; // valid values exist, we can look into corrections
      
      for (int j = 0; j < nBins; j++) {
	// what is the timestamp in the middle of this bin? (0.5 is for middle of bin)
	UInt_t timeStamp = fStartTime + (UInt_t)((j+0.5)*secondsPerBin);
      // get the temperature at this time; use average over whole SM for now (TO BE CHECKED LATER - if we can do better with finer grained info)
	Double_t dSMTemperature = GetTemperatureSM(iSM, timeStamp); 

	Double_t temperatureDiff = dSMTemperature - referenceTemperature ; // new - old
	if (fVerbosity > 0) {
	  cout << " referenceTemperature " << referenceTemperature
	       << " dSMTemperature " << dSMTemperature
	       << " temperatureDiff " << temperatureDiff
	       << endl;
	}
	// if the new temperature is higher than the old/reference one, then the gain has gone down 
	// if the new temperature is lower than the old/reference one, then the gain has gone up
	// dTempCoeff is a negative number describing how many % (hence 0.01 factor below) the gain 
	// changes with a positive degree change.  
	// i.e. the product temperatureDiff * dTempCoeff increase when the gain goes up
	// The correction we want to keep is what we should multiply our ADC value with as a function
	// of time, i.e. the inverse of the gain change..
	if (fabs(temperatureDiff)>fTemperatureResolution) { 
	  // significant enough difference that we need to consider it

	  // loop over all towers; effect of temperature change will depend on gain for this tower
	  for (iCol = 0; iCol < AliEMCALGeoParams::fgkEMCALCols; iCol++) {
	    for (iRow = 0; iRow < AliEMCALGeoParams::fgkEMCALRows; iRow++) {

	      // the correction should be inverse of modification in gain: (see discussion above)
	      // modification in gain: 1.0 + (temperatureDiff * dTempCoeff[iCol][iRow])*0.01; 
	      // 1/(1+x) ~= 1 - x for small x, i.e. we arrive at:
	      correction = 1.0 - (temperatureDiff * dTempCoeff[iCol][iRow])*0.01; 
	      dataCalibTimeDepCorrection->GetCorrection(iCol,iRow)->AddAt(correction, j);
	      if (fVerbosity > 1) {
		cout << " iSM " << iSM
		     << " iCol  " << iCol 
		     << " iRow  " << iRow 
		     << " j  " << j 
		     << " correction  " << correction 
		     << endl;
	      }
	    }
	  }

	} // if noteworthy temperature change
	else { // just set correction values to 1.0
	  correction = 1;
	  for (iCol = 0; iCol < AliEMCALGeoParams::fgkEMCALCols; iCol++) {
	    for (iRow = 0; iRow < AliEMCALGeoParams::fgkEMCALRows; iRow++) {
	      dataCalibTimeDepCorrection->GetCorrection(iCol,iRow)->AddAt(correction, j);
	    }
	  }
	} // else
      } // j, Bins

    } // if reference temperature exist 
    else { // could not do the needed check.. signal that in the return code
      nRemaining += AliEMCALGeoParams::fgkEMCALCols * AliEMCALGeoParams::fgkEMCALRows * nBins;
    }
  } // iSM

  return nRemaining;
}

