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

/* $Id: AliEMCALCalibTestBeam.cxx $ */

//_________________________________________________________________________
///*-- Author: 
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for EMCAL testbeam calibration                                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TGraphSmooth.h>
#include <TFile.h>
#include <TTree.h>
#include "AliEMCALCalibTestBeam.h"

/* first a bunch of constants.. */
const double fkSmoothBandwidth = 5.0; // 5 seems to be a good value..; see $ROOTSYS/tutorials/graphs/motorcycle.C for reference
const double fkSecToHour = 1.0/3600.0; // conversion factor from seconds to hours

// some global variables for APD handling
const double fkTempSlope = 0.017; // 1.7% per deg. C, seems about right for all APDs, from studies at Catania, and by Rachid
const double fkNormTemp = 20; // let's say that 20 degrees is normal, i.e. we'll do the corrections relative to that reference temperature  

// some global constants and ROOT objects for temperature handling
const int fkSensor = 1; // let's just look at temperature data from one sensor
// let's also have some ranges for valid readings of the temperatures..
const int fkMinValidTemp = 0;
const int fkMaxValidTemp = 100;

using namespace std;

ClassImp(AliEMCALCalibTestBeam)

//________________________________________________________________
AliEMCALCalibTestBeam::AliEMCALCalibTestBeam(const int runNumber)
{
  // Constructor
  // set pointers even though they will be re-assigned in Init()/Reset() later.. 
  fTempGraph = NULL; 
  fTempSpline = NULL;
  Init(runNumber); // collect and setup info; including the spline
}

//________________________________________________________________
AliEMCALCalibTestBeam::AliEMCALCalibTestBeam(const AliEMCALCalibTestBeam& calibtb) :
  fTimeStart(calibtb.GetTimeStart()),
  fTimeStop(calibtb.GetTimeStop()),
  fMinTemp(calibtb.GetMinTemp()),
  fMaxTemp(calibtb.GetMaxTemp()),
  fMinTime(calibtb.GetMinTime()),
  fMaxTime(calibtb.GetMaxTime())
{
  // copy constructor
  fTempGraph = calibtb.GetTemperatureGraph();
  fTempSpline = calibtb.GetTemperatureSpline();
}


//________________________________________________________________
AliEMCALCalibTestBeam &AliEMCALCalibTestBeam::operator =(const AliEMCALCalibTestBeam& calibtb)
{
  // assignment operator; use copy ctor
  if (&calibtb == this) return *this;

  new (this) AliEMCALCalibTestBeam(calibtb);
  return *this;

}

//________________________________________________________________
AliEMCALCalibTestBeam::~AliEMCALCalibTestBeam()
{
  // Destructor
  if (fTempGraph) delete fTempGraph;
  if (fTempSpline) delete fTempSpline;
}

//________________________________________________________________
void  AliEMCALCalibTestBeam::Reset() 
{ // clear variables
  ResetRunInfo();
  ResetTempInfo();
  return;
}
//________________________________________________________________
void  AliEMCALCalibTestBeam::ResetRunInfo() 
{
  // clear variables
  fTimeStart = 0;
  fTimeStop = 0;
  fNEvents = 0;
}

//________________________________________________________________
void  AliEMCALCalibTestBeam::ResetTempInfo() 
{
  // clear variables
  fMinTemp = 0;
  fMaxTemp = 0;
  fMinTime = 0;
  fMaxTime = 0;

  if (fTempGraph) delete fTempGraph;
  fTempGraph = new TGraph();

  if (fTempSpline) delete fTempSpline;
  fTempSpline = NULL;

  return;
}

//________________________________________________________________
void  AliEMCALCalibTestBeam::Print() const
{
  // print some info
  cout << endl << " AliEMCALCalibTestBeam::Print() " << endl;
  // basic variables, all 'publicly available' also
  cout << " VARIABLE DUMP: " << endl
       << " GetTimeStart() " << GetTimeStart() << endl
       << " GetTimeStop() " << GetTimeStop() << endl
       << " GetNTempVal() " << GetNTempVal() << endl
       << " GetMinTemp() " << GetMinTemp() << endl
       << " GetMaxTemp() " << GetMaxTemp() << endl;
  // run ranges
  cout << " RUN INFO: " << endl
       << " NEvents " << GetNEvents() << endl
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
double AliEMCALCalibTestBeam::GetLengthOfRunInHours() const
{ // run length in hours
  return (fTimeStop - fTimeStart)*fkSecToHour;
}
//________________________________________________________________ 
double AliEMCALCalibTestBeam::GetRangeOfTempMeasureInHours() const
{ // interval with temperature measurements
  return (fMaxTime - fMinTime)*fkSecToHour;
}
//________________________________________________________________ 
double AliEMCALCalibTestBeam::GetRangeOfTempMeasureInDegrees() const
{ // range of temperatures
  return (fMaxTemp - fMinTemp);
}

//________________________________________________________________
void AliEMCALCalibTestBeam::Init(const int runno)
{ // init/clear info
  Reset(); // start fresh

  // collect the needed information
  GetRunTime(runno); // when was the run exactly..
  GetTemperatureInfo(); // temperature readings during the run
  GetEMCALLogbookInfo(runno); // logbook notes during the run

  // make a smooth graph and spline, if we can
  if (fTempGraph->GetN() > 0) { // need to have some points 
    TGraphSmooth smooth;
    TGraph * graphSmoothed = smooth.SmoothKern(fTempGraph,"normal", fkSmoothBandwidth);
    // use a cubic spline on the smoothed graph..
    fTempSpline = new TSpline3("fTempSpline",graphSmoothed);
  }
  return;
}

//________________________________________________________________
double AliEMCALCalibTestBeam::GetCorrection(int secSinceRunStart) const
{ // calculate correction
  double temperature = GetTemperature(secSinceRunStart);
  return GetCorrection(temperature);
}

//________________________________________________________________
double AliEMCALCalibTestBeam::GetTemperature(int secSinceRunStart) const
{
  // first convert from seconds to hours..
  double timeHour = secSinceRunStart * fkSecToHour;
  return fTempSpline->Eval(timeHour);
}

//________________________________________________________________
double AliEMCALCalibTestBeam::GetCorrection(double temperature) const
{ // gain decreases with increasing temperature
  double gain = (1.0 - (temperature-fkNormTemp)*fkTempSlope);  
  // i.e. multiplicative correction to ADC increases with increasing temperature
  return 1.0/gain;
}

/* Next comes the method that does the work in picking up all the needed info..*/
//________________________________________________________________
void AliEMCALCalibTestBeam::GetRunTime(const int runno, 
				       const char *filename)
{ // get run info
  TFile *fR = new TFile(filename); 
  if (!fR) {
    printf("file %s could not be found. Perhaps you don't have afs working?\n",filename);
    printf("Then you can try to get it with:\n wget http://cern.ch/dsilverm/testbeam07/calib/daqLogbook.root");
  }

  TTree *daqLogbook = (TTree*)gDirectory->Get("daqLogbook");
  
  //Declaration of leaves types
  Int_t           runNumber;
  Int_t           timeStart;
  Int_t           runDuration;
  Int_t           totalEvents;
  Float_t         totalData;
  Float_t         averageDataRate;
  Float_t         averageEventsPerSecond;
  
  // Set branch addresses.
  daqLogbook->SetBranchAddress("runNumber",&runNumber);
  daqLogbook->SetBranchAddress("timeStart",&timeStart);
  daqLogbook->SetBranchAddress("runDuration",&runDuration);
  daqLogbook->SetBranchAddress("totalEvents",&totalEvents);
  daqLogbook->SetBranchAddress("totalData",&totalData);
  daqLogbook->SetBranchAddress("averageDataRate",&averageDataRate);
  daqLogbook->SetBranchAddress("averageEventsPerSecond",&averageEventsPerSecond);
  
   Long64_t nentries = daqLogbook->GetEntries();
   Long64_t nbytes = 0;

   for (Long64_t i=0; i<nentries;i++) {
     nbytes += daqLogbook->GetEntry(i);
     if (runNumber == runno) {
       cout << " getRunTime: runNumber " << runNumber
	    << " timeStart " << timeStart
	    << " runDuration (seconds) " << runDuration
	    << " totalEvents " << totalEvents
	    << " totalData (kByte) " << totalData
	    << endl;

       // assign the start and stop unix-timestamps 
       fTimeStart = timeStart;
       fTimeStop = timeStart + runDuration;
       fNEvents = totalEvents;
       return; // then we are done..; should only pick up info for one run maximum..
     }
   }

  return;
}

//________________________________________________________________
void AliEMCALCalibTestBeam::GetTemperatureInfo(const char* filename)
{ // get temperature info
  fTempGraph->Clear();

  TFile *fT = new TFile(filename);
  if (!fT) {
    printf("file %s could not be found. Perhaps you don't have afs working?\n",filename);
    printf("Then you can try to get it with:\n wget http://cern.ch/dsilverm/testbeam07/calib/temperature-merged-alldays.root");
  }
  TTree *tree = (TTree*)gDirectory->Get("tree");
  
  //Declaration of leaves types
  Int_t           unixTime;
  Float_t         temperature;
  Int_t           SensorNumber;
  
  // Set branch addresses.
  tree->SetBranchAddress("unixTime",&unixTime);
  tree->SetBranchAddress("temperature",&temperature);
  tree->SetBranchAddress("SensorNumber",&SensorNumber);
  
  Long64_t nentries = tree->GetEntries();  
  Long64_t nbytes = 0;

  // init counter and range info
  int np = 0; // number of points
  fMinTemp = fkMaxValidTemp;
  fMaxTemp = fkMinValidTemp;
  fMinTime = 2000000000; // some timestamp far in the future..
  fMaxTime = 0;

  for (Long64_t i=0; i<nentries;i++) {
    nbytes += tree->GetEntry(i);
    if (unixTime>=fTimeStart && unixTime<=fTimeStop && SensorNumber==fkSensor) {
      if (temperature>=fkMinValidTemp && temperature<=fkMaxValidTemp) {
	/*
	  cout << " i " << i
	  << " unixTime " << unixTime
	  << " temperature " << temperature
	  << " SensorNumber " << SensorNumber
	  << endl;
	*/
	if (fMinTemp > temperature) fMinTemp = temperature;
	if (fMaxTemp < temperature) fMaxTemp = temperature;
	if (fMinTime > unixTime) fMinTime = unixTime;
	if (fMaxTime < unixTime) fMaxTime = unixTime;

	fTempGraph->SetPoint(np, (unixTime-fTimeStart)*fkSecToHour, temperature);
	np++;
      } // valid temperature reading
    } // valid time reading
  }

  // what if no values were found? Then we reset the max/min values etc. to the defaults
  if (np==0) {
    ResetTempInfo(); // partial Reset()
  }
 
  cout << " Number of temperature points found: np= " << np 
       << " for sensor " << fkSensor 
       << " : maxtemp " << fMaxTemp
       << " mintemp " << fMinTemp
       << endl;
  return;
}

//________________________________________________________________
void AliEMCALCalibTestBeam::GetEMCALLogbookInfo(const int runno, 		
						const char *filename)
{ // get logbook info
  Int_t Run,DAQEvt,CCUSB,Counter,XPos,YPos,Momentum;
  char TriggerAndComment[1000];

  ifstream fin;
  fin.open(filename);
  if(!fin.good()) {
    printf("file %s could not be found. Perhaps you don't have afs working?\n",filename);
    printf("Then you can try to get it with:\n wget http://cern.ch/dsilverm/testbeam07/calib/EMCAL_Logbook_SPS_and_PS.csv");
  }

  char line[1000];
  int nlines=0;
  while (!fin.getline(line,1000).eof()) {
    if(!fin.good()) break;

    sscanf(line,"%i,%i,%i,%i,%i,%i,%i,%s",
	   &Run,&DAQEvt,&CCUSB,&Counter,&XPos,&YPos,&Momentum,
	   TriggerAndComment);

    if (Run == runno) { // we found the info
      cout << " EMCAL Logbook info line: " << line << endl;
    }

    nlines++;
  }
  return;
}


