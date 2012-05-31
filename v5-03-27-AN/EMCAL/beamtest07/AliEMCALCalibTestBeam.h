#ifndef ALIEMCALCALIBTESTBEAM_H
#define ALIEMCALCALIBTESTBEAM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliEMCALCalibTestBeam.h $ */

////////////////////////////////////////////////
//  class for EMCAL testbeam calibration                 //
////////////////////////////////////////////////

#include "TGraph.h"
#include "TSpline.h"

class AliEMCALCalibTestBeam {

 public:
  AliEMCALCalibTestBeam(const int runNumber=186);
  AliEMCALCalibTestBeam(const AliEMCALCalibTestBeam &calibtb);
  AliEMCALCalibTestBeam& operator= (const AliEMCALCalibTestBeam &calibtb);
  virtual ~AliEMCALCalibTestBeam();
  virtual void Print() const; 

  // simple getters
  TGraph * GetTemperatureGraph() const { return fTempGraph; } // 
  TSpline3 * GetTemperatureSpline() const { return fTempSpline; } //
  int GetTimeStart() const { return fTimeStart; } // 
  int GetTimeStop() const { return fTimeStop; } // 
  double GetLengthOfRunInHours() const; // 
  int GetNEvents() const { return fNEvents; } // 

  int GetNTempVal() const { return fTempGraph->GetN(); } // 
  double GetMinTemp() const { return fMinTemp; } // 
  double GetMaxTemp() const { return fMaxTemp; } // 
  int GetMinTime() const { return fMinTime; } // 
  int GetMaxTime() const { return fMaxTime; } // 
  double GetRangeOfTempMeasureInHours() const; // 
  double GetRangeOfTempMeasureInDegrees() const; // 

  // basic calibration info
  double GetCorrection(int secSinceRunStart) const; // 
  double GetTemperature(int secSinceStart) const; //
  double GetCorrection(double temperature) const; // overloaded; maybe a bit dangerous - be careful with the type of arguments you supply..

 private:

  void Reset();
  void ResetRunInfo(); 
  void ResetTempInfo(); 
  void Init(const int runno);
  void GetRunTime(const int runno, 
		  const char *filename="/afs/cern.ch/user/d/dsilverm/www/testbeam07/calib/daqLogbook.root");
  void GetTemperatureInfo(const char* filename="/afs/cern.ch/user/d/dsilverm/www/testbeam07/calib/temperature-merged-alldays.root");
  void GetEMCALLogbookInfo(const int runno, 		
			   const char *filename="/afs/cern.ch/user/d/dsilverm/www/testbeam07/calib/EMCAL_Logbook_SPS_and_PS.csv");
  //
  TGraph *fTempGraph; // graph of temperature values
  TSpline3 *fTempSpline; // spline of  temperature values
  //
  int fTimeStart; // start time
  int fTimeStop; // stop time
  int fNEvents; // number of events
  double fMinTemp; // minimum time
  double fMaxTemp; // maximum time
  int fMinTime; // minimum time
  int fMaxTime; // maximum time
  //
  ClassDef(AliEMCALCalibTestBeam,1)    // EMCAL Calibration data
};

#endif
