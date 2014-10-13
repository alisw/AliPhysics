#ifndef ALITPCCALIBSUMMARY_H
#define ALITPCCALIBSUMMARY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCcalibSummary.h 29498 2008-10-25 12:18:56Z marian $ */

////////////////////////////////////////////////
//  Class to make a TPC calibration summary information                          //
////////////////////////////////////////////////

 
#include "TNamed.h"
class AliTPCcalibDB;
class AliTPCcalibDButil;
class TTreeSRedirector;

class AliTPCcalibSummary : public TNamed {

public:
  AliTPCcalibSummary();
  ~AliTPCcalibSummary();
  void Process(const char * runList, Int_t first, Int_t last);
  void ProcessRun(Int_t irun, Int_t startTime=0, Int_t endTime=0);
  //
  void ProcessDrift(Int_t run, Int_t timeStamp);
  void ProcessDriftCE(Int_t run, Int_t timeStamp);
  void ProcessDriftAll(Int_t run, Int_t timeStamp);
  void ProcessKryptonTime(Int_t run, Int_t timeStamp);
  void ProcessCTP(Int_t run, Int_t timeStamp);
  void ProcessAlign(Int_t run, Int_t timeStamp);
  void ProcessGain(Int_t run, Int_t timeStamp);
  void ProcessCurrent(Int_t irun,Int_t itime);

  void ProcessDriftCERef();
  void ProcessPulserRef();
protected:
  AliTPCcalibDB     *fCalibDB;      //! pointer to the TPC calib manager
  AliTPCcalibDButil *fDButil;       //! pointer to the TPC calib db utils
  TTreeSRedirector * fPcstream;     //! streamer - to store output info
private:
  AliTPCcalibSummary(const AliTPCcalibSummary&);
  AliTPCcalibSummary &operator=(const AliTPCcalibSummary&);
  ClassDef(AliTPCcalibSummary,0)  // 
};

#endif

