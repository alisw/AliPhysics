#ifndef ALIMUONTRACKERHV_H
#define ALIMUONTRACKERHV_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

#include <vector>

class TMultiGraph;
class TMap;
class TGraph;

class AliMUONTrackerHV : public TObject
{
public:
  
  AliMUONTrackerHV(const char* runlist, const char* ocdbPath="raw://");
  AliMUONTrackerHV(Int_t runNumber, const char* ocdbPath="raw://");
  virtual ~AliMUONTrackerHV();

  void SetOCDB(const char* ocdbPath="raw://") { fOCDBPath = ocdbPath; }
  void SetRunList(Int_t runNumber);
  void SetRunList(const char* runlist);

  void HVoff(const char* logfile="lhc11de.log", const char* outputBaseName="hvoff");

  void Plot(const char* dcsname=0x0, Bool_t withPatch=kFALSE);
  
  void ReportTrips();
  
  void Scan(Int_t verbose=0);
  
private:

  void ReadIntegers(const char* filename, std::vector<int>& integers);

  TGraph* ShowValues(TMap* m, const char* name);
  
  Int_t CheckMap(TMap* hvMap, Int_t runNumber, Bool_t verbose);
  
  void TimeAxis(TMultiGraph* g);
  
  TMultiGraph* ShowHV(TMap* m, const char* dcsname);
  
private:
  std::vector<int> fRunList; // input run list
  TString fOCDBPath; // ocdb path (raw:// by default)
  
  ClassDef(AliMUONTrackerHV,1) // Utility class to inspect MUON Tracker HV values
};

#endif
