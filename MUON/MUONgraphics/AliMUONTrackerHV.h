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
class AliMpDCSNamer;

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

  void Plot(const char* dcsname=0x0, Bool_t withPatch=kFALSE, Bool_t plotIntermediate=kFALSE);

  void Print(Option_t* dcsname="") const;

  void ReportTrips(Bool_t includeLowOnes=kFALSE);
  
  void Scan(Int_t verbose=0);
  
  TGraph* Combine(TObjArray& graphs);

  TMultiGraph* CombineMulti(TObjArray& graphs);

  Int_t Compare(const TMap& hv1, const TMap& hv2, Bool_t verbose=kFALSE) const;
  
private:

  void ReadIntegers(const char* filename, std::vector<int>& integers);

  TGraph* GraphValues(TMap* m, const char* name);
  
  Int_t CheckMap(TMap* hvMap, Int_t runNumber, Bool_t verbose);
  
  void TimeAxis(TMultiGraph* g);
  
  TMultiGraph* GraphHV(TMap* m, const char* dcsname);
  
  AliMpDCSNamer* DCSNamer() const;

  AliMUONTrackerHV(const AliMUONTrackerHV& rhs); // not implemented on purpose
  AliMUONTrackerHV& operator=(const AliMUONTrackerHV& rhs); // not implemented on purpose
  
private:
  std::vector<int> fRunList; // input run list
  TString fOCDBPath; // ocdb path (raw:// by default)
  mutable AliMpDCSNamer* fDCSNamer; // helper to name things
  
  ClassDef(AliMUONTrackerHV,2) // Utility class to inspect MUON Tracker HV values
};

#endif
