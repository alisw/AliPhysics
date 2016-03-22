#ifndef ALIMUONTRACKERVOLTAGES_H
#define ALIMUONTRACKERVOLTAGES_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

#include <vector>
#include <cfloat>

class TMultiGraph;
class TMap;
class TGraph;
class AliMpDCSNamer;

class AliMUONTrackerVoltages : public TObject
{
public:

  AliMUONTrackerVoltages(const char* runlist, const char* ocdbPath="raw://");
  AliMUONTrackerVoltages(Int_t runNumber, const char* ocdbPath="raw://");
  virtual ~AliMUONTrackerVoltages();

  void SetOCDB(const char* ocdbPath="raw://") { fOCDBPath = ocdbPath; }
  void SetRunList(Int_t runNumber);
  void SetRunList(const char* runlist);

  virtual void Scan(Int_t verbose=0) = 0;

  void Plot(const char* dcsname=0x0, Bool_t withPatch=kFALSE, Bool_t plotIntermediate=kFALSE);

  void Print(Option_t* dcsname="") const;

  TGraph* Combine(TObjArray& graphs);

  TMultiGraph* CombineMulti(TObjArray& graphs);

protected:

  TMap* CreateMap(Int_t runNumber, Bool_t patched=kFALSE) const;

  TMultiGraph* Map2Graph(TMap* m, const char* dcsname);

  bool IsAlmostEqualRelative(float A, float B, float maxRelDiff = FLT_EPSILON) const;

  void ReadIntegers(const char* filename, std::vector<int>& integers);

  TGraph* GraphValues(TMap* m, const char* name);

  void TimeAxis(TMultiGraph* g);

  AliMpDCSNamer* DCSNamer() const;

  AliMUONTrackerVoltages(const AliMUONTrackerVoltages& rhs); // not implemented on purpose
  AliMUONTrackerVoltages& operator=(const AliMUONTrackerVoltages& rhs); // not implemented on purpose

  std::vector<int> fRunList; // input run list
  TString fOCDBPath; // ocdb path (raw:// by default)
  mutable AliMpDCSNamer* fDCSNamer; // helper to name things
  TString fOCDBObjectPath; // path to the object to retrieve from the OCDB

  ClassDef(AliMUONTrackerVoltages,1) // Utility class to inspect MUON Tracker HV values
};

#endif
