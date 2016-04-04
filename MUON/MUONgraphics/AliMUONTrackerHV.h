#ifndef ALIMUONTRACKERHV_H
#define ALIMUONTRACKERHV_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

#ifndef ALIMUONTRACKERVOLTAGES_H
#include "AliMUONTrackerVoltages.h"
#endif

#include <vector>

class TMultiGraph;
class TMap;
class TGraph;
class AliMpDCSNamer;

class AliMUONTrackerHV : public AliMUONTrackerVoltages
{
public:

  AliMUONTrackerHV(const char* runlist, const char* ocdbPath="raw://");
  AliMUONTrackerHV(Int_t runNumber, const char* ocdbPath="raw://");
  virtual ~AliMUONTrackerHV();

  void HVoff(const char* logfile="lhc11de.log", const char* outputBaseName="hvoff");

  void Plot(const char* dcsname=0x0, Bool_t withPatch=kFALSE, Bool_t plotIntermediate=kFALSE);

  void ReportTrips(Bool_t includeLowOnes=kFALSE);

  void Scan(Int_t verbose=0);

  Int_t CompareMaps(const TMap& hv1, const TMap& hv2, Bool_t verbose=kFALSE) const;

private:

  TMap* CreateMap(Int_t runNumber, Bool_t patched) const;

  TMultiGraph* GraphHV(TMap* m, const char* dcsname);

  AliMUONTrackerHV(const AliMUONTrackerHV& rhs); // not implemented on purpose
  AliMUONTrackerHV& operator=(const AliMUONTrackerHV& rhs); // not implemented on purpose

  ClassDef(AliMUONTrackerHV,3) // Utility class to inspect MUON Tracker HV values
};

#endif
