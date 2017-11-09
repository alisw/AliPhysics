/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

#ifndef ALIJRUNTABLE_H
#define ALIJRUNTABLE_H

#include <TString.h>
#include <vector>
class AliJRunTable {
public:
  enum { kUnknownPeriod, kLHC10b, kLHC10c, kLHC10d, kLHC10e,
    kLHC10h, kLHC11a,kLHC11b , kLHC11c , kLHC11d, kLHC11e, kLHC11h, kLHC12g,kLHC12h,
    kLHC13b, kLHC13c, kLHC13d, kLHC13e, kLHC13g, kLHC15o, kLHC16g1, kLHC16q, kLHC17f,
    kJNPeriod};
  enum { kRE, kMC };
  enum { kPP, kPbPb, kPA };
  AliJRunTable();
  void Init();
  void AddGoodRun(int period, int runnumber){ fGoodRun[period].push_back(runnumber); }
  TString GetPeriodName( int period=-1 ) const;
  int GetBeamType( int period ) const { return fBeamType[period]; }
  TString GetMCPeriod( int period ) const { return fMCPeriod[period]; }
  int GetPeriodCode(TString perstr) const;
  int GetRunNumberToPeriod( int runnumber );
  void SetPeriodInformation(int period, TString name, int beamtype, int datatype, int energy, int run0, int run1, TString MCPeriod);
  int GetBeamEnergy(int period) const { return fEnergy[period]; }
  
  int GetRunNumberFromString( const char * tstr );
  TString GetPeriodFromString( const char * tstr ) const;
  TString GetMCPeriodFromString( const char * tstr) const;
  int GetPeriodFromMCPeriod( const char * tstr );
  const char * GetBeamStr( int ib=-1 ) const;  // Never used anywhere, should it be removed?
  const char * GetDataStr( int ib=-1 ) const;  // Never used anywhere, should it be removed?
  const char * GetEnergyStr( int ib=-1 ) const; // Never used anywhere, should it be removed?
  bool ParseString( const char * tstr );
  
  int GetPeriod() const { return fCPeriod; }
  TString GetPeriodMCName() const { return fCPeriodMCName; }
  int GetBeamType() const { return fBeamType[fCPeriod]; }
  TString GetBeamTypeStr() const { return IsPP()?"pp":IsPA()?"pA":IsHeavyIon()?"PbPb":""; }
  void SetRunNumber( int run ){ fCRunNumber=run;fCPeriod=GetRunNumberToPeriod(run); }
  int GetRunNumber() const { return fCRunNumber; }
  
  bool IsHeavyIon() const { return fBeamType[fCPeriod]==kPbPb; }
  bool IsPA() const { return fBeamType[fCPeriod]==kPA; }
  bool IsPP() const { return fBeamType[fCPeriod]==kPP; }
  
  AliJRunTable *  operator()( int run ){ SetRunNumber(run);return this; }
  static AliJRunTable& GetSpecialInstance();
  static const AliJRunTable& GetInstance();
  
private:
  int fBeamType[kJNPeriod];  // comment needed
  int fDataType[kJNPeriod];  // comment needed
  int fEnergy[kJNPeriod]; // comment needed
  long fRunRange[kJNPeriod][2]; // comment needed
  TString fMCPeriod[kJNPeriod]; // comment needed
  TString fPeriodName[kJNPeriod]; // comment needed
  std::vector<int> fGoodRun[kJNPeriod]; // comment needed
  
  int fCPeriod; // comment needed
  int fCRunNumber; // comment needed
  TString fCPeriodMCName; // comment needed
};

#endif
