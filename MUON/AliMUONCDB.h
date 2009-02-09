#ifndef ALIMUONCDB_H
#define ALIMUONCDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup base
/// \class AliMUONCDB
/// \brief Helper class to experience the OCDB
/// 
//  Author Laurent Aphecetche

#include <TObject.h>
#include <TString.h>
#include "AliCDBRunRange.h"

class AliMUONVStore;
class TMap;
class AliMUONVCalibParam;
class AliMUONTriggerLut;
class AliMUONTriggerEfficiencyCells;
class AliMUONRegionalTriggerConfig;
class AliMUONGlobalCrateConfig;

class AliMUONCDB : public TObject
{
public:
  /// Ctor. change the path for testing the Shuttle preprocessor, to
  /// "local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB"
  AliMUONCDB(const char* cdbpath = "local://$ALICE_ROOT/OCDB");
  virtual ~AliMUONCDB();
  
  void SetMaxNofChannelsToGenerate(Int_t n);
  
  Int_t MakeNeighbourStore(AliMUONVStore& neighbourStore);

  Int_t MakeHVStore(TMap& aliasMap, Bool_t defaultValues);
  Int_t MakeTriggerDCSStore(TMap& aliasMap, Bool_t defaultValues);
  Int_t MakePedestalStore(AliMUONVStore& pedestalStore, Bool_t defaultValues);
  Int_t MakeCapacitanceStore(AliMUONVStore& capaStore, Bool_t defaultValues);
  Int_t MakeCapacitanceStore(AliMUONVStore& capaStore, const char* file);
  Int_t MakeGainStore(AliMUONVStore& gainStore, Bool_t defaultValues);
  
  Int_t MakeLocalTriggerMaskStore(AliMUONVStore& ltm) const;  
  Int_t MakeRegionalTriggerConfigStore(AliMUONRegionalTriggerConfig& rtm) const;
  Int_t MakeGlobalTriggerConfigStore(AliMUONGlobalCrateConfig& gtm) const;
  
  AliMUONTriggerLut* MakeTriggerLUT(const char* file="$(ALICE_ROOT)/MUON/data/lutAptLpt1Hpt1p7.root") const;
  AliMUONTriggerEfficiencyCells* MakeTriggerEfficiency(const char* file="$ALICE_ROOT/MUON/data/efficiencyCells.dat") const;

  /// Compute the difference between two (compatible) stores
  AliMUONVStore* Diff(AliMUONVStore& store1, AliMUONVStore& store2, 
                        const char* opt="abs");
    
  void Plot(const AliMUONVStore& store, const char* name, Int_t nbins=512);

  void WriteToCDB(const char* calibpath, TObject* object, 
                  Int_t startRun, Int_t endRun, Bool_t defaultValues);
  void WriteToCDB(const char* calibpath, TObject* object, 
                  Int_t startRun, Int_t endRun, const char* filename);
  void WriteToCDB(TObject* object, const char* calibpath, Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity(),
		  const char* comment="", const char* responsible="AliMUONCDB tester class");

  void WriteTrigger(Bool_t defaultValues=kTRUE, Int_t startRun=0,Int_t endRun=AliCDBRunRange::Infinity());
  void WriteTracker(Bool_t defaultValues=kTRUE, Int_t startRun=0,Int_t endRun=AliCDBRunRange::Infinity());
  
  void WriteNeighbours(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteHV(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteTriggerDCS(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WritePedestals(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteGains(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteCapacitances(Bool_t defaultValues, Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteCapacitances(const char* file, Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  
  void WriteLocalTriggerMasks(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteRegionalTriggerConfig(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteGlobalTriggerConfig(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  
  void WriteTriggerLut(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteTriggerEfficiency(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  
private:
  /// Not implemented
  AliMUONCDB(const AliMUONCDB& rhs);
  /// Not implemented
  AliMUONCDB& operator=(const AliMUONCDB& rhs);
  
private:
  TString fCDBPath; //!< where to write stuff
  Int_t fMaxNofChannelsToGenerate; //!< to limit the number of generated channels (debug)
  
  ClassDef(AliMUONCDB,0) // Helper class to experience OCDB
};

#endif
