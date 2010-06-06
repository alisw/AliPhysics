#ifndef ALIMUONCDB_H
#define ALIMUONCDB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \namespace AliMUONCDB
/// \brief Helper functions to experience the OCDB
/// 
//  Author Laurent Aphecetche

#include <TString.h>
#include "AliCDBRunRange.h"

class AliMUONVStore;
class TMap;
class TClonesArray;
class AliMUONVCalibParam;
class AliMUONTriggerLut;
class AliMUONTriggerEfficiencyCells;
class AliMUONRegionalTriggerConfig;
class AliMUONGlobalCrateConfig;
class AliMUONRejectList;
class AliMUONRecoParam;

namespace AliMUONCDB
{
  Bool_t CheckOCDB(Bool_t pathOnly = kFALSE);
  Bool_t CheckMapping(Bool_t segmentationOnly = kFALSE);
  
  Bool_t LoadField();
  Bool_t LoadMapping(Bool_t segmentationOnly = kFALSE);
  AliMUONRecoParam* LoadRecoParam();
  TClonesArray* LoadAlignmentData();
  
  Int_t MakeNeighbourStore(AliMUONVStore& neighbourStore);

  Int_t MakeHVStore(TMap& aliasMap, Bool_t defaultValues);
  Int_t MakeTriggerDCSStore(TMap& aliasMap, Bool_t defaultValues);
  Int_t MakePedestalStore(AliMUONVStore& pedestalStore, Bool_t defaultValues);
  Int_t MakeCapacitanceStore(AliMUONVStore& capaStore, Bool_t defaultValues);
  Int_t MakeCapacitanceStore(AliMUONVStore& capaStore, const char* file);
  Int_t MakeGainStore(AliMUONVStore& gainStore, Bool_t defaultValues);
  Int_t MakeOccupancyMapStore(AliMUONVStore& occupancyMap, Bool_t defaultValues);
  AliMUONRejectList* MakeRejectListStore(Bool_t defaultValues);
  
  Int_t MakeLocalTriggerMaskStore(AliMUONVStore& ltm);  
  Int_t MakeRegionalTriggerConfigStore(AliMUONRegionalTriggerConfig& rtm);
  Int_t MakeGlobalTriggerConfigStore(AliMUONGlobalCrateConfig& gtm);
  
  AliMUONTriggerLut* MakeTriggerLUT(const char* file="$(ALICE_ROOT)/MUON/data/lutAptLpt1Hpt1p7.root");
  AliMUONTriggerEfficiencyCells* MakeTriggerEfficiency(const char* file="$ALICE_ROOT/MUON/data/efficiencyCells.dat");

  AliMUONVStore* Diff(AliMUONVStore& store1, AliMUONVStore& store2, const char* opt="abs");
  
  void Plot(const AliMUONVStore& store, const char* name, Int_t nbins=512);

  void ShowConfig();
  
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
  void WritePedestals(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteGains(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteCapacitances(Bool_t defaultValues, Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteCapacitances(const char* file, Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteOccupancyMap(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteRejectList(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteConfig(Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());

  void WriteLocalTriggerMasks(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteRegionalTriggerConfig(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteGlobalTriggerConfig(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  
  void WriteTriggerDCS(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteTriggerLut(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteTriggerEfficiency(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
}

#endif
