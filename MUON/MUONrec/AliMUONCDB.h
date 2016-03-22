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
#include <vector>

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
class TH1;
class AliMergeableCollection;

namespace AliMUONCDB
{
  Bool_t CheckOCDB(Bool_t pathOnly = kFALSE);
  Bool_t CheckMapping(Bool_t segmentationOnly = kFALSE);

  Double_t MeanHVValueForDCSAlias(TMap& hvMap, const char* hvChannel);

  void CheckHV(Int_t runNumber, Int_t verbose=0);
  void CheckHV_ALIROOT_6402(const char* runlist, Bool_t verbose=kFALSE);
  Bool_t CheckHV_ALIROOT_6402(Int_t runNumber, Bool_t verbose=kFALSE);
  Bool_t IsSt1DCSAliasRemapped(const TString& name);
  void PatchHV(TMap& hvMap, TList* messages, Bool_t onlySt1remapped=kFALSE);

  Bool_t LoadField();
  Bool_t LoadMapping(Bool_t segmentationOnly = kFALSE);
  AliMUONRecoParam* LoadRecoParam();
  TClonesArray* LoadAlignmentData();

  void AddDCSValue ( TMap& aliasMap, Int_t imeas, const char* smt, const char* sInOut, Int_t rpc, Float_t value );

  Int_t MakeHVStore(TMap& aliasMap, Bool_t defaultValues);
  Int_t MakeLVStore(TMap& aliasMap, Bool_t defaultValues, time_t refTime);

  Int_t MakeTriggerDCSStore(TMap& aliasMap);
  Int_t MakePedestalStore(AliMUONVStore& pedestalStore, Bool_t defaultValues);
  Int_t MakeOccupancyMapStore(AliMUONVStore& occupancyMap, Bool_t defaultValues);
  Int_t MakeBusPatchEvolution(AliMergeableCollection& bpevo, int timeResolution=60);

  AliMUONRejectList* MakeRejectListStore(Bool_t defaultValues);

  Int_t MakeLocalTriggerMaskStore(AliMUONVStore& ltm);
  Int_t MakeRegionalTriggerConfigStore(AliMUONRegionalTriggerConfig& rtm);
  Int_t MakeGlobalTriggerConfigStore(AliMUONGlobalCrateConfig& gtm);

  AliMUONTriggerLut* MakeTriggerLUT(const char* file="$(ALICE_ROOT)/MUON/data/lutAptLpt1Hpt1p7.root");
  AliMUONTriggerEfficiencyCells* MakeTriggerEfficiency(const char* file="$ALICE_ROOT/MUON/data/efficiencyCells.dat");

  AliMUONVStore* Diff(AliMUONVStore& store1, AliMUONVStore& store2, const char* opt="abs");

  TH1** Plot(const AliMUONVStore& store, const char* name, Int_t nbins=512);

  void ReadIntegers(const char* filename, std::vector<int>& integers);

  void ShowConfig(Bool_t withStatusMap=kFALSE);

  void ShowFaultyBusPatches(const char* runlist,
                            double occLimit=0.1,
                            const char* outputBaseName="faulty.buspatches",
                            const char* ocdbPath="raw://");

  void WriteToCDB(const char* calibpath, TObject* object,
                  Int_t startRun, Int_t endRun, Bool_t defaultValues);
  void WriteToCDB(const char* calibpath, TObject* object,
                  Int_t startRun, Int_t endRun, const char* filename);
  void WriteToCDB(TObject* object, const char* calibpath, Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity(),
                  const char* comment="", const char* responsible="AliMUONCDB tester class");

  void WriteMapping(Int_t startRun=0,Int_t endRun=AliCDBRunRange::Infinity());

  void WriteTrigger(Bool_t defaultValues=kTRUE, Int_t startRun=0,Int_t endRun=AliCDBRunRange::Infinity());
  void WriteTracker(Bool_t defaultValues=kTRUE, Int_t startRun=0,Int_t endRun=AliCDBRunRange::Infinity());

  void WriteHV(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteHV(const char* inputFile, Int_t runNumber);

  void WriteLV(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity(), time_t refTime=1449969676);
  void WritePedestals(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteOccupancyMap(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteRejectList(Bool_t defaultValues, Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteConfig(Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteBPEVO(Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());

  void WriteLocalTriggerMasks(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteRegionalTriggerConfig(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteGlobalTriggerConfig(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());

  void WriteTriggerDCS(Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteTriggerLut(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
  void WriteTriggerEfficiency(Int_t startRun=0, Int_t endRun=AliCDBRunRange::Infinity());
}

#endif
