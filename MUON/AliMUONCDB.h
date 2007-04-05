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

class TList;
class AliMUONV1DStore;
class AliMUONV2DStore;
class TMap;
class AliMUONVCalibParam;
class AliMUONTriggerLut;
class AliMUONTriggerEfficiencyCells;

#define ALIMUONCDBINFINITY 99999999

class AliMUONCDB : public TObject
{
public:
  /// Ctor. change the path for testing the Shuttle preprocessor, to
  /// "local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB"
  AliMUONCDB(const char* cdbpath = "local://$ALICE_ROOT");
  virtual ~AliMUONCDB();
  
  Int_t MakeNeighbourStore(AliMUONV2DStore& neighbourStore);

  Int_t MakeHVStore(TMap& aliasMap, Bool_t defaultValues);
  Int_t MakePedestalStore(AliMUONV2DStore& pedestalStore, Bool_t defaultValues);
  Int_t MakeCapacitanceStore(AliMUONV1DStore& capaStore, Bool_t defaultValues);
  Int_t MakeGainStore(AliMUONV2DStore& gainStore, Bool_t defaultValues);
  
  Int_t MakeLocalTriggerMaskStore(AliMUONV1DStore& ltm) const;
  Int_t MakeRegionalTriggerMaskStore(AliMUONV1DStore& rtm) const;
  Int_t MakeGlobalTriggerMaskStore(AliMUONVCalibParam& gtm) const;
  AliMUONTriggerLut* MakeTriggerLUT(const char* file="$(ALICE_ROOT)/MUON/data/lutAptLpt1Hpt1p7.root") const;
  AliMUONTriggerEfficiencyCells* MakeTriggerEfficiency(const char* file="$ALICE_ROOT/MUON/data/efficiencyCells.dat") const;

  /// Compute the difference between two (compatible) stores
  AliMUONV2DStore* Diff(AliMUONV2DStore& store1, AliMUONV2DStore& store2, 
                        const char* opt="abs");
    
  void Plot(const AliMUONV2DStore& store, const char* name, Int_t nbins=512);

  void WriteToCDB(const char* calibpath, TObject* object, 
                  Int_t startRun, Int_t endRun, Bool_t defaultValues);

  void WriteTrigger(Int_t startRun=0,Int_t endRun=ALIMUONCDBINFINITY);
  void WriteTracker(Bool_t defaultValues=kTRUE, Int_t startRun=0,Int_t endRun=ALIMUONCDBINFINITY);
  
  void WriteNeighbours(Int_t startRun=0, Int_t endRun=ALIMUONCDBINFINITY);
  void WriteHV(Bool_t defaultValues, Int_t startRun, Int_t endRun);
  void WritePedestals(Bool_t defaultValues, Int_t startRun, Int_t endRun=ALIMUONCDBINFINITY);
  void WriteGains(Bool_t defaultValues, Int_t startRun, Int_t endRun=ALIMUONCDBINFINITY);
  void WriteCapacitances(Bool_t defaultValues, Int_t startRun=0, Int_t endRun=ALIMUONCDBINFINITY);
  
  void WriteLocalTriggerMasks(Int_t startRun=0, Int_t endRun=ALIMUONCDBINFINITY);
  void WriteRegionalTriggerMasks(Int_t startRun=0, Int_t endRun=ALIMUONCDBINFINITY);
  void WriteGlobalTriggerMasks(Int_t startRun=0, Int_t endRun=ALIMUONCDBINFINITY);
  void WriteTriggerLut(Int_t startRun=0, Int_t endRun=ALIMUONCDBINFINITY);
  void WriteTriggerEfficiency(Int_t startRun=0, Int_t endRun=ALIMUONCDBINFINITY);
  
private:
    AliMUONCDB(const AliMUONCDB& rhs);
  AliMUONCDB& operator=(const AliMUONCDB& rhs);
  
private:
  TString fCDBPath; //!< where to write stuff
  TList* fManuList; //!< full list of manus
  static const Int_t fgkMaxNofChannelsPerManu; //!< 64
  
  ClassDef(AliMUONCDB,0) // Helper class to experience OCDB
};

#endif
