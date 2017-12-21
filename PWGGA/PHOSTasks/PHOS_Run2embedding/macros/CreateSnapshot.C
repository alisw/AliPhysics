#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TString.h>
#include <TIterator.h>
#include <TGrid.h>
#include <AliRawReader.h>
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliLog.h"
#include "TStopwatch.h"
#endif

void CreateSnapshot(const char* snapshotName=0, const char* rawdata=0);


const char* kCDBExclude[] = {
  "EMCAL/Config/Preprocessor"
  ,"EMCAL/Config/Temperature"
  ,"TPC/Config/Goofie"
  ,"TPC/Config/HighVoltage"
  ,"TPC/Config/HighVoltageStat"
  ,"TPC/Config/Preprocessor"
  ,"TPC/Config/Temperature"
  ,"HLT/ConfigHLT/BarrelHighMultiplicity"	
  ,"HLT/ConfigHLT/BarrelMultiplicityTrigger"	
  ,"HLT/ConfigHLT/BarrelPt_v01"	
  ,"HLT/ConfigHLT/BarrelPt_v02"	
  ,"HLT/ConfigHLT/BarrelPt_v03"	
  ,"HLT/ConfigHLT/EmcalClusterEnergyTrigger"	
  ,"HLT/ConfigHLT/esdLayout"	
  ,"HLT/ConfigHLT/GlobalOfflineVertexer"	
  ,"HLT/ConfigHLT/GlobalTrackMatcher"	
  ,"HLT/ConfigHLT/HLTGlobalTrigger"	
  ,"HLT/ConfigHLT/ITSMultiplicityTrigger"	
  ,"HLT/ConfigHLT/PhosClusterEnergyTrigger"	
  ,"HLT/ConfigHLT/PhosMipTrigger"	
  ,"HLT/ConfigHLT/PrimaryVertexFinder"	
  ,"HLT/ConfigHLT/V0Finder"	
  ,"HLT/ConfigITS/ITSClusterFinderSDD"	
  ,"HLT/ConfigITS/ITSClusterFinderSPD"	
  ,"HLT/ConfigITS/ITSClusterFinderSSD"	
  ,"HLT/ConfigITS/ITSTracker"	
  ,"HLT/ConfigMUON/DecisionComponent"	
  ,"HLT/ConfigMUON/FieldIntegrals"	
  ,"HLT/ConfigMUON/HitReconstructor"	
  ,"HLT/ConfigMUON/MansoTrackerFSM"	
  ,"HLT/ConfigMUON/TriggerReconstructor"	
  ,"HLT/ConfigPHOS/PHOSClusterizer"	
  ,"HLT/ConfigPHOS/PHOSDigitMaker"	
  ,"HLT/ConfigPHOS/PHOSRawAnalyzer"	
  ,"HLT/ConfigSample/SampleComponent1"	
  ,"HLT/ConfigSample/SampleESDAnalysis"	
  ,"HLT/ConfigTPC/TPCCAGlobalMerger"	
  ,"HLT/ConfigTPC/TPCCATracker"	
  ,"HLT/ConfigTPC/TPCDataCompressor"	
  ,"HLT/ConfigTPC/TPCDataCompressorHuffmanTables"	
  ,"HLT/ConfigTPC/TPCHWClusterFinder"	
  ,"HLT/ConfigTPC/TPCHWClusterTransform"	
  ,"HLT/ConfigTPC/TPCTrackHisto"	
  ,"HLT/ConfigVZERO/VZEROReconstruction"	
  ,"HLT/ConfigZDC/ZDCESDReco"
};

Bool_t IsExcluded(const char* objName)
{
  // check if the object should be excluded from the snapshot
  const int kNExclude = sizeof(kCDBExclude)/sizeof(char*);
  TString nameS = objName;
  for (int i=kNExclude;i--;) if (nameS==kCDBExclude[i]) return kTRUE;
  return kFALSE;
}

const Char_t *snapshotName[2] = {
  "OCDBsim.root",
  "OCDBrec.root"
};

CreateSnapshot(Int_t mode)
{

  gROOT->LoadMacro("Sim/OCDBConfig.C");
  
  // run number
  Int_t runNumber = -1;
  if (gSystem->Getenv("CONFIG_RUN"))
    runNumber = atoi(gSystem->Getenv("CONFIG_RUN"));
  if (runNumber <= 0) {
    printf("Invalid run number: %d \n", runNumber);
    abort();
  }

  OCDBConfig(kOCDBDefault, runNumber, mode);
  CreateSnapshot(snapshotName[mode]);
  
}

void CreateSnapshot(const char* snapshotName, const char* rawdata)
{
  TStopwatch sw;
  sw.Start();
  AliCDBManager* man = AliCDBManager::Instance();
  //
  // Make sure the run number is set
  int run = man->GetRun();
  if (run<0) {
    printf("I-CreateSnapshot: Run number is not set in CDBManager\n");
    if (gSystem->Getenv("DC_RUN")) {
      run = atoi(gSystem->Getenv("DC_RUN"));
      printf("I-CreateSnapshot: Run number %d extracted from $DC_RUN env.var\n",run);
    }
    else {
      TString rawdataS = rawdata;
      if (rawdataS.IsNull()) {
	printf("F-CreateSnapshot: No raw data ir $DC_RUN provided to extract run number\n");
	exit(1);
      }
      if (!gGrid && rawdataS.BeginsWith("alien://")) {
	TGrid::Connect("alien://");
	if (!gGrid) {
	  printf("F-CreateSnapshot: Failed to connect to alien\n"); 
	  exit(1);
	}
      }
      AliRawReader* reader = AliRawReader::Create(rawdataS.Data());
      if (!reader) {
	printf("F-CreateSnapshot: Failed to open %s\n",rawdataS.Data()); 
	exit(1);
      }
      reader->NextEvent();
      int run = reader->GetRunNumber();   
      delete reader;
      printf("I-CreateSnapshot: Extracted run number %d from %s\n",run,rawdataS.Data());
    }
    if (run<0) {
      printf("F-CreateSnapshot: Attempts to extract positve run number failed\n");
      exit(1);
    }
    man->SetRun(run);
  }
  //
  if (!man->IsDefaultStorageSet()) {
    printf("F-Default storage is not set\n");
    exit(1);
  }
  //
  AliCDBStorage *defStorage = man->GetDefaultStorage();
  TObjArray* arrCDBID = defStorage->GetQueryCDBList();
  AliCDBId* cdbID = 0;
  const TMap* stMap = man->GetStorageMap();
  man->SetCacheFlag(kTRUE);
  //
  printf("I-CreateSnapshot: Processing default storage\n");
  TIter nxt(arrCDBID);
  while ((cdbID=(AliCDBId*)nxt())) { // loop over default storage
    TString path = cdbID->GetPath();
    if (stMap->GetValue(path)) {
      printf("I-CreateSnapshot: %s has specific storage set, skip it in DefaultStorage\n",
	     path.Data());
      continue; // defined in the specific storage
    }
    if (IsExcluded(path.Data())) {
      printf("I-CreateSnapshot: object %s is in the exclusion list\n",path.Data());
      continue;
    }
    man->Get(path.Data());
  }
  // 
  TIter nextSt(stMap);
  TObjString *str;
  printf("CreateSnapshot: Processing specific storages\n");
  while ((str=(TObjString*)nextSt())) { // exclusion is not applied to specific objects
    TString calType = str->GetString();
    if (calType=="default") continue;
    man->Get(calType.Data());
    //
  }
  //
  TString snapshotNameS = snapshotName;
  if (snapshotNameS.IsNull()) {
    snapshotNameS = getenv("OCDB_SNAPSHOT_FILENAME");
    if (snapshotNameS.IsNull()) snapshotNameS = "OCDB";
  }
  if (!snapshotNameS.EndsWith(".root")) snapshotNameS += ".root";
  man->DumpToSnapshotFile(snapshotNameS.Data(),kFALSE);
  //
  sw.Stop();
  sw.Print();
}
