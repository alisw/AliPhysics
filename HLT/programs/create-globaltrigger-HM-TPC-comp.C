#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliHLTTriggerMenu.h"
#include "AliHLTReadoutList.h"
#include "AliHLTGlobalTriggerConfig.h"
#include "TObjString.h"
#include "TString.h"
#include "TSystem.h"
#include "Riostream.h"
using std::cerr;
using std::endl;
#endif

/**
 * Generates a default CDB entry for the trigger menu in the given CDB storage
 * \param cdbPath  The path to the default CDB storage.
 */
void create_globaltrigger_HM_TPC_comp(
		 const char* cdbPath = "local:///opt/HLT/HCDB",
		 Int_t version = -1,
		 Int_t firstRun = 0,
		 Int_t lastRun = AliCDBRunRange::Infinity()
		 )
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libAliHLTUtil.so");
  gSystem->Load("libAliHLTMUON.so");
  gSystem->Load("libAliHLTTRD.so");
  gSystem->Load("libAliHLTTrigger.so");
  
  // Setup the CDB default storage and run number.
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if (cdbManager == NULL) {
    cerr << "ERROR: Global CDB manager object does not exist." << endl;
    return;
  }

  AliCDBStorage* storage = cdbManager->GetStorage(cdbPath);
  if (storage == NULL) {
    cerr << "ERROR: Could not get storage for: " << cdbPath << endl;
    return;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  // Create the trigger menu:
  
  AliHLTGlobalTriggerConfig config("HM-PHYSICS-V0001");
  
  config.AddSymbol("domainALL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:***\")");
  config.AddSymbol("domainHLT", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:HLT ,DAQRDOUT:HLT\")");
  config.AddSymbol("domainNoTPC", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:***,-DAQRDOUT:TPC\")");
  
  config.AddItem(
		 3, // priority group.
		 "C0LSR-ABCE-NOPF-TPC",
		 "domainALL",
		 "H-TPC_LASER-V0001.001-ALL-ALL"
		 );

  config.AddItem(
		 3, // priority group.
		 "C0LSR-ABCE-NOPF-FASTNOTRD",
		 "domainALL",
		 "H-TPC_LASER-V0001.001-ALL-ALL"
		 );

  config.AddItem(
		 3, // priority group.
		 "C0LSR-ABCE-NOPF-CENTNOTRD",
		 "domainALL",
		 "H-TPC_LASER-V0001.001-ALL-ALL"
		 );

  config.AddItem(
		 2, // priority group.
		 "true",
		 "domainALL",
		 "H-MONITOR-V0001.001-ALL-ALL",
		 1.0  // scaled down to 1 % rate.
		 );
  /*
  config.AddItem(
		 1, // priority group.
		 "BarrelMultiplicityTrigger",
		 "domainNoTPC",
		 "H-BarrelMultiplicity-V0001.001-ALL-NoTPC"
		 );

  config.AddItem(
		 1, // priority group.
		 "ITSMultiplicityTrigger",
		 "domainNoTPC",
		 "H-ITSMultiplicityTrigger-V0001.001-ALL-NoTPC"
		 );
  */
  
  config.AddItem(
		 10, // priority group.
		 "START_OF_DATA",
		 "START_OF_DATA | domainHLT",
		 "Start-Of-Data"
		 );

  config.AddItem(
		 9, // priority group.
		 "END_OF_DATA",
		 "END_OF_DATA | domainHLT",
		 "End-Of-Data"
		 );

  config.AddItem(
		 8, // priority group.
		 "CALIBRATION",
		 "CALIBRATION | domainHLT",
		 "Calibration-Trigger"
		 );

  config.AddItem(
		 7, // priority group.
		 "SOFTWARE",
		 "SOFTWARE | domainHLT",
		 "Software-Trigger"
		 );

  ///////////////////////////////////////////////////////////////////////////////////////////
  // Setup defaults in case there is no global trigger.
  // For non-triggered events store everything from HLT and readout all detectors.
  config.SetDefaultTriggerDescription("No HLT global trigger");
  AliHLTTriggerDomain defaultDomain("*******:***,-DAQRDOUT:TPC");
  config.SetDefaultTriggerDomain(defaultDomain);
  config.SetDefaultResult(false);
  
  TObject* menu = AliHLTGlobalTriggerConfig::Menu()->Clone();
  menu->Print();

  ///////////////////////////////////////////////////////////////////////////////////////////
  // Write the trigger menu object to the CDB.
  AliCDBId id("HLT/ConfigHLT/HLTGlobalTrigger", firstRun, lastRun, version);
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("ALICE HLT alice-hlt-operations@cern.ch");
  metaData->SetComment("HM-PHYSICS-V0001");
  storage->Put(menu, id, metaData);
}

