// $Id: $
/**
 * \ingroup trigger_menus
 * \file HM_PHYSICS_V0002.C
 * \brief Macro for generating the trigger menu for min-bias p+p triggering.
 *
 * This macro generates the HM-PHYSICS-V0002 global trigger configuration.
 *
 * You can run this macro with defaults using the following shell command:
 * \code
 *   > aliroot -b -q $ALICE_ROOT/HLT/trigger/HM_PHYSICS_V0002.C
 * \endcode
 *
 * Configuration:
 *  Only triggering on min-bias CTP triggers, pass through all other triggers.
 *  Scale down the trigger rate and store HLT ESD for rejected events.
 *  Triggering on high pT tracks as defined by CDB entries:
 *   H_._Barrel_pT_Single_._V0001.001
 *   H_._Barrel_pT_Single_._V0002.001
 *   H_._Barrel_pT_Single_._V0003.001
 *
 * \note The above mentioned CDB entries must already exist. They are not created
 *    by this macro.
 *
 * \author Artur Szostak <artursz@iafrica.com>
 */

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
 * (local by default).
 * \param cdbPath  The path to the default CDB storage.
 */
void HM_PHYSICS_V0002(
		      const char* cdbPath = "local://$ALICE_ROOT/OCDB",
		      Int_t version = 0,
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
  
  AliHLTGlobalTriggerConfig config("HM-PHYSICS-V0002");
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  // Predefine some explicit domain definitions used in the trigger menu.
  config.AddSymbol("domainTObject"   , "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"ROOTTOBJ:HLT \")");
  config.AddSymbol("domainESD"       , "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"ALIESDV0:HLT \")");
  config.AddSymbol("domainHistogram" , "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"ROOTHIST:HLT \")");
  config.AddSymbol("domainSPDCluster", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"CLUSTERS:ISPD\")");
  config.AddSymbol("domainSDDCluster", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"CLUSTERS:ISDD\")");
  config.AddSymbol("domainSSDCluster", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"CLUSTERS:ISSD\")");
  config.AddSymbol("domainTPCCluster", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"CLUSTERS:TPC \")");

  config.AddSymbol("domainHLTOUT", "AliHLTTriggerDomain", "", 
		   "domainTObject    | "
		   "domainESD        | "
		   "domainHistogram  | "
		   "domainSPDCluster | "
		   "domainSDDCluster | "
		   "domainSSDCluster | "
		   "domainTPCCluster"
		   );

  config.AddSymbol("domainALLDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:***\")");
  config.AddSymbol("domainHLTDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:HLT\")");
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  // Setup the trigger items in 5 priority groups.
  // Group 5 is for the HLT triggers.
  // Group 4 is for scaling down min-bias CTP interaction triggers.
  // Group 2 and 3 are used to readout only HLT ESDs for rejected min bias events.
  // The last group (1) handles special software triggers.
  
  config.AddItem(
		 5, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-Barrel_pT_Single-V0001.001",
		 "domainHLTOUT | domainALLDDL", 
		 5,  // scaledown factor 1/5
		 "H-Barrel_pT_Single-V0001.001-ALL-ALL"
		 );

  config.AddItem(
		 5, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-Barrel_pT_Single-V0002.001",
		 "domainHLTOUT | domainALLDDL",
		 2,  // scaledown factor 1/2
		 "H-Barrel_pT_Single-V0002.001-ALL-ALL"
		 );

  config.AddItem(
		 5, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-Barrel_pT_Single-V0003.001",
		 "domainHLTOUT | domainALLDDL", 
		 "H-Barrel_pT_Single-V0003.001-ALL-ALL"
		 );

  config.AddItem(
		 4, // priority group.
		 "CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD",
		 "domainHLTOUT | domainALLDDL",
		 "H-MINBIAS_SCALE_DOWN-V0003.001-CENTRAL-ALL",
		 15.  // scaledown factor 0.15
		 );

  // Readout only 50% of HLT ESDs for min bias.
  config.AddItem(
		 3, // priority group.
		 "CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD",
		 "domainESD | domainHLTDDL",
		 2,  // scaledown factor 1/2
		 "Rejected min-bias with HLT ESD readout",
		 false  // default global trigger decision result
		 );

  // Reject completely the other 50% min bias.
  config.AddItem(
		 2, // priority group.
		 "CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD",
		 "domainHLTDDL",  // Only HLT DDL to deliver at least the decision.
		 0,  // No prescalar (no scaledown)
		 "Rejected min-bias",
		 false  // default global trigger decision result
		 );

  config.AddItem(
		 1, // priority group.
		 "SOFTWARE || CALIBRATION || START_OF_DATA || END_OF_DATA",
		 "SOFTWARE | CALIBRATION | START_OF_DATA | END_OF_DATA | domainHLTDDL",
		 "H-SoftwareTrigger-V0001.001-ALL-ALL"
		 );

  ///////////////////////////////////////////////////////////////////////////////////////////
  // Setup defaults in case there is no global trigger.
  // For non-triggered events store everything from HLT and readout all detectors.
  config.SetDefaultTriggerDescription("No HLT global trigger");
  AliHLTTriggerDomain defaultDomain("*******:***");
  AliHLTReadoutList readoutlist;
  readoutlist.Enable(AliHLTReadoutList::kALLDET);
  defaultDomain.Add(readoutlist);
  config.SetDefaultTriggerDomain(defaultDomain);
  config.SetDefaultResult(true);
  
  TObject* menu = AliHLTGlobalTriggerConfig::Menu()->Clone();
  menu->Print();
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  // Write the trigger menu object to the CDB.
  AliCDBId id("HLT/ConfigHLT/HLTGlobalTrigger", firstRun, lastRun, version);
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("ALICE HLT Artur.Szostak@cern.ch");
  metaData->SetComment("HM-PHYSICS-V0002");
  storage->Put(menu, id, metaData);
}
