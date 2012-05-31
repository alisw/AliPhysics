// $Id: $
/**
 * \ingroup trigger_menus
 * \file HM_TEST_V0001.C
 * \brief Macro for generating the test trigger menu for p+p triggering.
 *
 * This macro generates the HM-TEST-V0001 global trigger configuration.
 *
 * You can run this macro with defaults using the following shell command:
 * \code
 *   > aliroot -b -q $ALICE_ROOT/HLT/trigger/HM_TEST_V0001.C
 * \endcode
 *
 * This is a test setup of the trigger menu for p+p containing the following triggers:
 *  
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
void HM_TEST_V0001(
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
  
  AliHLTGlobalTriggerConfig config("HM-TEST-V0001");
  
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

  // DDLs for all detectors
  config.AddSymbol("domainALLDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:***\")");
  
  // DDLs for HLT only
  config.AddSymbol("domainHLTDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:HLT\")");
  
  // DDLs for the fast muon cluster detectors.
  config.AddSymbol("domainFastDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:MCH,DAQRDOUT:MTR,DAQRDOUT:SPD,DAQRDOUT:V00,DAQRDOUT:ZDC\")");
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  // Setup the trigger items in several priority groups.
  // First group (20) is for all the HLT triggers. The next 3 groups (12 - 10) are for
  // scaling down min-bias CTP interaction triggers but reading out only 50% of the HLT ESDs.
  // The last group (1) handles special software triggers.
  
  // At least one track in barrel with pT > 1 GeV/c, downscale by half.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-Barrel_pT_Single-V0001.001",
		 "domainHLTOUT | domainALLDDL",
		 2,  // scaledown factor 1/2
		 "H-Barrel_pT_Single-V0001.001-ALL-ALL"
		 );

  // At least one track in barrel with pT > 3 GeV/c, no downscale.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-Barrel_pT_Single-V0002.001",
		 "domainHLTOUT | domainALLDDL",
		 "H-Barrel_pT_Single-V0002.001-ALL-ALL"
		 );

  // At least one track in barrel with pT > 5 GeV/c, no downscale.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-Barrel_pT_Single-V0003.001",
		 "domainHLTOUT | domainALLDDL",
		 "H-Barrel_pT_Single-V0003.001-ALL-ALL"
		 );

  // Track multiplicity > 30 or 40, no downscale.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-Barrel_Multiplicity-V0001.001",
		 "domainHLTOUT | domainALLDDL",
		 "H-Barrel_Multiplicity-V0001.001-ALL-ALL"
		 );

  // At least one track in muon spectrometer, pT > 1 GeV/c, no downscale.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-MUON_pT_Single-V0001.001",
		 "domainHLTOUT | domainFastDDL",
		 "H-MUON_pT_Single-V0001.001-MUON-ALL"
		 );

  // At least one track in muon spectrometer, pT > 2 GeV/c, no downscale.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-MUON_pT_Single-V0001.002",
		 "domainHLTOUT | domainFastDDL",
		 "H-MUON_pT_Single-V0001.002-MUON-ALL"
		 );

  // Unlike sign track pair in muon spectrometer, at least one track pT > 1 GeV/c, no downscale.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-MUON_pT_UnlikePair-V0001.001",
		 "domainHLTOUT | domainFastDDL",
		 "H-MUON_pT_UnlikePair-V0001.001-MUON-ALL"
		 );

  // Unlike sign track pair in muon spectrometer, at least one track pT > 2 GeV/c, no downscale.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-MUON_pT_UnlikePair-V0001.002",
		 "domainHLTOUT | domainFastDDL",
		 "H-MUON_pT_UnlikePair-V0001.002-MUON-ALL"
		 );

  // Unlike sign track pair in muon spectrometer, invariant mass > 2.5 GeV/c, no downscale.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-MUON_Mass_UnlikePair-V0001.001",
		 "domainHLTOUT | domainFastDDL",
		 "H-MUON_Mass_UnlikePair-V0001.001-MUON-ALL"
		 );

  // Unlike sign track pair in muon spectrometer, invariant mass > 7 GeV/c, no downscale.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-MUON_Mass_UnlikePair-V0001.002",
		 "domainHLTOUT | domainFastDDL",
		 "H-MUON_Mass_UnlikePair-V0001.002-MUON-ALL"
		 );


  // Trigger on cluster energy in EMCAL.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-ClusterEnergy-V0001.001",
		 "domainHLTOUT | domainALLDDL",
		 "H-ClusterEnergy-V0001.001-ALL-ALL"
		 );

  // Trigger on cluster energy in PHOS.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-ClusterEnergy-V0001.002",
		 "domainHLTOUT | domainALLDDL",
		 "H-ClusterEnergy-V0001.002-ALL-ALL"
		 );

  // Trigger based on invariant mass cut for D0.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-D0_Mass-V0001.001",
		 "domainHLTOUT | domainALLDDL",
		 "H-D0_Mass-V0001.001-ALL-ALL"
		 );

  // Trigger using kT algorithm.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-JET_ET-V0001.001",
		 "domainHLTOUT | domainALLDDL",
		 "H-JET_ET-V0001.001-ALL-ALL"
		 );

  // Trigger using anti-kT algorithm.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-JET_ET-V0002.001",
		 "domainHLTOUT | domainALLDDL",
		 "H-JET_ET-V0002.001-ALL-ALL"
		 );

  // Trigger using fixed seed cone algorithm.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-JET_ET-V0003.001",
		 "domainHLTOUT | domainALLDDL",
		 "H-JET_ET-V0003.001-ALL-ALL"
		 );

  // Trigger if kT and anti-kT triggered.
  config.AddItem(
		 20, // priority group.
		 "(CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD) && H-JET_ET-V0004.001",
		 "domainHLTOUT | domainALLDDL",
		 "H-JET_ET-V0004.001-ALL-ALL"
		 );

  // Scaled down min bias trigger.
  config.AddItem(
		 12, // priority group.
		 "CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD",
		 "domainHLTOUT | domainALLDDL",
		 2,  // scaledown factor 1/2
		 "H-MinBias_Scale_Down-V0001.001-ALL-ALL"
		 );

  // Readout only 50% of HLT ESDs for min bias.
  config.AddItem(
		 11, // priority group.
		 "CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD",
		 "domainESD | domainHLTDDL",
		 2,  // scaledown factor 1/2
		 "Rejected min-bias with HLT ESD readout",
		 false  // default global trigger decision result
		 );

  // Reject completely the other 50% min bias.
  config.AddItem(
		 10, // priority group.
		 "CINT1WU-B-NOPF-ALL || CINT1-B-NOPF-ALLNOTRD",
		 "domainHLTDDL",  // Only HLT DDL to deliver at least the decision.
		 0,  // no scaledown
		 "Rejected min-bias",
		 false  // default global trigger decision result
		 );

  config.AddItem(
		 1, // priority group.
		 "SOFTWARE || CALIBRATION || START_OF_DATA || END_OF_DATA",
		 "domainHLTOUT | domainALLDDL",
		 "H-SoftwareTrigger-V0001.001-ALL-ALL"
		 );

  ///////////////////////////////////////////////////////////////////////////////////////////
  // Setup defaults in case there is no global trigger.
  // For non-triggered events always readout all detectors. This is the catch all for rare
  // and background CTP triggers.
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
  metaData->SetComment("HM-TEST-V0001");
  storage->Put(menu, id, metaData);
}
