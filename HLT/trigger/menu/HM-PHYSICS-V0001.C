// $Id$
/**
 * @ingroup trigger_manus
 * @file HM_PHYSICS_V0001.C
 * @brief Macro for generating the trigger menu CDB entry for the HLT global trigger.
 *
 * This macro generates the HM-PHYSICS-V0001 global trigger configuration
 * and the trigger component configurations.
 * It is an experts macro so make sure you know what you are doing.
 *
 * You can run this macro with defaults using the following shell command:
 * @code
 *   > aliroot -b -q $ALICE_ROOT/HLT/trigger/HM_PHYSICS_V0001.C
 * @endcode
 *
 * Initial configuration:
 *  1. Trigger items based on TPC tracks
 *  ------------------------------------
 *  a) H-TRACK_MULTIPLICITY-V0002.001-CENTRAL-ALL   min 10 tracks
 *     Trigger name: BarrelMultiplicityTrigger
 *  b) H-TRACK_MULTIPLICITY-V0002.002-CENTRAL-ALL   min 100 tracks
 *     Trigger name: BarrelHighMultiplicity
 *  c) H-TRACK_MULTIPLICITY-V0002.003-CENTRAL-ALL   min 1 tracks pt > 0.5 GeV
 *     Trigger name: BarrelPt_v01
 *  d) H-TRACK_MULTIPLICITY-V0002.004-CENTRAL-ALL   min 1 tracks pt > 1 GeV
 *     Trigger name: BarrelPt_v02
 *  e) H-TRACK_MULTIPLICITY-V0002.005-CENTRAL-ALL   min 1 tracks pt > 5 GeV
 *     Trigger name: BarrelPt_v03
 *
 *  2. Min bias trigger
 *  -----------------------------
 *  every 100th event
 *
 * ChangeLog:
 *  2009-12-04    adding explicite trigger domain entries for blocks to be
 *                included in 'domainHLTOUT' in order to avoid wildcards in
 *                the definition of domain entries
 * @author Matthias.Richter@ift.uib.no
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliHLTTriggerMenu.h"
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
void HM_PHYSICS_V0001(
		      const char* cdbPath = "local://$ALICE_ROOT/OCDB",
		      Int_t version = 0,
		      Int_t firstRun = 0,
		      Int_t lastRun = AliCDBRunRange::Infinity()
		      ) {
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

  // /////////////////////////////////////////////////////////////////////////////////////////	
  // Create the trigger menu.
  AliHLTGlobalTriggerConfig config("HM-COSMICS-V0001");

  config.AddSymbol("domainAll", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:***\")");
  
  // /////////////////////////////////////////////////////////////////////////////////////////	
  // the domain definitions for the global HLT output and the HLT DDLs
  //config.AddSymbol("domainHLTOUT", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:HLT \")");
  config.AddSymbol("domainHLTDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:HLT\\0\")");

  // some explicite domain entries
  config.AddSymbol("domainTObject"   , "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"ROOTTOBJ:HLT \")");
  config.AddSymbol("domainESD"       , "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"ALIESDV0:HLT \")");
  config.AddSymbol("domainHistogram" , "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"ROOTHIST:HLT \")");

  config.AddSymbol("domainSPDCluster", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"CLUSTERS:ISPD\")");
  config.AddSymbol("domainSDDCluster", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"CLUSTERS:ISDD\")");
  config.AddSymbol("domainSSDCluster", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"CLUSTERS:ISSD\")");
  config.AddSymbol("domainTPCCluster", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"CLUSTERS:TPC \")");

  // an explicite HLTOUT selection which avoids wildcards
  config.AddSymbol("domainHLTOUT", "AliHLTTriggerDomain", "", 
		   "domainTObject    | "
		   "domainESD        | "
		   "domainHistogram  | "
		   "domainSPDCluster | "
		   "domainSDDCluster | "
		   "domainSSDCluster | "
		   "domainTPCCluster"
		   );

  // /////////////////////////////////////////////////////////////////////////////////////////	
  // -- DETECTOR READOUT DOMAINS
  config.AddSymbol("domainSPDDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:SPD\\0\")");
  config.AddSymbol("domainSDDDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:SDD\\0\")");
  config.AddSymbol("domainSSDDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:SSD\\0\")");
  config.AddSymbol("domainTPCDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:TPC\\0\")");
  config.AddSymbol("domainTRDDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:TRD\\0\")");
  config.AddSymbol("domainTOFDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:TOF\\0\")");
  config.AddSymbol("domainHMPDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:HMP\\0\")");
  config.AddSymbol("domainPHSDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:PHS\\0\")");
  config.AddSymbol("domainCPVDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:CPV\\0\")");
  config.AddSymbol("domainPMDDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:PMD\\0\")");
  config.AddSymbol("domainMCHDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:MCH\\0\")");
  config.AddSymbol("domainMTRDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:MTR\\0\")");
  config.AddSymbol("domainFMDDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:FMD\\0\")");
  config.AddSymbol("domainT00DDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:T00\\0\")");
  config.AddSymbol("domainV00DDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:V00\\0\")");
  config.AddSymbol("domainZDCDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:ZDC\\0\")");
  config.AddSymbol("domainACODDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:ACO\\0\")");
  config.AddSymbol("domainCTPDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:TRI\\0\")");
  config.AddSymbol("domainEMCDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:EMC\\0\")");

  // /////////////////////////////////////////////////////////////////////////////////////////	
  // -- DETECTOR READOUT DOMAINS - SPECIAL
  config.AddSymbol("domainALLDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:***\\0\")");
  
  // /////////////////////////////////////////////////////////////////////////////////////////	
  // NOTE: always make sure that the global HLT output and the HLT DDLs are included
  // in the readout, i.e. add domainHLTOUT|domainHLTDDL to the trigger domain

  // -- BarrelMultiplicity
  config.AddItem(
		 "BarrelMultiplicityTrigger", 
		 "domainHLTOUT|domainALLDDL", 
		 "H-TRACK_MULTIPLICITY-V0002.001-CENTRAL-ALL"
		 );

  config.AddItem(
		 "BarrelHighMultiplicity", 
		 "domainHLTOUT|domainALLDDL", 
		 "H-TRACK_MULTIPLICITY-V0002.002-CENTRAL-ALL"
		 );

  config.AddItem(
		 "BarrelPt_v01", 
		 "domainHLTOUT|domainALLDDL", 
		 "H-TRACK_MULTIPLICITY-V0002.003-CENTRAL-ALL"
		 );

  config.AddItem(
		 "BarrelPt_v01", 
		 "domainHLTOUT|domainALLDDL", 
		 "H-TRACK_MULTIPLICITY-V0002.004-CENTRAL-ALL"
		 );

  config.AddItem(
		 "BarrelPt_v01", 
		 "domainHLTOUT|domainALLDDL", 
		 "H-TRACK_MULTIPLICITY-V0002.005-CENTRAL-ALL"
		 );

  // -- Min bias trigger
  config.AddItem(
		 "true",
		 "domainALLDDL|domainHLTOUT", 
		 10, "H-MINBIAS_SCALE_DOWN-V0002.001-CENTRAL-ALL"
		 );

  
  // /////////////////////////////////////////////////////////////////////////////////////////	
  // default domain in case there is no global trigger
  // readout the output of the reconstruction
  // this refers to the domain domainHLTOUT|domainHLTDDL
  config.SetDefaultTriggerDescription("No HLT global trigger");

  // HLT payload also stored for not triggered events
  config.DefaultTriggerDomain().Add("*******", "HLT ");
  AliHLTReadoutList readoutlist;
  readoutlist.Enable(AliHLTReadoutList::kHLT);
  config.DefaultTriggerDomain().Add(readoutlist);
  
  
  TObject* menu = AliHLTGlobalTriggerConfig::Menu()->Clone();
  menu->Print();
  
  // /////////////////////////////////////////////////////////////////////////////////////////	
  // Write the trigger menu object to the CDB.
  AliCDBId id("HLT/ConfigHLT/HLTGlobalTrigger", firstRun, lastRun, version);
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("ALICE HLT Matthias.Richter@cern.ch");
  metaData->SetComment("HM-PHYSICS-V0001");
  storage->Put(menu, id, metaData);


  // /////////////////////////////////////////////////////////////////////////////////////////	
  // /////////////////////////////////////////////////////////////////////////////////////////	
  // /////////////////////////////////////////////////////////////////////////////////////////	
  //
  // component configurations
  gROOT->LoadMacro("$ALICE_ROOT/HLT/exa/makeComponentConfigurationObject.C");

  // /////////////////////////////////////////////////////////////////////////////////////////	
  // configuration of BarrelMultiplicityTrigger instances
  makeComponentConfigurationObject("HLT/ConfigHLT/BarrelMultiplicityTrigger", "-mintracks 10", cdbPath, firstRun, lastRun);
  makeComponentConfigurationObject("HLT/ConfigHLT/BarrelHighMultiplicity", "-mintracks 100"  , cdbPath, firstRun, lastRun);
  makeComponentConfigurationObject("HLT/ConfigHLT/BarrelPt_v01", "-mintracks 1 -minpt 0.5"   , cdbPath, firstRun, lastRun);
  makeComponentConfigurationObject("HLT/ConfigHLT/BarrelPt_v02", "-mintracks 1 -minpt 1.0"   , cdbPath, firstRun, lastRun);
  makeComponentConfigurationObject("HLT/ConfigHLT/BarrelPt_v03", "-mintracks 1 -minpt 5.0"   , cdbPath, firstRun, lastRun);
}
