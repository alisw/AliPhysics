//$Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**
 * \ingroup macros
 * \file CreateTriggerMenuCDBEntry.C
 * \brief Macro for generating the trigger menu CDB entry for the HLT global trigger.
 *
 * This macro is used to generate a default CDB entry for the HLT global trigger.
 * The trigger menu is used by the HLT global trigger as a configuration object.
 * It is an experts macro so make sure you know what you are doing.
 *
 * You can run this macro with defaults using the following shell command:
 * \code
 *   > aliroot -b -q $ALICE_ROOT/HLT/trigger/CreateTriggerMenuCDBEntry.C
 * \endcode
 *
 * \author Artur Szostak <artursz@iafrica.com>
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
void CreateTriggerMenuCDBEntry(
		const char* cdbPath = "local://$ALICE_ROOT/OCDB",
		Int_t version = 0,
		Int_t firstRun = 0,
		Int_t lastRun = AliCDBRunRange::Infinity()
	)
{
	gSystem->Load("libAliHLTUtil.so");
	gSystem->Load("libAliHLTTRD.so");
	gSystem->Load("libAliHLTMUON.so");
	gSystem->Load("libAliHLTTrigger.so");

	// Setup the CDB default storage and run number.
	AliCDBManager* cdbManager = AliCDBManager::Instance();
	if (cdbManager == NULL)
	{
		cerr << "ERROR: Global CDB manager object does not exist." << endl;
		return;
	}
	AliCDBStorage* storage = cdbManager->GetStorage(cdbPath);
	if (storage == NULL)
	{
		cerr << "ERROR: Could not get storage for: " << cdbPath << endl;
		return;
	}

	///////////////////////////////////////////////////////////////////////////////////////////	
	// Create the trigger menu.
	AliHLTGlobalTriggerConfig config("Default Global Trigger Config");

	config.AddSymbol("domainAll", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:***\")");

	///////////////////////////////////////////////////////////////////////////////////////////	
	// the domain definitions for the global HLT output and the HLT DDLs
	config.AddSymbol("domainHLTOUT", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"*******:HLT \")");
	config.AddSymbol("domainHLTDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:HLT\\0\")");
	config.AddSymbol("domainALLDDL", "AliHLTTriggerDomain", "", "AliHLTTriggerDomain(\"DAQRDOUT:***\\0\")");

	///////////////////////////////////////////////////////////////////////////////////////////	
	// NOTE: always make sure that the global HLT output and the HLT DDLs are included
	// in the readout, i.e. add domainHLTOUT|domainHLTDDL to the trigger domain
	config.AddItem("BarrelMultiplicityTrigger", "BarrelMultiplicityTrigger|domainHLTOUT|domainALLDDL", "charged barrel track multiplicity triggered");
	config.AddItem("MuonSpectroTrigger", "MuonSpectroTrigger|domainHLTOUT|domainALLDDL", "Muon spectrometer triggered");
	
	///////////////////////////////////////////////////////////////////////////////////////////	
	// default domain in case there is no global trigger
	// readout the output of the reconstruction
	// this refers to the domain domainHLTOUT|domainHLTDDL
	config.SetDefaultTriggerDescription("No HLT global trigger");
	config.DefaultTriggerDomain().Add("*******", "HLT ");
	AliHLTReadoutList readoutlist;
	readoutlist.Enable(AliHLTReadoutList::kHLT);
	config.DefaultTriggerDomain().Add(readoutlist);

	TObject* menu = AliHLTGlobalTriggerConfig::Menu()->Clone();
	menu->Print();

	///////////////////////////////////////////////////////////////////////////////////////////	
	// Write the trigger menu object to the CDB.
	AliCDBId id("HLT/ConfigHLT/HLTGlobalTrigger", firstRun, lastRun, version);
	AliCDBMetaData* metaData = new AliCDBMetaData();
	metaData->SetResponsible("HLT");
	metaData->SetComment("Default trigger menu for HLT global trigger.");
	storage->Put(menu, id, metaData);
}
