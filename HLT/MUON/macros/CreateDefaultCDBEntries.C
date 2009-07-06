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

/* $Id: $ */

/**
 * \ingroup macros
 * \file CreateDefaultCDBEntries.C
 * \brief Macro for generating default CDB entries for dHLT.
 *
 * This macro is used to generate some default CDB entries for the dHLT.
 * It is an experts macro so make sure you know what you are doing.
 *
 * The simplest way to run this macro with defaults is to copy "rootlogon.C" from
 * $ALICE_ROOT/HLT/MUON/macros into your current working directory, then from
 * the shell command prompt run the following command:
 * \code
 *   > aliroot -b -q -l $ALICE_ROOT/HLT/MUON/macros/CreateDefaultCDBEntries.C+
 * \endcode
 *
 * \author Artur Szostak <artursz@iafrica.com>
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliHLTMUONConstants.h"
#include "TMap.h"
#include "TObjString.h"
#include "TString.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

/**
 * Generates default CDB entries in the given CDB storage (local by default)
 * \param cdbPath  The path to the local storage.
 */
void CreateDefaultCDBEntries(
		const char* cdbPath = "local://$ALICE_ROOT/OCDB",
		float zmiddle = -975.,
		float bfieldintegral = -3.,
		int dccut = 50,
		float roiParamAchamber7 = 0.020,
		float roiParamBchamber7 = 2.2,
		float roiParamAchamber8 = 0.023,
		float roiParamBchamber8 = 2.3,
		float roiParamAchamber9 = 0.049,
		float roiParamBchamber9 = 4.8,
		float roiParamAchamber10 = 0.045,
		float roiParamBchamber10 = 4.2,
		float lowptcut = 0.779,
		float highptcut = 1.755,
		float lowmasscut = 1.5,
		float highmasscut = 6.,
		float chamber7postion = -1276.5,
		float chamber8postion = -1307.5,
		float chamber9postion = -1406.6,
		float chamber10postion = -1437.6,
		float chamber11postion = -1603.5,
		float chamber13postion = -1703.5
	)
{
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
	
	Int_t version = 0;
	Int_t firstRun = 0;
	Int_t lastRun = AliCDBRunRange::Infinity();
	
	const char* path = NULL;
	AliCDBMetaData* metaData = NULL;
	TMap* params = NULL;
	AliCDBId id;
	
	// Create and store the configuration parameters for the trigger reconstructor.
	params = new TMap;
	params->SetOwner(kTRUE);
	char valstr[1024];
	sprintf(valstr, "%8.8f", zmiddle);
	params->Add(new TObjString("zmiddle"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", bfieldintegral);
	params->Add(new TObjString("bfieldintegral"), new TObjString(valstr));
	
	path = AliHLTMUONConstants::TriggerReconstructorCDBPath();
	id = AliCDBId(path, firstRun, lastRun, version);
	metaData = new AliCDBMetaData();
	metaData->SetResponsible("dimuon HLT");
	metaData->SetComment("Trigger reconstructor configuration parameters for dimuon HLT.");
	storage->Put(params, id, metaData);
	
	// Create and store the configuration parameters for the hit reconstructor.
	params = new TMap;
	params->SetOwner(kTRUE);
	sprintf(valstr, "%d", dccut);
	params->Add(new TObjString("dccut"), new TObjString(valstr));
	
	path = AliHLTMUONConstants::HitReconstructorCDBPath();
	id = AliCDBId(path, firstRun, lastRun, version);
	metaData = new AliCDBMetaData();
	metaData->SetResponsible("dimuon HLT");
	metaData->SetComment("Hit reconstructor DC cut parameter for dimuon HLT.");
	storage->Put(params, id, metaData);
	
	// Create and store the configuration parameters for the Manso tracker.
	params = new TMap;
	params->SetOwner(kTRUE);
	sprintf(valstr, "%8.8f", zmiddle);
	params->Add(new TObjString("zmiddle"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", bfieldintegral);
	params->Add(new TObjString("bfieldintegral"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", roiParamAchamber7);
	params->Add(new TObjString("roi_paramA_chamber7"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", roiParamBchamber7);
	params->Add(new TObjString("roi_paramB_chamber7"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", roiParamAchamber8);
	params->Add(new TObjString("roi_paramA_chamber8"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", roiParamBchamber8);
	params->Add(new TObjString("roi_paramB_chamber8"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", roiParamAchamber9);
	params->Add(new TObjString("roi_paramA_chamber9"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", roiParamBchamber9);
	params->Add(new TObjString("roi_paramB_chamber9"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", roiParamAchamber10);
	params->Add(new TObjString("roi_paramA_chamber10"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", roiParamBchamber10);
	params->Add(new TObjString("roi_paramB_chamber10"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", chamber7postion);
	params->Add(new TObjString("chamber7postion"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", chamber8postion);
	params->Add(new TObjString("chamber8postion"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", chamber9postion);
	params->Add(new TObjString("chamber9postion"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", chamber10postion);
	params->Add(new TObjString("chamber10postion"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", chamber11postion);
	params->Add(new TObjString("chamber11postion"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", chamber13postion);
	params->Add(new TObjString("chamber13postion"), new TObjString(valstr));
	
	path = AliHLTMUONConstants::MansoTrackerFSMCDBPath();
	id = AliCDBId(path, firstRun, lastRun, version);
	metaData = new AliCDBMetaData();
	metaData->SetResponsible("dimuon HLT");
	metaData->SetComment("Manso tracker FSM component configuration parameters for dimuon HLT.");
	storage->Put(params, id, metaData);
	
	// Create and store the configuration parameters for the trigger decision cuts.
	params = new TMap;
	params->SetOwner(kTRUE);
	sprintf(valstr, "%8.8f", lowptcut);
	params->Add(new TObjString("lowptcut"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", highptcut);
	params->Add(new TObjString("highptcut"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", lowmasscut);
	params->Add(new TObjString("lowmasscut"), new TObjString(valstr));
	sprintf(valstr, "%8.8f", highmasscut);
	params->Add(new TObjString("highmasscut"), new TObjString(valstr));
	
	path = AliHLTMUONConstants::DecisionComponentCDBPath();
	id = AliCDBId(path, firstRun, lastRun, version);
	metaData = new AliCDBMetaData();
	metaData->SetResponsible("dimuon HLT");
	metaData->SetComment("Trigger decision cuts for dimuon HLT.");
	storage->Put(params, id, metaData);
}
