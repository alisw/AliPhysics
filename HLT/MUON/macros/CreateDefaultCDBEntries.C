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
void CreateDefaultCDBEntries(const char* cdbPath = "local://$ALICE_ROOT")
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
	
	Int_t verison = 0;
	Int_t firstRun = 0;
	Int_t lastRun = AliCDBRunRange::Infinity();
	
	const char* path = NULL;
	AliCDBMetaData* metaData = NULL;
	TMap* params = NULL;
	AliCDBId id;
	
	// Create and store the configuration parameters for the hit reconstructor.
	params = new TMap;
	params->SetOwner(kTRUE);
	params->Add(new TObjString("dccut"), new TObjString("50"));
	
	path = AliHLTMUONConstants::HitReconstructorCDBPath();
	id = AliCDBId(path, firstRun, lastRun, verison);
	metaData = new AliCDBMetaData();
	metaData->SetResponsible("dimuon HLT");
	metaData->SetComment("Hit reconstructor DC cut parameter for dimuon HLT.");
	storage->Put(params, id, metaData);
	
	// Create and store the configuration parameters for the trigger decision cuts.
	params = new TMap;
	params->SetOwner(kTRUE);
	params->Add(new TObjString("lowptcut"), new TObjString("1"));
	params->Add(new TObjString("highptcut"), new TObjString("2"));
	params->Add(new TObjString("lowmasscut"), new TObjString("2.5"));
	params->Add(new TObjString("highmasscut"), new TObjString("7"));
	
	path = AliHLTMUONConstants::DecisionComponentCDBPath();
	id = AliCDBId(path, firstRun, lastRun, verison);
	metaData = new AliCDBMetaData();
	metaData->SetResponsible("dimuon HLT");
	metaData->SetComment("Trigger decision cuts for dimuon HLT.");
	storage->Put(params, id, metaData);
}
