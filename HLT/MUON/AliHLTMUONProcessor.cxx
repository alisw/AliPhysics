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

///
/// @file   AliHLTMUONProcessor.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 May 2008
/// @brief  Implementation of the abstract base dHLT processor component.
///
/// This component is the abstract base class of dHLT specific components.
/// It implements some common methods used by all the dHLT components.
///

#include "AliHLTMUONProcessor.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDEStore.h"

ClassImp(AliHLTMUONProcessor)


int AliHLTMUONProcessor::FetchMappingStores(
		const char* cdbPath, Int_t run, bool useDefault
	) const
{
	/// Fetches the DDL and detector element store objects for MUON mapping.
	/// \param cdbPath  The CDB path to use. If set to NULL and the path has
	///      not been set in the CDB manager then the default path
	///      "local://$ALICE_ROOT" is used if the 'useDefault' flag is also true.
	/// \param run  The run number to use. If set to -1 and the run number has
	///      not been set in the CDB manager then a value of zero is used if
	///      the 'useDefault' flag is also true.
	/// \param useDefault  If set to true then a default CDB path and run number
	///      is used if they have not been set and 'cdbPath' == NULL or
	///      'run' == NULL.
	/// \return Zero if the object could be loaded. Otherwise an error code is
	///      returned, which is compatible with the HLT framework.
	/// \note AliMpDDLStore::Instance() and AliMpDEStore::Instance() must be used
	///      to fetch the objects after this method returns a code equal to zero.
	
	Bool_t warn = kFALSE;
	
	// Check if the objects are already loaded. If they are then exit early,
	// otherwise we need to try load the objects.
	if (AliMpDDLStore::Instance(warn) != NULL and AliMpDEStore::Instance(warn) != NULL)
		return 0;
	
	const char* defaultPath = "local://$ALICE_ROOT";
	Int_t defaultRun = 0;
	
	AliCDBManager* cdbManager = AliCDBManager::Instance();
	if (cdbManager == NULL)
	{
		HLTError("CDB manager instance does not exist.");
		return -EIO;
	}
	
	// Setup the CDB path.
	const char* cdbPathUsed = "unknown (not set)";
	if (cdbPath != NULL)
	{
		cdbManager->SetDefaultStorage(cdbPath);
		cdbPathUsed = cdbPath;
	}
	else
	{
		AliCDBStorage* store = cdbManager->GetDefaultStorage();
		if (store == NULL)
		{
			if (useDefault)
			{
				cdbManager->SetDefaultStorage(defaultPath);
				cdbPathUsed = defaultPath;
			}
		}
		else
		{
			cdbPathUsed = store->GetURI().Data();
		}
	}
	
	// Now setup the run number.
	if (run != -1)
	{
		cdbManager->SetRun(run);
	}
	else
	{
		if (useDefault) cdbManager->SetRun(defaultRun);
	}
	Int_t runUsed = cdbManager->GetRun();
	
	// Now we can try load the DDL and DE store objects.
	if (not AliMpCDB::LoadDDLStore(warn))
	{
		HLTError("Failed to load DDL or detector element store specified"
			 " for CDB path '%s' and run no. %d.",
			cdbPathUsed, runUsed
		);
		return -ENOENT;
	}
	
	if (AliMpDDLStore::Instance(warn) == NULL or AliMpDEStore::Instance(warn) == NULL)
	{
		HLTError("Could not find or load the DDL or detector element store instance.");
		return -EIO;
	}
	
	return 0;
}

