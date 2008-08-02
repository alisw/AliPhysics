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
#include "AliMUONRecoParam.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDEStore.h"
#include "TMap.h"
#include "TObjString.h"
#include "TString.h"

ClassImp(AliHLTMUONProcessor)


int AliHLTMUONProcessor::SetCDBPathAndRunNo(
		const char* cdbPath, Int_t run, bool useDefault
	) const
{
	/// Sets the CDB path and run number to read from.
	/// \param cdbPath  The CDB path to use. If set to NULL and the path has
	///      not been set in the CDB manager then the default path
	///      "local://$ALICE_ROOT" is used if the 'useDefault' flag is also true.
	/// \param run  The run number to use. If set to -1 and the run number has
	///      not been set in the CDB manager then a value of zero is used if
	///      the 'useDefault' flag is also true.
	/// \param useDefault  If set to true then a default CDB path and/or run number
	///      is used if they have not been set and 'cdbPath' == NULL or
	///      'run' == -1.
	/// \return Zero if the object could be loaded. Otherwise an error code,
	///      compatible with the HLT framework, is returned.
	
	const char* defaultPath = "local://$ALICE_ROOT";
	Int_t defaultRun = 0;
	
	AliCDBManager* cdbManager = AliCDBManager::Instance();
	if (cdbManager == NULL)
	{
		HLTError("CDB manager instance does not exist.");
		return -EIO;
	}
	
	// Setup the CDB path.
	if (cdbPath != NULL)
	{
		cdbManager->SetDefaultStorage(cdbPath);
	}
	else if (not cdbManager->IsDefaultStorageSet() and useDefault)
	{
		cdbManager->SetDefaultStorage(defaultPath);
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
	
	return 0;
}


int AliHLTMUONProcessor::FetchMappingStores() const
{
	/// Fetches the DDL and detector element store objects for MUON mapping.
	/// \return Zero if the objects could be loaded. Otherwise an error code,
	///      which is compatible with the HLT framework, is returned.
	/// \note AliMpDDLStore::Instance() and AliMpDEStore::Instance() must be used
	///      to fetch the objects after this method returns a code equal to zero.
	
	Bool_t warn = kFALSE;
	
	// Check if the objects are already loaded. If they are then exit early,
	// otherwise we need to try load the objects.
	if (AliMpDDLStore::Instance(warn) != NULL and AliMpDEStore::Instance(warn) != NULL)
		return 0;
	
	AliCDBManager* cdbManager = AliCDBManager::Instance();
	if (cdbManager == NULL)
	{
		HLTError("CDB manager instance does not exist.");
		return -EIO;
	}
	
	const char* cdbPathUsed = "unknown (not set)";
	AliCDBStorage* store = cdbManager->GetDefaultStorage();
	if (store != NULL) cdbPathUsed = store->GetURI().Data();
	
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


int AliHLTMUONProcessor::FetchTMapFromCDB(const char* pathToEntry, TMap*& map) const
{
	/// Fetches a TMap object from the CDB.
	/// [in] \param pathToEntry  The relative path to the entry in the CDB to fetch.
	/// [out] \param map  This will be filled with the TMap object found if
	///      a successful status code is returned. Otherwise it will be unchanged.
	/// \return Zero if the object could be found. Otherwise an error code,
	///      which is compatible with the HLT framework, is returned.
	
	assert(AliCDBManager::Instance() != NULL);
	
	AliCDBStorage* store = AliCDBManager::Instance()->GetDefaultStorage();
	if (store == NULL)
	{
		HLTError("Could not get the the default storage for the CDB.");
		return -EIO;
	}

	Int_t version = store->GetLatestVersion(pathToEntry, GetRunNo());
	Int_t subVersion = store->GetLatestSubVersion(pathToEntry, GetRunNo(), version);
	AliCDBEntry* entry = AliCDBManager::Instance()->Get(pathToEntry, GetRunNo(), version, subVersion);
	if (entry == NULL)
	{
		HLTError("Could not get the CDB entry for \"%s\".", pathToEntry);
		return -EIO;
	}
	
	TObject* obj = entry->GetObject();
	if (obj == NULL)
	{
		HLTError("Configuration object for \"%s\" is missing.", pathToEntry);
		return -ENOENT;
	}
	
	if (obj->IsA() != TMap::Class())
	{
		HLTError("Wrong type for configuration object in \"%s\". Found a %s but we need a TMap.",
			pathToEntry, obj->ClassName()
		);
		return -EPROTO;
	}
	map = dynamic_cast<TMap*>(obj);
	
	return 0;
}


int AliHLTMUONProcessor::GetValueFromTMap(
		TMap* map, const char* paramName, TString& value,
		const char* pathToEntry, const char* prettyName
	) const
{
	/// Tries to find the string value associated with a certain parameter in a TMap.
	/// [in] \param map  The TMap object to search in.
	/// [in] \param paramName  The name of the parameter to search for.
	/// [out] \param value  Will be filled with the object found.
	/// [in] \param pathToEntry  The relative path to the entry in the CDB.
	///      Used when printing error messages. If set to NULL then a string of
	///      "(unknown)" is used. (default is NULL).
	/// [in] \param prettyName  Should be the name of the parameter which will
	///      be used when printing error messages. If this is set to NULL then
	///      the paramName will be used instead (default is NULL).
	/// \return Zero if the object could be found. Otherwise an error code,
	///      which is compatible with the HLT framework, is returned.
	
	if (pathToEntry == NULL) pathToEntry = "(unknown)";
	if (prettyName == NULL) prettyName = paramName;
	
	TPair* pair = static_cast<TPair*>(map->FindObject(paramName));
	if (pair == NULL)
	{
		HLTError("Configuration object for \"%s\" does not contain the %s value.",
			pathToEntry, prettyName
		);
		return -ENOENT;
	}
	TObject* valueObj = pair->Value();
	if (valueObj->IsA() != TObjString::Class())
	{
		HLTError("The %s parameter found in configuration object \"%s\""
			" is not a TObjString. Found an object of type %s instead.",
			prettyName, pathToEntry, valueObj->ClassName()
		);
		return -EPROTO;
	}
	value = dynamic_cast<TObjString*>(valueObj)->GetString();
	
	return 0;
}


int AliHLTMUONProcessor::GetIntFromTMap(
		TMap* map, const char* paramName, Int_t& value,
		const char* pathToEntry, const char* prettyName
	) const
{
	/// Tries to find a certain parameter in the TMap object and convert it to
	/// an integer value.
	/// [in] \param map  The TMap object to search in.
	/// [in] \param paramName  The name of the parameter to search for.
	/// [out] \param value  Will be filled with the integer value for the parameter,
	///       if it was found and it was an integer value.
	/// [in] \param pathToEntry  The relative path to the entry in the CDB.
	///      Used when printing error messages. If set to NULL then a string of
	///      "(unknown)" is used. (default is NULL).
	/// [in] \param prettyName  Should be the name of the parameter which will
	///      be used when printing error messages. If this is set to NULL then
	///      the paramName will be used instead (default is NULL).
	/// \return Zero if the object could be found and is valid. Otherwise an
	///       error code, which is compatible with the HLT framework, is returned.
	
	if (pathToEntry == NULL) pathToEntry = "(unknown)";
	if (prettyName == NULL) prettyName = paramName;
	
	TString valueStr;
	int result = GetValueFromTMap(map, paramName, valueStr, pathToEntry, prettyName);
	if (result != 0) return result;
	
	if (not valueStr.IsDigit())
	{
		HLTError("The %s parameter found in configuration object \"%s\""
			"is not a valid integer number string; found \"%s\".",
			prettyName, pathToEntry, valueStr.Data()
		);
		return -EPROTO;
	}
	value = valueStr.Atoi();
	
	return 0;
}


int AliHLTMUONProcessor::GetPositiveIntFromTMap(
		TMap* map, const char* paramName, Int_t& value,
		const char* pathToEntry, const char* prettyName
	) const
{
	/// Tries to find a certain parameter in the TMap object and convert it to
	/// a positive integer value.
	/// [in] \param map  The TMap object to search in.
	/// [in] \param paramName  The name of the parameter to search for.
	/// [out] \param value  Will be filled with the integer value for the parameter,
	///       if it was found and it was a positive integer value.
	/// [in] \param pathToEntry  The relative path to the entry in the CDB.
	///      Used when printing error messages. If set to NULL then a string of
	///      "(unknown)" is used. (default is NULL).
	/// [in] \param prettyName  Should be the name of the parameter which will
	///      be used when printing error messages. If this is set to NULL then
	///      the paramName will be used instead (default is NULL).
	/// \return Zero if the object could be found and is valid. Otherwise an
	///       error code, which is compatible with the HLT framework, is returned.
	
	if (pathToEntry == NULL) pathToEntry = "(unknown)";
	if (prettyName == NULL) prettyName = paramName;
	
	TString valueStr;
	int result = GetValueFromTMap(map, paramName, valueStr, pathToEntry, prettyName);
	if (result != 0) return result;
	
	if (not valueStr.IsDigit())
	{
		HLTError("The %s parameter found in configuration object \"%s\""
			"is not a valid integer number string; found \"%s\".",
			prettyName, pathToEntry, valueStr.Data()
		);
		return -EPROTO;
	}
	Int_t val = valueStr.Atoi();
	if (val < 0)
	{
		HLTError("The %s parameter found in configuration object \"%s\""
			"is not a positive integer number; found \"%d\".",
			prettyName, pathToEntry, val
		);
		return -EPROTO;
	}
	value = val;
	
	return 0;
}


int AliHLTMUONProcessor::GetFloatFromTMap(
		TMap* map, const char* paramName, Double_t& value,
		const char* pathToEntry, const char* prettyName
	) const
{
	/// Tries to find a certain parameter in the TMap object and convert it to
	/// an floating point value.
	/// [in] \param map  The TMap object to search in.
	/// [in] \param paramName  The name of the parameter to search for.
	/// [out] \param value  Will be filled with the floating point value for the
	///       parameter, if it was found and it was a floating point value.
	/// [in] \param pathToEntry  The relative path to the entry in the CDB.
	///      Used when printing error messages. If set to NULL then a string of
	///      "(unknown)" is used. (default is NULL).
	/// [in] \param prettyName  Should be the name of the parameter which will
	///      be used when printing error messages. If this is set to NULL then
	///      the paramName will be used instead (default is NULL).
	/// \return Zero if the object could be found and is valid. Otherwise an
	///       error code, which is compatible with the HLT framework, is returned.
	
	if (pathToEntry == NULL) pathToEntry = "(unknown)";
	if (prettyName == NULL) prettyName = paramName;
	
	TString valueStr;
	int result = GetValueFromTMap(map, paramName, valueStr, pathToEntry, prettyName);
	if (result != 0) return result;
	
	if (not valueStr.IsFloat())
	{
		HLTError("The %s parameter found in configuration object \"%s\""
			"is not a valid floating point number string; found \"%s\".",
			prettyName, pathToEntry, valueStr.Data()
		);
		return -EPROTO;
	}
	value = valueStr.Atof();
	
	return 0;
}


int AliHLTMUONProcessor::GetPositiveFloatFromTMap(
		TMap* map, const char* paramName, Double_t& value,
		const char* pathToEntry, const char* prettyName
	) const
{
	/// Tries to find a certain parameter in the TMap object and convert it to
	/// an positive floating point value.
	/// [in] \param map  The TMap object to search in.
	/// [in] \param paramName  The name of the parameter to search for.
	/// [out] \param value  Will be filled with the floating point value for the
	///       parameter, if it was found and it was a positive floating point value.
	/// [in] \param pathToEntry  The relative path to the entry in the CDB.
	///      Used when printing error messages. If set to NULL then a string of
	///      "(unknown)" is used. (default is NULL).
	/// [in] \param prettyName  Should be the name of the parameter which will
	///      be used when printing error messages. If this is set to NULL then
	///      the paramName will be used instead (default is NULL).
	/// \return Zero if the object could be found and is valid. Otherwise an
	///       error code, which is compatible with the HLT framework, is returned.
	
	if (pathToEntry == NULL) pathToEntry = "(unknown)";
	if (prettyName == NULL) prettyName = paramName;
	
	TString valueStr;
	int result = GetValueFromTMap(map, paramName, valueStr, pathToEntry, prettyName);
	if (result != 0) return result;
	
	if (not valueStr.IsFloat())
	{
		HLTError("The %s parameter found in configuration object \"%s\""
			"is not a valid floating point number string; found \"%s\".",
			prettyName, pathToEntry, valueStr.Data()
		);
		return -EPROTO;
	}
	Double_t val = valueStr.Atof();
	if (val < 0)
	{
		HLTError("The %s parameter found in configuration object \"%s\""
			"is not a positive floating point number; found \"%d\".",
			prettyName, pathToEntry, val
		);
		return -EPROTO;
	}
	value = val;
	
	return 0;
}


int AliHLTMUONProcessor::LoadRecoParamsFromCDB(AliMUONRecoParam*& params) const
{
	/// Fetches the reconstruction parameters object from the CDB for MUON.
	/// [out] \param params  This will be filled with the reconstruction
	///      parameters object found if a successful status code is returned.
	///      Otherwise it will be unchanged.
	/// \return Zero if the object could be found. Otherwise an error code,
	///      which is compatible with the HLT framework, is returned.
	
	assert(AliCDBManager::Instance() != NULL);

	AliCDBStorage* store = AliCDBManager::Instance()->GetDefaultStorage();
	if (store == NULL)
	{
		HLTError("Could not get the the default storage for the CDB.");
		return -EIO;
	}

	const char* pathToEntry = "MUON/Calib/RecoParam";
	Int_t version = store->GetLatestVersion(pathToEntry, GetRunNo());
	Int_t subVersion = store->GetLatestSubVersion(pathToEntry, GetRunNo(), version);
	AliCDBEntry* entry = AliCDBManager::Instance()->Get(pathToEntry, GetRunNo(), version, subVersion);
	if (entry == NULL)
	{
		HLTError("Could not get the CDB entry for \"%s\".", pathToEntry);
		return -EIO;
	}
	
	TObject* obj = entry->GetObject();
	if (obj == NULL)
	{
		HLTError("Reconstruction parameters object for \"%s\" is missing.", pathToEntry);
		return -ENOENT;
	}
	
	TObjArray* objarr = dynamic_cast<TObjArray*>(obj);
	if (objarr != NULL)
	{
		obj = objarr->Last();
	}
	
	AliMUONRecoParam* par = dynamic_cast<AliMUONRecoParam*>(obj);
	if (par == NULL)
	{
		HLTError("No AliMUONRecoParam class found for entry \"%s\". Found a %s class instead.",
			pathToEntry, obj->ClassName()
		);
		return -EPROTO;
	}
	
	params = par;
	return 0;
}

