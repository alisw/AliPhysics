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

// $Id: $

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
#include "AliHLTMUONConstants.h"
#include "AliHLTMisc.h"
#include "AliMUONRecoParam.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliMagF.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDEStore.h"
#include "TGeoGlobalMagField.h"
#include "TMap.h"
#include "TObjString.h"
#include "TString.h"
#include <string>
#include <cstdlib>
#include <fstream>


ClassImp(AliHLTMUONProcessor)


AliHLTMUONProcessor::AliHLTMUONProcessor() :
	AliHLTProcessor(),
	fWarnForUnexpecedBlock(false),
	fDelaySetup(false),
	fDumpDataOnError(false),
	fDumpPath("./")
{
	/// Default constructor.
}


int AliHLTMUONProcessor::DoInit(int argc, const char** argv)
{
	/// Parses common dHLT component arguments.

	// Set the default values for various arguments comming from the command line.
	fDelaySetup = false;
	fDumpDataOnError = false;
	fDumpPath = "./";
	const char* cdbPath = NULL;
	Int_t run = -1;

	for (int i = 0; i < argc; i++)
	{
		// Ignore the argument if the child class indicates to do so.
		if (IgnoreArgument(argv[i])) continue;
	
		if (strcmp(argv[i], "-cdbpath") == 0)
		{
			if (cdbPath != NULL)
			{
				HLTWarning("CDB path was already specified. Will"
					" replace previous value given by -cdbpath."
				);
			}
			if (argc <= i+1)
			{
				HLTError("The CDB path was not specified." );
				return -EINVAL;
			}
			cdbPath = argv[i+1];
			i++;
			continue;
		}
	
		if (strcmp(argv[i], "-run") == 0)
		{
			if (run != -1)
			{
				HLTWarning("Run number was already specified. Will"
					" replace previous value given by -run."
				);
			}
			if (argc <= i+1)
			{
				HLTError("The run number was not specified.");
				return -EINVAL;
			}
			
			char* cpErr = NULL;
			run = Int_t( strtol(argv[i+1], &cpErr, 0) );
			if (cpErr == NULL or *cpErr != '\0' or run < 0)
			{
				HLTError("Cannot convert '%s' to a valid run number."
					" Expected a positive integer value.", argv[i+1]
				);
				return -EINVAL;
			}
			
			i++;
			continue;
		}
		
		if (strcmp(argv[i], "-delaysetup") == 0)
		{
			fDelaySetup = true;
			continue;
		}
		
		if (strcmp(argv[i], "-dumponerror") == 0)
		{
			fDumpDataOnError = true;
			continue;
		}
		
		if (strcmp(argv[i], "-dumppath") == 0)
		{
			if (fDumpPath != NULL)
			{
				HLTWarning("The dump path was already specified. Will"
					" replace previous value given by -dumppath."
				);
			}
			if (argc <= i+1)
			{
				HLTError("The dump path was not specified.");
				return -EINVAL;
			}
			fDumpPath = argv[i+1];
			i++;
			continue;
		}
	}
	
	if (cdbPath != NULL or run != -1)
	{
		int result = SetCDBPathAndRunNo(cdbPath, run);
		if (result != 0)
		{
			// Error messages already generated in SetCDBPathAndRunNo.
			return result;
		}
	}

	return 0;
}


bool AliHLTMUONProcessor::ArgumentAlreadyHandled(int& i, const char* argi) const
{
	/// This method can be used by the derivind child class to check if a particular
	/// argument in argv was already processed.

	if (strcmp(argi, "-cdbpath") == 0)
	{
		if (IgnoreArgument(argi)) return false;
		i++;  // argument takes one parameter
		return true;
	}

	if (strcmp(argi, "-run") == 0)
	{
		if (IgnoreArgument(argi)) return false;
		i++;  // argument takes one parameter
		return true;
	}
	
	if (strcmp(argi, "-delaysetup") == 0)
	{
		if (IgnoreArgument(argi)) return false;
		return true;
	}
	
	if (strcmp(argi, "-dumponerror") == 0)
	{
		if (IgnoreArgument(argi)) return false;
		return true;
	}
	
	if (strcmp(argi, "-dumppath") == 0)
	{
		if (IgnoreArgument(argi)) return false;
		i++;  // argument takes one parameter
		return true;
	}

	return false;
}


int AliHLTMUONProcessor::SetCDBPathAndRunNo(
		const char* cdbPath, Int_t run, bool useDefault
	) const
{
	/// Sets the CDB path and run number to read from.
	/// \param cdbPath  The CDB path to use. If set to NULL and the path has
	///      not been set in the CDB manager then the default path
	///      "local://$ALICE_ROOT/OCDB" is used if the 'useDefault' flag is also true.
	/// \param run  The run number to use. If set to -1 and the run number has
	///      not been set in the CDB manager then a value of zero is used if
	///      the 'useDefault' flag is also true.
	/// \param useDefault  If set to true then a default CDB path and/or run number
	///      is used if they have not been set and 'cdbPath' == NULL or
	///      'run' == -1.
	/// \return Zero if the object could be loaded. Otherwise an error code,
	///      compatible with the HLT framework, is returned.
	
	const char* defaultPath = "local://$ALICE_ROOT/OCDB";
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
	if (AliHLTMisc::Instance().LoadOCDBEntry("MUON/Calib/MappingData", runUsed) == NULL)
	{
		HLTError("Could not find entry in CDB path '%s/MUON/Calib/MappingData' and run no. %d.",
			cdbPathUsed, runUsed
		);
		return -ENOENT;
	}
	if (AliHLTMisc::Instance().LoadOCDBEntry("MUON/Calib/Pedestals", runUsed) == NULL)
	{
		HLTError("Could not find entry in CDB path '%s/MUON/Calib/Pedestals' and run no. %d.",
			cdbPathUsed, runUsed
		);
		return -ENOENT;
	}
	if (not AliMpCDB::LoadDDLStore(warn))
	{
		HLTError("Failed to load DDL or detector element store specified"
			 " for CDB path '%s' and run no. %d.",
			cdbPathUsed, runUsed
		);
		return -EIO;
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
	/// \param [in] pathToEntry  The relative path to the entry in the CDB to fetch.
	/// \param [out] map  This will be filled with the TMap object found if
	///      a successful status code is returned. Otherwise it will be unchanged.
	/// \return Zero if the object could be found. Otherwise an error code,
	///      which is compatible with the HLT framework, is returned.
	
	TObject* obj = LoadAndExtractOCDBObject(pathToEntry);
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
	map = static_cast<TMap*>(obj);
	
	return 0;
}


int AliHLTMUONProcessor::GetValueFromTMap(
		TMap* map, const char* paramName, TString& value,
		const char* pathToEntry, const char* prettyName
	) const
{
	/// Tries to find the string value associated with a certain parameter in a TMap.
	/// \param [in] map  The TMap object to search in.
	/// \param [in] paramName  The name of the parameter to search for.
	/// \param [out] value  Will be filled with the object found.
	/// \param [in] pathToEntry  The relative path to the entry in the CDB.
	///      Used when printing error messages. If set to NULL then a string of
	///      "(unknown)" is used. (default is NULL).
	/// \param [in] prettyName  Should be the name of the parameter which will
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
	value = static_cast<TObjString*>(valueObj)->GetString();
	
	return 0;
}


int AliHLTMUONProcessor::GetIntFromTMap(
		TMap* map, const char* paramName, Int_t& value,
		const char* pathToEntry, const char* prettyName
	) const
{
	/// Tries to find a certain parameter in the TMap object and convert it to
	/// an integer value.
	/// \param [in] map  The TMap object to search in.
	/// \param [in] paramName  The name of the parameter to search for.
	/// \param [out] value  Will be filled with the integer value for the parameter,
	///       if it was found and it was an integer value.
	/// \param [in] pathToEntry  The relative path to the entry in the CDB.
	///      Used when printing error messages. If set to NULL then a string of
	///      "(unknown)" is used. (default is NULL).
	/// \param [in] prettyName  Should be the name of the parameter which will
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
	/// \param [in] map  The TMap object to search in.
	/// \param [in] paramName  The name of the parameter to search for.
	/// \param [out] value  Will be filled with the integer value for the parameter,
	///       if it was found and it was a positive integer value.
	/// \param [in] pathToEntry  The relative path to the entry in the CDB.
	///      Used when printing error messages. If set to NULL then a string of
	///      "(unknown)" is used. (default is NULL).
	/// \param [in] prettyName  Should be the name of the parameter which will
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
	/// \param [in] map  The TMap object to search in.
	/// \param [in] paramName  The name of the parameter to search for.
	/// \param [out] value  Will be filled with the floating point value for the
	///       parameter, if it was found and it was a floating point value.
	/// \param [in] pathToEntry  The relative path to the entry in the CDB.
	///      Used when printing error messages. If set to NULL then a string of
	///      "(unknown)" is used. (default is NULL).
	/// \param [in] prettyName  Should be the name of the parameter which will
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
	/// \param [in] map  The TMap object to search in.
	/// \param [in] paramName  The name of the parameter to search for.
	/// \param [out] value  Will be filled with the floating point value for the
	///       parameter, if it was found and it was a positive floating point value.
	/// \param [in] pathToEntry  The relative path to the entry in the CDB.
	///      Used when printing error messages. If set to NULL then a string of
	///      "(unknown)" is used. (default is NULL).
	/// \param [in] prettyName  Should be the name of the parameter which will
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


int AliHLTMUONProcessor::FetchFieldIntegral(Double_t& bfieldintegral) const
{
	// Fetches the correct dipole magnetic field integral to use.
	
	Float_t currentL3 = 0;
	Float_t currentDip = 0;
	
	if (TGeoGlobalMagField::Instance() == NULL or
	    (TGeoGlobalMagField::Instance() != NULL and TGeoGlobalMagField::Instance()->GetField() == NULL)
	   )
	{
		HLTWarning("The magnetic field has not been set in TGeoGlobalMagField."
			" Will try and load the GRP entry directly."
		);
		
		AliGRPManager grpman;
		if (not grpman.ReadGRPEntry() or grpman.GetGRPData() == NULL)
		{
			HLTError("GRP entry could not be loaded.");
			return -EIO;
		}
		
		const AliGRPObject* grp = grpman.GetGRPData();
		Char_t polarityL3 = grp->GetL3Polarity();
		Char_t polarityDip = grp->GetDipolePolarity();
		currentL3 = grp->GetL3Current(AliGRPObject::kMean);
		currentDip = grp->GetDipoleCurrent(AliGRPObject::kMean);
		if (polarityL3 == AliGRPObject::GetInvalidChar())
		{
			HLTError("L3 polarity in GRP is invalid.");
			return -EPROTO;
		}
		if (polarityDip == AliGRPObject::GetInvalidChar())
		{
			HLTError("Dipole polarity in GRP is invalid.");
			return -EPROTO;
		}
		if (currentL3 == AliGRPObject::GetInvalidFloat())
		{
			HLTError("L3 current in GRP is invalid.");
			return -EPROTO;
		}
		if (currentDip == AliGRPObject::GetInvalidFloat())
		{
			HLTError("Dipole current in GRP is invalid.");
			return -EPROTO;
		}
		if (grp->IsPolarityConventionLHC())
		{
			currentL3 *= (polarityL3 ? -1:1);
			currentDip *= (polarityDip ? -1:1);
		}
		else
		{
			currentL3 *= (polarityL3 ? -1:1);
			currentDip *= (polarityDip ? 1:-1);
		}
	}
	else
	{
		TVirtualMagField* vfield = TGeoGlobalMagField::Instance()->GetField();
		AliMagF* field = dynamic_cast<AliMagF*>(vfield);
		if (vfield->IsA() != AliMagF::Class() or field == NULL)
		{
			HLTError(Form(
				"The magnetic field is not of type AliMagF."
				" Do not know how to handle class of type '%s'.",
				vfield->ClassName()
			));
			return -EPROTO;
		}
		currentL3 = field->GetCurrentSol();
		currentDip = field->GetCurrentDip();
	}
	
	const char* path = AliHLTMUONConstants::FieldIntegralsCDBPath();
	TMap* map = NULL;
	int result = FetchTMapFromCDB(path, map);
	if (result != 0) return result;
	const char* paramName = Form("L3_current=%0.2e;Dipole_current=%0.2e", currentL3, currentDip);
	Double_t value;
	result = GetFloatFromTMap(map, paramName, value, path);
	if (result != 0) return result;
	bfieldintegral = value;
	return 0;
}


int AliHLTMUONProcessor::LoadRecoParamsFromCDB(AliMUONRecoParam*& params) const
{
	/// Fetches the reconstruction parameters object from the CDB for MUON.
	/// \param [out] params  This will be filled with the reconstruction
	///      parameters object found if a successful status code is returned.
	///      Otherwise it will be unchanged.
	/// \return Zero if the object could be found. Otherwise an error code,
	///      which is compatible with the HLT framework, is returned.
	
	const char* pathToEntry = "MUON/Calib/RecoParam";
	TObject* obj = LoadAndExtractOCDBObject(pathToEntry);
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


void AliHLTMUONProcessor::DumpBuffer(
		const void* buffer, AliHLTUInt32_t size, const char* filename
	) const
{
	/// Dumps the data contained in a buffer to file as is.

	using std::fstream;

	fstream file(filename, fstream::out | fstream::trunc | fstream::binary);
	if (file.good())
	{
		file.write(reinterpret_cast<const char*>(buffer), size);
		if (file.fail())
		{
			HLTError("Could not write data block to file %s during"
				" dumping operation!",
				filename
			);
		}
	}
	else
	{
		HLTError("Could not open file %s for dumping data block!", filename);
	}
}


void AliHLTMUONProcessor::DumpBlock(
		const AliHLTComponentBlockData* block, const char* fileNamePrefix
	) const
{
	/// Dumps the data block and meta information to file.

	std::string filename = fDumpPath;
	filename += fileNamePrefix;
	filename += "-blockmeta.bin";
	DumpBuffer(block, sizeof(AliHLTComponentBlockData), filename.c_str());
	filename = fDumpPath;
	filename += fileNamePrefix;
	filename += "-data.bin";
	DumpBuffer(block->fPtr, block->fSize, filename.c_str());
}


void AliHLTMUONProcessor::DumpEvent(
		const AliHLTComponentEventData& evtData,
		const AliHLTComponentBlockData* blocks,
		AliHLTComponentTriggerData& trigData,
		AliHLTUInt8_t* outputPtr,
		AliHLTUInt32_t& size,
		AliHLTComponentBlockDataList& outputBlocks
	) const
{
	/// Dumps the event information to files in the dump path given by the
	/// method DumpPath, which can be set by the command line argument -dumppath.

	using std::fstream;
	char strbuf[1024];

	std::string filename = fDumpPath;
	sprintf(strbuf, "dump_event-0x%16.16llX.log", evtData.fEventID);
	filename += strbuf;
	fstream logfile(filename.c_str(), fstream::out | fstream::trunc);
	if (logfile.fail())
	{
		HLTError("Could not open log file '%s' for dump information.", filename.c_str());
		return;
	}

	filename = fDumpPath;
	sprintf(strbuf, "dump_event-0x%16.16llX-eventdata.bin", evtData.fEventID);
	filename += strbuf;
	logfile << "Dumping event data structure to file: " << filename << std::endl;
	DumpBuffer(&evtData, sizeof(AliHLTComponentEventData), filename.c_str());

	filename = fDumpPath;
	sprintf(strbuf, "dump_event-0x%16.16llX-triggerdata.bin", evtData.fEventID);
	filename += strbuf;
	logfile << "Dumping trigger data structure to file: " << filename << std::endl;
	DumpBuffer(&trigData, sizeof(AliHLTComponentTriggerData), filename.c_str());

	for (unsigned int n = 0; n < evtData.fBlockCnt; n++)
	{
		sprintf(strbuf, "dump_event-0x%16.16llX-block-0x%8.8X", evtData.fEventID, n);
		filename = strbuf;
		sprintf(strbuf, "0x%8.8X", blocks[n].fSpecification);
		logfile << "Found block with data type = " << DataType2Text(blocks[n].fDataType)
			<< ", specification = " << strbuf << ". Dumping to file: "
			<< filename << "-data.bin" << std::endl;
		DumpBlock(&blocks[n], filename.c_str());
	}

	filename = fDumpPath;
	sprintf(strbuf, "dump_event-0x%16.16llX-output-buffer.bin", evtData.fEventID);
	filename += strbuf;
	logfile << "Dumping output buffer to file: " << filename << std::endl;
	DumpBuffer(outputPtr, size, filename.c_str());

	for (size_t i = 0; i < outputBlocks.size(); i++)
	{
		sprintf(strbuf, "dump_event-0x%16.16llX-output-block-0x%8.8X", evtData.fEventID, int(i));
		filename = strbuf;
		sprintf(strbuf, "0x%8.8X", outputBlocks[i].fSpecification);
		logfile << "Generated output data block with type = "
			<< DataType2Text(outputBlocks[i].fDataType)
			<< ", specification = " << strbuf << ". Dumping to file: "
			<< filename << "-data.bin" << std::endl;
		DumpBlock(&outputBlocks[i], filename.c_str());
	}
}


void AliHLTMUONProcessor::DumpEvent(
		const AliHLTComponentEventData& evtData,
		AliHLTComponentTriggerData& trigData
	) const
{
	/// Dumps the event information to files in the dump path given by the
	/// method DumpPath, which can be set by the command line argument -dumppath.

	using std::fstream;
	char strbuf[1024];

	std::string filename = fDumpPath;
	sprintf(strbuf, "dump_event-0x%16.16llX.log", evtData.fEventID);
	filename += strbuf;
	fstream logfile(filename.c_str(), fstream::out | fstream::trunc);
	if (logfile.fail())
	{
		HLTError("Could not open log file '%s' for dump information.", filename.c_str());
		return;
	}

	filename = fDumpPath;
	sprintf(strbuf, "dump_event-0x%16.16llX-eventdata.bin", evtData.fEventID);
	filename += strbuf;
	logfile << "Dumping event data structure to file: " << filename << std::endl;
	DumpBuffer(&evtData, sizeof(AliHLTComponentEventData), filename.c_str());

	filename = fDumpPath;
	sprintf(strbuf, "dump_event-0x%16.16llX-triggerdata.bin", evtData.fEventID);
	filename += strbuf;
	logfile << "Dumping trigger data structure to file: " << filename << std::endl;
	DumpBuffer(&trigData, sizeof(AliHLTComponentTriggerData), filename.c_str());

	for (int i = 0; i < GetNumberOfInputBlocks(); i++)
	{
		const AliHLTComponentBlockData* block = GetInputBlock(i);
		sprintf(strbuf, "dump_event-0x%16.16llX-block-0x%8.8X", evtData.fEventID, i);
		filename = strbuf;
		sprintf(strbuf, "0x%8.8X", block->fSpecification);
		logfile << "Found block with data type = " << DataType2Text(block->fDataType)
			<< ", specification = " << strbuf << ". Dumping to file: "
			<< filename << "-data.bin" << std::endl;
		DumpBlock(block, filename.c_str());
	}
}

