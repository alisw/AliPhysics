/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.17  2007/08/28 16:03:30  acolla
Restored to v1.14:


Function Bool_t GetHLTStatus() added to preprocessor interface. It will return
the status of HLT read from the run logbook.

Revision 1.16  2007/08/22 09:20:50  hristov
Updated QA classes (Yves)

Revision 1.14  2007/05/30 06:35:21  jgrosseo
Adding functionality to the Shuttle/TestShuttle:
o) Function to retrieve list of sources from a given system (GetFileSources with id=0)
o) Function to retrieve list of IDs for a given source      (GetFileIDs)
These functions are needed for dealing with the tag files that are saved for the GRP preprocessor
Example code has been added to the TestProcessor in TestShuttle

Revision 1.13  2007/04/12 08:26:33  jgrosseo
updated commment

Revision 1.12  2007/04/05 08:05:55  acolla
Conversion from online to offline detector name in StoreReferenceFile

Revision 1.11  2007/04/04 10:29:18  jgrosseo
1) Storing of files to the Grid is now done _after_ your preprocessors succeeded. This is transparent, which means that you can still use the same functions (Store, StoreReferenceData) to store files to the Grid. However, the Shuttle first stores them locally and transfers them after the preprocessor finished. The return code of these two functions has changed from UInt_t to Bool_t which gives you the success of the storing.
In case of an error with the Grid, the Shuttle will retry the storing later, the preprocessor does not need to be run again.

2) The meaning of the return code of the preprocessor has changed. 0 is now success and any other value means failure. This value is stored in the log and you can use it to keep details about the error condition.

3) New function StoreReferenceFile to _directly_ store a file (without opening it) to the reference storage.

4) The memory usage of the preprocessor is monitored. If it exceeds 2 GB it is terminated.

5) New function AliPreprocessor::ProcessDCS(). If you do not need to have DCS data in all cases, you can skip the processing by implemting this function and returning kFALSE under certain conditions. E.g. if there is a certain run type.
If you always need DCS data (like before), you do not need to implement it.

6) The run type has been added to the monitoring page

Revision 1.9  2007/02/28 10:42:58  acolla
Run type field added in SHUTTLE framework. Run type is read from "run type" logbook and retrieved by
AliPreprocessor::GetRunType() function.

Revision 1.7  2006/11/06 14:24:21  jgrosseo
reading of run parameters from the logbook
online offline naming conversion

Revision 1.6  2006/10/02 12:57:48  jgrosseo
Small interface change of function StoreReferenceData in Shuttle

Revision 1.5  2006/09/04 17:42:34  hristov
Changes required by Effective C++

Revision 1.4  2006/08/08 14:20:49  jgrosseo
Update to shuttle classes (Alberto)

- Possibility to set the full object's path in the Preprocessor's and
Shuttle's  Store functions
- Possibility to extend the object's run validity in the same classes
("startValidity" and "validityInfinite" parameters)
- Implementation of the StoreReferenceData function to store reference
data in a dedicated CDB storage.

Revision 1.3  2006/07/11 12:42:43  jgrosseo
adding parameters for extended validity range of data produced by preprocessor

Revision 1.2  2006/06/06 16:36:49  jgrosseo
minor changes in AliShuttleInterface and AliPreprocessor

Revision 1.1  2006/06/02 14:14:36  hristov
Separate library for CDB (Jan)

Revision 1.2  2006/03/07 07:52:34  hristov
New version (B.Yordanov)

Revision 1.3  2005/11/17 17:47:34  byordano
TList changed to TObjArray

Revision 1.2  2005/11/17 14:43:22  byordano
import to local CVS

Revision 1.1.1.1  2005/10/28 07:33:58  hristov
Initial import as subdirectory in AliRoot

Revision 1.1.1.1  2005/09/12 22:11:40  byordano
SHUTTLE package

Revision 1.2  2005/08/29 21:15:47  byordano
some docs added

*/

// Description:
// This class is the CDBPreProcessor interface,
// supposed to be implemented by any detector
// interested in immediate processing of data 
// which is retrieved from DCS.
// For every particular run set of aliases and
// their corespoding value sets are returned.
// Usage schema:
//	1) virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime)	
//	This method is called at the begining of data retrieval.
//	run: run number
//	startTime: when the run started
//	endTime: when the run finished	
//
//	2) virtual void Process()
//
//	This method is called and passed a list of retrieved values from DCS
//
//


#include "AliPreprocessor.h"

#include <TString.h>
#include <TMap.h>
#include <TObjString.h>

#include "AliLog.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliShuttleInterface.h"

ClassImp(AliPreprocessor)

//______________________________________________________________________________________________
AliPreprocessor::AliPreprocessor(const char* detector, AliShuttleInterface* shuttle) :
  TNamed(detector, ""),
  fRun(-1),
  fStartTime(0),
  fEndTime(0),
  fShuttle(shuttle),
  fRunTypes()
{
	SetTitle(Form("AliPreprocessor for %s subdetector.", detector));

  if (!fShuttle)
  {
    AliFatal("Initialized without Shuttle instance.");
    return;
  }

  fShuttle->RegisterPreprocessor(this);
  
  fRunTypes.SetOwner(kTRUE);
}

//______________________________________________________________________________________________
AliPreprocessor::~AliPreprocessor()
{
}

//______________________________________________________________________________________________
void AliPreprocessor::Initialize(Int_t run, UInt_t startTime,	UInt_t endTime)
{
  // Sets the information of the run which is currently processed
  // can be overriden for special behaviour, make sure that you call base class
  // function

  fRun = run;
  fStartTime = startTime;
  fEndTime = endTime;
}

//______________________________________________________________________________________________
Bool_t AliPreprocessor::Store(const char* pathLevel2, const char* pathLevel3, TObject* object,
		AliCDBMetaData* metaData, Int_t validityStart, Bool_t validityInfinite)
{
  // Stores a CDB object in the storage for offline reconstruction. Objects that are not needed for
  // offline reconstruction, but should be stored anyway (e.g. for debugging) should NOT be stored
  // using this function. Use StoreReferenceData instead!
  //
  // This function should be called at the end of the preprocessor cycle
  //
  // The parameters are
  //   1, 2) the 2nd and 3rd level of the object's path. The first level is the detector name which is provided
  //         by the Preprocessor and converted to the Offline name. Thus the object's path is "DET/level2/level3"
  //   3) the object to be stored
  //   4) the metaData to be associated with the object
  //   5) the validity start run number w.r.t. the current run,
  //      if the data is valid only for this run leave the default 0
  //   6) specifies if the calibration data is valid for infinity (this means until updated),
  //      typical for calibration runs, the default is kFALSE
  //
  // The call is delegated to AliShuttleInterface

  const char* offlineDetName = AliShuttleInterface::GetOfflineDetName(GetName());
  if(!offlineDetName) return 0;

  return fShuttle->Store(AliCDBPath(offlineDetName, pathLevel2, pathLevel3), object,
  		metaData, validityStart, validityInfinite);
}

//______________________________________________________________________________________________
Bool_t AliPreprocessor::StoreReferenceData(const char* pathLevel2, const char* pathLevel3, TObject* object,
		AliCDBMetaData* metaData)
{
  // Stores a CDB object in the storage for reference data. This objects will not be available during
  // offline reconstrunction. Use this function for reference data only!
  //
  // This function should be called at the end of the preprocessor cycle
  //
  // The parameters are
  //   1, 2) the 2nd and 3rd level of the object's path. The first level is the detector name which is provided
  //         by the Preprocessor and converted to the Offline name. Thus the object's path is "DET/level2/level3"
  //   3) the object to be stored
  //   4) the metaData to be associated with the object
  //
  // The call is delegated to AliShuttleInterface

  const char* offlineDetName = AliShuttleInterface::GetOfflineDetName(GetName());
  if(!offlineDetName) return 0;

  return fShuttle->StoreReferenceData(AliCDBPath(offlineDetName, pathLevel2, pathLevel3), object,
  		metaData);
}
    
//______________________________________________________________________________________________
Bool_t AliPreprocessor::StoreReferenceFile(const char* localFile, const char* gridFileName)
{
	//
	// Stores a file directly (without opening it) in the reference storage in the Grid
	//
	// The file is stored under the following location: 
	// <base folder of reference storage>/<DET>/<RUN#>_<gridFileName>
	// where <gridFileName> is the second parameter given to the function
	//
	// The call is delegated to AliShuttleInterface
	
	const char* offlineDetName = AliShuttleInterface::GetOfflineDetName(GetName());
  	if(!offlineDetName) return 0;
	return fShuttle->StoreReferenceFile(GetName(), localFile, gridFileName);
}
    
//______________________________________________________________________________________________
Bool_t AliPreprocessor::StoreRunMetadataFile(const char* localFile, const char* gridFileName)
{
	//
	// Stores Run metadata file to the Grid, in the run folder
	//
	// Only GRP can call this function.
		
  	TString detName(GetName());
	if (detName != "GRP") 
	{
		Log("StoreRunMetadataFile - Only GRP can call this function.");
		return kFALSE;
	}
	return fShuttle->StoreRunMetadataFile(localFile, gridFileName);
}

//______________________________________________________________________________________________
const char* AliPreprocessor::GetFile(Int_t system, const char* id, const char* source)
{
  // This function retrieves a file from the given system (kDAQ, kDCS, kHLT) with the given file id
  // and from the given source in the system.
  // The function returnes the path to the local file.
  //
  // The call is delegated to AliShuttleInterface

  return fShuttle->GetFile(system, GetName(), id, source);
}

//______________________________________________________________________________________________
TList* AliPreprocessor::GetFileSources(Int_t system, const char* id)
{
  // Returns a list of sources in a given system that saved a file with the given id
  // if id is not given all sources are returned
  //
  // The call is delegated to AliShuttleInterface

  return fShuttle->GetFileSources(system, GetName(), id);
}

//______________________________________________________________________________________________
TList* AliPreprocessor::GetFileIDs(Int_t system, const char* source)
{
  // Returns a list of ids in a given system that saved a file with the given source
  //
  // The call is delegated to AliShuttleInterface

  return fShuttle->GetFileIDs(system, GetName(), source);
}

//______________________________________________________________________________________________
void AliPreprocessor::Log(const char* message, UInt_t level)
{
  // Adds a log message to the Shuttle log of this preprocessor
  //
  // The call is delegated to AliShuttleInterface

  fShuttle->Log(GetName(), message, level);
}

//______________________________________________________________________________________________
const char* AliPreprocessor::GetRunParameter(const char* param)
{
  // Return run parameter read from run logbook
  //
  // The call is delegated to AliShuttleInterface

  return fShuttle->GetRunParameter(param);
}

//______________________________________________________________________________________________
AliCDBEntry* AliPreprocessor::GetFromOCDB(const char* pathLevel2, const char* pathLevel3)
{
  // Return object from OCDB valid for current run
  //
  // The call is delegated to AliShuttleInterface

  const char* offlineDetName = AliShuttleInterface::GetOfflineDetName(GetName());
  if (!offlineDetName) return 0;

  return dynamic_cast<AliCDBEntry*>
  	(fShuttle->GetFromOCDB(GetName(), AliCDBPath(offlineDetName, pathLevel2, pathLevel3)));
}

//______________________________________________________________________________________________
AliCDBEntry* AliPreprocessor::GetGeometryFromOCDB()
{
  // Return geometry from OCDB (GRP/Geometry/Data) valid for current run
  //
  // The call is delegated to AliShuttleInterface

  return dynamic_cast<AliCDBEntry*>
  	(fShuttle->GetFromOCDB(GetName(), AliCDBPath("GRP", "Geometry", "Data")));
}

//______________________________________________________________________________________________
const char* AliPreprocessor::GetRunType()
{
  // Return run type string read from "run type" logbook
  //
  // The call is delegated to AliShuttleInterface

  return fShuttle->GetRunType();
}

//______________________________________________________________________________________________
Bool_t AliPreprocessor::GetHLTStatus()
{
  // Return HLT status (ON or OFF)
  // Converts the HLT status from the status string read in the run logbook (not just a bool)
  // The call is delegated to AliShuttleInterface

  return fShuttle->GetHLTStatus();

}
    
//______________________________________________________________________________________________
const char* AliPreprocessor::GetTriggerConfiguration()
{
	// Returns the trigger configuration which is read from a table in the DAQ logbook
	// The call is delegated to AliShuttleInterface
	// Only GRP can call this function.
	
	TString detName(GetName());
	if (detName != "GRP") 
	{
		Log("GetTriggerConfiguration - Only GRP can call this function.");
		return 0;
	}
	
	return fShuttle->GetTriggerConfiguration();
}
//______________________________________________________________________________________________
const char* AliPreprocessor::GetCTPTimeParams()
{
	// Returns the CTP timing parameters read from a table in the DAQ logbook
	// The call is delegated to AliShuttleInterface
	// Only GRP can call this function.
	
	TString detName(GetName());
	if (detName != "GRP") 
	{
		Log("GetCTPTimeParams - Only GRP can call this function.");
		return 0;
	}
	
	return fShuttle->GetCTPTimeParams();
}
//______________________________________________________________________________________________
const char* AliPreprocessor::GetTriggerDetectorMask()
{
	// Returns the trigger detector mask which is read from a table in the DAQ logbook
	// The call is delegated to AliShuttleInterface
	// Only TRI can call this function.
	
	TString detName(GetName());
	if (detName != "TRI") 
	{
		Log("GetTriggerDetectorMask - Only TRI can call this function.");
		return 0;
	}
	
	return fShuttle->GetTriggerDetectorMask();
}

//______________________________________________________________________________________________
void AliPreprocessor::AddRunType(const char* runType)
{
	// adds the given run type to the list of run types that are processed
	// this function should be called in the constructor of the derived preprocessor
	
	if (!runType)
		return;
	
	fRunTypes.Add(new TObjString(runType));
}
    
//______________________________________________________________________________________________
Bool_t AliPreprocessor::AliPreprocessor::ProcessRunType()
{
	// searches for the current run type in the list of run types that are processed by this
	// preprocessor. The list is populated by AddRunType
	
	const char* runType = GetRunType();

	Log(Form("Checking if run type %s is in the list of run types to be processed by this preprocessor...", runType));
	
	if (fRunTypes.GetEntries() == 0)
		Log("WARNING! There are no run types defined. This preprocessor will never run.");

	if (fRunTypes.FindObject(runType))
	{
		Log("Run type found. Processing this run.");
		return kTRUE;
	}
	
	Log("Run type not found. Skipping this run.");
	return kFALSE;
}
//______________________________________________________________________________________________
UInt_t AliPreprocessor::GetStartTimeDCSQuery()
{
	// Return Start Time for the DCS query
	//
	// The call is delegated to AliShuttleInterface

	return fShuttle->GetStartTimeDCSQuery();
}
//______________________________________________________________________________________________
UInt_t AliPreprocessor::GetEndTimeDCSQuery()
{
	// Return End Time for the DCS query
	//
	// The call is delegated to AliShuttleInterface

	return fShuttle->GetEndTimeDCSQuery();
}

//______________________________________________________________________________________________
const char* AliPreprocessor::GetForeignFile(const char* detector, Int_t system, const char* id, const char* source)
{
	// This function retrieves a file produced from the given detector (different from the 
	// current one) from the given system (kDAQ, kDCS, kHLT) 
	// with the given file id
	// and from the given source in the system.
	// The function returnes the path to the local file.
	//
	// The call is delegated to AliShuttleInterface
	// The function can be called only from the GRP preprocessor
	
	TString detName(GetName());
	if (detName != "GRP") 
	{
		Log("GetForeignFile - Only GRP can call this function.");
		return 0;
	}
	return fShuttle->GetFile(system, detector, id, source);
}

//______________________________________________________________________________________________
TList* AliPreprocessor::GetForeignFileSources(const char* detector, Int_t system, const char* id)
{
	// Returns the list of sources in a given system that produced the file 
	// with the given id for the foreign detector specified
	//
	// The call is delegated to AliShuttleInterface
	// The function can be called only from the GRP preprocessor

	TString detName(GetName());
	if (detName != "GRP") 
	{
		Log("GetForeignFileSources - Only GRP can call this function.");
		return 0;
	}
	return fShuttle->GetFileSources(system, detector, id);
}

