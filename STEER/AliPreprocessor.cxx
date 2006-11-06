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
#include <TList.h>
#include <TMap.h>

#include "AliLog.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBPath.h"
#include "AliShuttleInterface.h"

ClassImp(AliPreprocessor)

//______________________________________________________________________________________________
AliPreprocessor::AliPreprocessor(const char* detector, AliShuttleInterface* shuttle) :
  TNamed(detector, ""),
  fRun(-1),
  fStartTime(0),
  fEndTime(0),
  fShuttle(shuttle)
{
	SetTitle(Form("AliPreprocessor for %s subdetector.", detector));

  if (!fShuttle)
  {
    AliFatal("Initialized without Shuttle instance.");
    return;
  }

  fShuttle->RegisterPreprocessor(this);
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
UInt_t AliPreprocessor::Store(const char* pathLevel2, const char* pathLevel3, TObject* object,
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
UInt_t AliPreprocessor::StoreReferenceData(const char* pathLevel2, const char* pathLevel3, TObject* object,
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
  //
  // The call is delegated to AliShuttleInterface

  return fShuttle->GetFileSources(system, GetName(), id);
}

//______________________________________________________________________________________________
void AliPreprocessor::Log(const char* message)
{
  // Adds a log message to the Shuttle log of this preprocessor
  //
  // The call is delegated to AliShuttleInterface

  fShuttle->Log(GetName(), message);
}

//______________________________________________________________________________________________
const char* AliPreprocessor::GetRunParameter(const char* param)
{
  // Return run parameter read from run logbook
  //
  // The call is delegated to AliShuttleInterface

  return fShuttle->GetRunParameter(param);
}
