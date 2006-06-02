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

AliPreprocessor::AliPreprocessor(const char* detector, AliShuttleInterface* shuttle) :
  TNamed(detector, ""),
  fShuttle(shuttle)
{
	SetTitle(Form("AliPreprocessor for %s subdetector.", detector));

  if (!fShuttle)
    AliFatal("Initialized without Shuttle instance.");
}

AliPreprocessor::~AliPreprocessor()
{
}

void AliPreprocessor::Initialize(Int_t run, UInt_t startTime,	UInt_t endTime)
{
  // Sets the information of the run which is currently processed
  // can be overriden for special behaviour, make sure that you call base class
  // function

  fRun = run;
  fStartTime = startTime;
  fEndTime = endTime;
}

Int_t AliPreprocessor::Store(TObject* object, AliCDBMetaData* metaData)
{
  // Stores the CDB object
  // This function should be called at the end of the preprocessor cycle
  //
  // The call is delegated to AliShuttleInterface

  return fShuttle->Store(GetName(), object, metaData);
}

const char* AliPreprocessor::GetFile(Int_t system, const char* id, const char* source)
{
  // This function retrieves a file from the given system (kDAQ, kDCS, kHLT) with the given file id
  // and from the given source in the system.
  // The function returnes the path to the local file.
  //
  // The call is delegated to AliShuttleInterface

  return fShuttle->GetFile(system, GetName(), id, source);
}

TList* AliPreprocessor::GetFileSources(Int_t system, const char* id)
{
  // Returns a list of sources in a given system that saved a file with the given id
  //
  // The call is delegated to AliShuttleInterface

  return fShuttle->GetFileSources(system, GetName(), id);
}

void AliPreprocessor::Log(const char* message)
{
  // Adds a log message to the Shuttle log of this preprocessor
  //
  // The call is delegated to AliShuttleInterface

  fShuttle->Log(GetName(), message);
}
