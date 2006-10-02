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
Revision 1.4  2006/08/08 14:19:07  jgrosseo
Update to shuttle classes (Alberto)

- Possibility to set the full object's path in the Preprocessor's and
Shuttle's  Store functions
- Possibility to extend the object's run validity in the same classes
("startValidity" and "validityInfinite" parameters)
- Implementation of the StoreReferenceData function to store reference
data in a dedicated CDB storage.

Revision 1.3  2006/07/11 12:44:32  jgrosseo
adding parameters for extended validity range of data produced by preprocessor

Revision 1.2  2006/06/06 14:20:05  jgrosseo
o) updated test preprocessor (alberto)
o) added comments to example macro
o) test shuttle implements new interface

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

//
// test implementation of the AliShuttleInterface, to be used for local tests of preprocessors
//
// reads files from the local disk
// stores to local CDB
// logs to the screen
//

#include "AliTestShuttle.h"
#include "AliLog.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBMetaData.h"
#include "AliCDBPath.h"
#include "AliCDBId.h"
#include "AliPreprocessor.h"

#include <TMap.h>
#include <TList.h>
#include <TString.h>
#include <TObjString.h>

ClassImp(AliTestShuttle)

//______________________________________________________________________________________________
AliTestShuttle::AliTestShuttle(Int_t run, UInt_t startTime, UInt_t endTime) :
  fRun(run),
  fStartTime(startTime),
  fEndTime(endTime),
  fInputFiles(0),
  fPreprocessors(0),
  fDcsAliasMap(0)
{
  // constructor

  fInputFiles = new TMap;
  fPreprocessors = new TObjArray;
}

//______________________________________________________________________________________________
AliTestShuttle::~AliTestShuttle()
{
  // destructor

  delete fInputFiles;
  fInputFiles = 0;

  delete fPreprocessors;
  fPreprocessors = 0;

  delete fDcsAliasMap;
  fDcsAliasMap = 0;
}

//______________________________________________________________________________________________
UInt_t AliTestShuttle::Store(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData,
				Int_t validityStart, Bool_t validityInfinite)
{
  // Stores the CDB object
  // This function should be called at the end of the preprocessor cycle
  //
  // This implementation just stores it on the local disk, the full AliShuttle
  // puts it to the Grid FileCatalog

  Int_t startRun = fRun - validityStart;
  if(startRun < 0) {
	AliError("First valid run happens to be less than 0! Setting it to 0...");
	startRun=0;
  }

  Int_t endRun = -1;
  if(validityInfinite) {
	endRun = AliCDBRunRange::Infinity();
  } else {
	endRun = fRun;
  }

  AliCDBId id(path, startRun, endRun);

  return AliCDBManager::Instance()->Put(object, id, metaData);
}

//______________________________________________________________________________________________
UInt_t AliTestShuttle::StoreReferenceData(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData)
{
  // Stores the object as reference data
  // This function should be called at the end of the preprocessor cycle
  //
  // This implementation just stores it on the local disk, the full AliShuttle
  // puts it to the Grid FileCatalog

  AliCDBId id(path, fRun, fRun);

  return AliCDBManager::Instance()->GetStorage("local://ReferenceStorage")->Put(object, id, metaData);
}

//______________________________________________________________________________________________
const char* AliTestShuttle::GetFile(Int_t system, const char* detector, const char* id, const char* source)
{
  // This function retrieves a file from the given system (kDAQ, kDCS, kHLT) with the given file id
  // and from the given source in the system.
  // The function returnes the path to the local file.
  //
  // test implementation of GetFile
  // takes files from the local disks, files are passen in a TMap in the constructor

  TString key;
  key.Form("%s-%s-%s", fkSystemNames[system], detector, id);
  TPair* sourceListPair = dynamic_cast<TPair*> (fInputFiles->FindObject(key.Data()));
  TMap* sourceList = 0;
  if (sourceListPair)
    sourceList = dynamic_cast<TMap*> (sourceListPair->Value());
  if (!sourceList)
  {
    AliError(Form("Could not find any file in %s with id %s (%s)", fkSystemNames[system], id, key.Data()));
    return 0;
  }

  TPair* fileNamePair = dynamic_cast<TPair*> (sourceList->FindObject(source));
  TObjString* fileName = dynamic_cast<TObjString*> (fileNamePair->Value());
  if (!fileName)
  {
    AliError(Form("Could not find files from source %s in %s with id %s",
			source, fkSystemNames[system], id));
    return 0;
  }

  return fileName->GetString().Data();
}

//______________________________________________________________________________________________
TList* AliTestShuttle::GetFileSources(Int_t system, const char* detector, const char* id)
{
  // Returns a list of sources in a given system that saved a file with the given id
  //
  // test implementation of GetFileSources
  // takes files from the local disks, files are passen in a TMap in the constructor

  TString key;
  key.Form("%s-%s-%s", fkSystemNames[system], detector, id);
  TPair* sourceListPair = dynamic_cast<TPair*> (fInputFiles->FindObject(key.Data()));
  TMap* sourceList = 0;
  if (sourceListPair)
    sourceList = dynamic_cast<TMap*> (sourceListPair->Value());
  if (!sourceList)
  {
    AliError(Form("Could not find any file in %s with id %s (%s)", fkSystemNames[system], id, key.Data()));
    return 0;
  }

  TIterator* iter = sourceList->GetTable()->MakeIterator();
  TObject* obj = 0;
  TList* list = new TList;
  while ((obj = iter->Next()))
  {
    TPair* pair = dynamic_cast<TPair*> (obj);
    if (pair)
      list->Add(pair->Key());
  }

  delete iter;

  return list;
}

//______________________________________________________________________________________________
void AliTestShuttle::Log(const char* detector, const char* message)
{
  // test implementation of Log
  // just prints to the screen

  AliInfo(Form("%s: %s", detector, message));
}

//______________________________________________________________________________________________
void AliTestShuttle::AddInputFile(Int_t system, const char* detector, const char* id, const char* source, const char* fileName)
{
  // This function adds a file to the list of input files

  TString key;
  key.Form("%s-%s-%s", fkSystemNames[system], detector, id);
  TPair* sourceListPair = dynamic_cast<TPair*> (fInputFiles->FindObject(key.Data()));
  TMap* sourceList = 0;
  if (sourceListPair)
    sourceList = dynamic_cast<TMap*> (sourceListPair->Value());
  if (!sourceList)
  {
    sourceList = new TMap;
    fInputFiles->Add(new TObjString(key), sourceList);
  }

  sourceList->Add(new TObjString(source), new TObjString(fileName));
}

//______________________________________________________________________________________________
void AliTestShuttle::Process()
{
  // This function tests all preprocessors that are registered to it
  // All preprocessors get the same dcs alias map and have access to the same list of files.

  for (Int_t i=0; i<fPreprocessors->GetEntries(); ++i)
  {
    AliPreprocessor* preprocessor = dynamic_cast<AliPreprocessor*> (fPreprocessors->At(i));
    if (preprocessor)
    {
      preprocessor->Initialize(fRun, fStartTime, fEndTime);
      preprocessor->Process(fDcsAliasMap);
    }
  }
}

//______________________________________________________________________________________________
void AliTestShuttle::RegisterPreprocessor(AliPreprocessor* preprocessor)
{
  // registers a preprocessor

  fPreprocessors->Add(preprocessor);
}
