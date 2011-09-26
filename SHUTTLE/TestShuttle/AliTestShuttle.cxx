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
Revision 1.14  2007/08/06 12:25:47  acolla
Function Bool_t GetHLTStatus added to preprocessor. It returns the status of HLT
read from the run logbook.
TestShuttle setup updated.
TRD data point configuration updated.

Revision 1.13  2007/05/30 06:35:21  jgrosseo
Adding functionality to the Shuttle/TestShuttle:
o) Function to retrieve list of sources from a given system (GetFileSources with id=0)
o) Function to retrieve list of IDs for a given source      (GetFileIDs)
These functions are needed for dealing with the tag files that are saved for the GRP preprocessor
Example code has been added to the TestProcessor in TestShuttle

Revision 1.12  2007/04/27 07:06:48  jgrosseo
GetFileSources returns empty list in case of no files, but successful query
No mails sent in testmode

Revision 1.11  2007/04/04 10:33:36  jgrosseo
1) Storing of files to the Grid is now done _after_ your preprocessors succeeded. This is transparent, which means that you can still use the same functions (Store, StoreReferenceData) to store files to the Grid. However, the Shuttle first stores them locally and transfers them after the preprocessor finished. The return code of these two functions has changed from UInt_t to Bool_t which gives you the success of the storing.
In case of an error with the Grid, the Shuttle will retry the storing later, the preprocessor does not need to be run again.

2) The meaning of the return code of the preprocessor has changed. 0 is now success and any other value means failure. This value is stored in the log and you can use it to keep details about the error condition.

3) New function StoreReferenceFile to _directly_ store a file (without opening it) to the reference storage.

4) The memory usage of the preprocessor is monitored. If it exceeds 2 GB it is terminated.

5) New function AliPreprocessor::ProcessDCS(). If you do not need to have DCS data in all cases, you can skip the processing by implemting this function and returning kFALSE under certain conditions. E.g. if there is a certain run type.
If you always need DCS data (like before), you do not need to implement it.

6) The run type has been added to the monitoring page

Revision 1.10  2007/02/28 10:41:01  acolla
Run type field added in SHUTTLE framework. Run type is read from "run type" logbook and retrieved by
AliPreprocessor::GetRunType() function.
Added some ldap definition files.

Revision 1.8  2007/02/13 11:22:25  acolla
Shuttle getters and setters of main/local OCDB/Reference storages, temp and log
folders moved to AliShuttleInterface

Revision 1.6  2006/11/06 14:22:47  jgrosseo
major update (Alberto)
o) reading of run parameters from the logbook
o) online offline naming conversion
o) standalone DCSclient package

Revision 1.5  2006/10/02 12:58:52  jgrosseo
Small interface change in StoreReferenceData

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
#include "AliLTUConfig.h"
#include "AliDAQ.h"
#include "AliTriggerInput.h"

#include <TMap.h>
#include <TList.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <TObjArray.h>

ClassImp(AliTestShuttle)

//______________________________________________________________________________________________
AliTestShuttle::AliTestShuttle(Int_t run, UInt_t startTime, UInt_t endTime) :
  fRun(run),
  fStartTime(startTime),
  fEndTime(endTime),
  fTimeCreated(startTime),
  fDCSQueryOffset(0),
  fInputFiles(0),
  fRunParameters(0),
  fRunType(),
  fPreprocessors(0),
  fDcsAliasMap(0),
  fTriggerConfiguration(""),
  fTriggerDetectorMask(""),
  fCTPtiming(""),
  fltuConfig(0x0)
{
  // constructor

  fInputFiles = new TMap;
  fRunParameters = new TMap;
  fPreprocessors = new TObjArray;

  fInputFiles->SetOwner(1);
  fRunParameters->SetOwner(1);
  fPreprocessors->SetOwner(1);

  fltuConfig = new TObjArray();
  fltuConfig->SetOwner(1);
}

//______________________________________________________________________________________________
AliTestShuttle::~AliTestShuttle()
{
  // destructor

  delete fInputFiles;
  fInputFiles = 0;

  delete fRunParameters;
  fRunParameters = 0;

  delete fPreprocessors;
  fPreprocessors = 0;

  delete fDcsAliasMap;
  fDcsAliasMap = 0;

  delete fltuConfig;
  fltuConfig = 0;
}

//______________________________________________________________________________________________
Bool_t AliTestShuttle::Store(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData,
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

  return AliCDBManager::Instance()->GetStorage(fgkMainCDB)->Put(object, id, metaData);
}

//______________________________________________________________________________________________
Bool_t AliTestShuttle::StoreReferenceData(const AliCDBPath& path, TObject* object, AliCDBMetaData* metaData)
{
  // Stores the object as reference data
  // This function should be called at the end of the preprocessor cycle
  //
  // This implementation just stores it on the local disk, the full AliShuttle
  // puts it to the Grid FileCatalog

  AliCDBId id(path, fRun, fRun);

  return AliCDBManager::Instance()->GetStorage(fgkMainRefStorage)->Put(object, id, metaData);
}

//______________________________________________________________________________________________
Bool_t AliTestShuttle::StoreReferenceFile(const char* detector, const char* localFile, const char* gridFileName)
{
	//
	// Stores reference file directly (without opening it). 
	//
	// This implementation just stores it on the local disk, the full AliShuttle 
	// puts it to the Grid FileCatalog
	
	AliCDBManager* man = AliCDBManager::Instance();
	AliCDBStorage* sto = man->GetStorage(fgkLocalRefStorage);
	
	TString localBaseFolder = sto->GetBaseFolder();
	
	TString targetDir = GetRefFilePrefix(localBaseFolder.Data(), detector);	
	
	return CopyFileLocally(targetDir, localFile, gridFileName);
}

//______________________________________________________________________________________________
Bool_t AliTestShuttle::StoreRunMetadataFile(const char* localFile, const char* gridFileName)
{
	//
	// Stores Run metadata file to the Grid, in the run folder
	//
	// Only GRP can call this function.
	
	AliCDBManager* man = AliCDBManager::Instance();
	AliCDBStorage* sto = man->GetStorage(fgkLocalRefStorage);
	
	TString localBaseFolder = sto->GetBaseFolder();
	
	// Build Run level folder
	// folder = /alice/data/year/lhcPeriod/runNb/Raw
	
	TTimeStamp startTime(fStartTime);
		
	TString year =  Form("%d",startTime.GetDate());
	year = year(0,4);
		
	TString lhcPeriod = GetRunParameter("LHCperiod");
	
	if (lhcPeriod.Length() == 0) 
	{
		Log("SHUTTLE","StoreRunMetaDataFile - LHCPeriod not found in logbook!");
		return 0;
	}
	
	// TODO: currently SHUTTLE cannot write in /alice/data/ !!!!!
	//TString targetDir = Form("%s/GRP/RunMetadata/alice/data/%s/%s/%d/Raw", 
	//			localBaseFolder.Data(), year.Data(), 
	//			lhcPeriod.Data(), fRun);
	
	TString targetDir = Form("%s/GRP/RunMetadata/alice/simulation/%s/%s/%d/Raw", 
				localBaseFolder.Data(), year.Data(), 
				lhcPeriod.Data(), fRun);
					
	return CopyFileLocally(targetDir, localFile, gridFileName);
}

//______________________________________________________________________________________________
Bool_t AliTestShuttle::CopyFileLocally(TString& targetDir, const char* localFile, const char* gridFileName)
{
	//
	// Stores file locally. Called by StoreReferenceFile and StoreRunMetadataFile
	//
	
	//try to open folder, if it does not exist
	void* dir = gSystem->OpenDirectory(targetDir.Data());
	if (dir == NULL) {
		if (gSystem->mkdir(targetDir.Data(), kTRUE)) {
			Log("SHUTTLE", Form("StoreFileLocally - Can't open directory <%s>", targetDir.Data()));
			return kFALSE;
		}

	} else {
		gSystem->FreeDirectory(dir);
	}

	TString target = Form("%s/%s", targetDir.Data(), gridFileName);
	
	Int_t result = gSystem->GetPathInfo(localFile, 0, (Long64_t*) 0, 0, 0);
	if (result)
	{
		Log("SHUTTLE", Form("StoreFileLocally - %s does not exist", localFile));
		return kFALSE;
	}

	result = gSystem->CopyFile(localFile, target);

	if (result == 0)
	{
		Log("SHUTTLE", Form("StoreFileLocally - File %s stored locally to %s", localFile, target.Data()));
		return kTRUE;
	}
	else
	{
		Log("SHUTTLE", Form("StoreFileLocally - Could not store file %s to %s!. Error code = %d", 
				localFile, target.Data(), result));
		return kFALSE;
	}	



}

//______________________________________________________________________________________________
const char* AliTestShuttle::GetRefFilePrefix(const char* base, const char* detector)
{
	//
	// Get folder name of reference files 
	//

	TString offDetStr(GetOfflineDetName(detector));
	static TString dir;
	if (offDetStr == "ITS" || offDetStr == "MUON" || offDetStr == "PHOS")
	{
		dir.Form("%s/%s/%s", base, offDetStr.Data(), detector);
	} else {
		dir.Form("%s/%s", base, offDetStr.Data());
	}
	
	return dir.Data();
}

//______________________________________________________________________________________________
const char* AliTestShuttle::GetFile(Int_t system, const char* detector, const char* id, const char* source)
{
  // This function retrieves a file from the given system (kDAQ, kDCS, kHLT, kDQM) with the given file id
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

  TObjString* fileName = 0;
  TPair* fileNamePair = dynamic_cast<TPair*> (sourceList->FindObject(source));
  if (fileNamePair)
  	fileName = dynamic_cast<TObjString*> (fileNamePair->Value());
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
  if (id)
    key.Form("%s-%s-%s", fkSystemNames[system], detector, id);
  else
    key.Form("%s-%s", fkSystemNames[system], detector);
  
  TList* list = new TList;
  
  TIterator* iter = fInputFiles->MakeIterator();
  TObject* obj = 0;
  while ((obj = iter->Next()))
  {
	TObjString* objStr = dynamic_cast<TObjString*> (obj);
	if (objStr)
	{
		Bool_t found = kFALSE;
 		if (id)
		{
			found = (objStr->String().CompareTo(key) == 0);
		}
		else
			found = objStr->String().BeginsWith(key);
		
		if (found)
		{
			TPair* sourceListPair = dynamic_cast<TPair*> (fInputFiles->FindObject(objStr->String().Data()));
			TMap* sourceList = dynamic_cast<TMap*> (sourceListPair->Value());
	
			TIterator* iter2 = sourceList->GetTable()->MakeIterator();
			TObject* obj2 = 0;
			while ((obj2 = iter2->Next()))
			{
				TPair* pair = dynamic_cast<TPair*> (obj2);
				if (pair)
				{
					if (!list->FindObject(pair->Key()))
						list->Add(new TObjString(pair->Key()->GetName()));
				}
			}
			
			delete iter2;
		}
	}
  }
  
  if (list->GetEntries() == 0)
    AliInfo(Form("Could not find any file in %s with id %s (%s)", fkSystemNames[system], id, key.Data()));
  
  return list;
}

//______________________________________________________________________________________________
TList* AliTestShuttle::GetFileIDs(Int_t system, const char* detector, const char* source)
{
  // Returns a list of ids in a given system that saved a file with the given source
  //
  // test implementation of GetFileSources
  // takes files from the local disks, files are passen in a TMap in the constructor


  TString key;
  key.Form("%s-%s", fkSystemNames[system], detector);
  
  TList* list = new TList;
  
  TIterator* iter = fInputFiles->MakeIterator();
  TObject* obj = 0;
  while ((obj = iter->Next()))
  {
	TObjString* objStr = dynamic_cast<TObjString*> (obj);
	if (objStr)
	{
		if (objStr->String().BeginsWith(key))
		{
			Bool_t found = kFALSE;
		
			TPair* sourceListPair = dynamic_cast<TPair*> (fInputFiles->FindObject(objStr->String().Data()));
			TMap* sourceList = dynamic_cast<TMap*> (sourceListPair->Value());
	
			TIterator* iter2 = sourceList->GetTable()->MakeIterator();
			TObject* obj2 = 0;
			while ((obj2 = iter2->Next()))
			{
				TPair* pair = dynamic_cast<TPair*> (obj2);
				if (pair)
				{
					if (strcmp(pair->Key()->GetName(), source) == 0)
						found = kTRUE;
				}
			}
			
			delete iter2;
			
			if (found)
			{
				TObjArray* tokens = objStr->String().Tokenize("-");
				if (tokens->GetEntries() == 3)
				{
					TObjString* id = dynamic_cast<TObjString*> (tokens->At(2));
					if (id && !list->FindObject(id->String()))
						list->Add(new TObjString(id->String()));
				}
				
				delete tokens;
	
			}
		}
	}
  }
  
  if (list->GetEntries() == 0)
    AliInfo(Form("Could not find any file in %s with source %s (%s)", fkSystemNames[system], source, key.Data()));
  
  return list;
}

//______________________________________________________________________________________________
void AliTestShuttle::Log(const char* detector, const char* message, UInt_t level)
{
  // test implementation of Log
  // just prints to the screen

  TString fullMessage = detector;
  fullMessage.Append(Form(": %s",message));
  AliLog::Message(level, fullMessage, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);
}

//______________________________________________________________________________________________
void AliTestShuttle::AddInputFile(Int_t system, const char* detector, const char* id, const char* source, const char* fileName)
{
  //
  // This function adds a file to the list of input files
  // the list is stored in fInputFiles 
  // fInputFiles: TMap (key -> value)
  //    <system>-<detector>-<id> -> TMap (key -> value)
  //                                <source> -> <filename>
  //  
  
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
Bool_t AliTestShuttle::AddInputCDBEntry(AliCDBEntry* entry)
{
  // This function adds an object in the OCDB to be later retrieved with GetFromOCDB

	AliCDBStorage *sto = AliCDBManager::Instance()->GetStorage(fgkMainCDB);
	if (!sto)
	{
		Log("SHUTTLE", "GetFromOCDB - Cannot activate main OCDB for query!");
		return 0;
	}

	return sto->Put(entry);
}

//______________________________________________________________________________________________
AliCDBEntry* AliTestShuttle::GetFromOCDB(const char* detector, const AliCDBPath& path)
{
// returns obiect from OCDB valid for current run

	AliCDBStorage *sto = AliCDBManager::Instance()->GetStorage(fgkMainCDB);
	if (!sto)
	{
		Log("SHUTTLE", "GetFromOCDB - Cannot activate main OCDB for query!");
		return 0;
	}

	return (AliCDBEntry*) sto->Get(path, fRun);
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
      if (preprocessor->ProcessRunType())
	      {
		      preprocessor->Initialize(fRun, fStartTime, fEndTime);
		      preprocessor->Process(fDcsAliasMap);
	      }
      
    }
  }
}

//______________________________________________________________________________________________
void AliTestShuttle::RegisterPreprocessor(AliPreprocessor* preprocessor)
{
  // registers a preprocessor

	const char* detName = preprocessor->GetName();
	if(strcmp("DET", detName) != 0) {
		if(GetDetPos(detName) < 0)
			AliFatal(Form("********** !!!!! Invalid detector name: %s !!!!! **********", detName));
	}

  	fPreprocessors->Add(preprocessor);
}

//______________________________________________________________________________________________
void AliTestShuttle::AddInputRunParameter(const char* key, const char* value){
// set a run parameter (in reality it will be read from the DAQ logbook)

	TObjString* keyObj = new TObjString(key);
	if (fRunParameters->Contains(key)) {
		AliWarning(Form("Parameter %s already existing and it will be replaced.", key));
		delete fRunParameters->Remove(keyObj);

	}
	fRunParameters->Add(keyObj, new TObjString(value));
	AliDebug(2, Form("Number of parameters: %d", fRunParameters->
	GetEntries()));
}

//______________________________________________________________________________________________
const char* AliTestShuttle::GetRunType()
{
	//
	// get a run parameter
	//

	return fRunType.Data();
}

//______________________________________________________________________________________________
const char* AliTestShuttle::GetRunParameter(const char* key){
// get a run parameter

	TObjString* value = dynamic_cast<TObjString*> (fRunParameters->GetValue(key));
	if(!value) {
		AliError(Form("No such parameter: %s", key));
		return 0;
	}
	return value->GetName();
}

//______________________________________________________________________________________________
void AliTestShuttle::SetShuttleTempDir(const char* tmpDir)
{
// sets Shuttle temp directory

	fgkShuttleTempDir = gSystem->ExpandPathName(tmpDir);
}

//______________________________________________________________________________________________
void AliTestShuttle::SetShuttleLogDir(const char* logDir)
{
// sets Shuttle log directory

	fgkShuttleLogDir = gSystem->ExpandPathName(logDir);
}

//______________________________________________________________________________________________
const char* AliTestShuttle::GetTriggerConfiguration()
{
	//returns trigger configuration
	if (fTriggerConfiguration.Length()>0){
		return fTriggerConfiguration;
	}
	return NULL;
}
//______________________________________________________________________________________________
const char* AliTestShuttle::GetCTPTimeParams()
{
	//returns trigger configuration
	if (fCTPtiming.Length()>0){
		return fCTPtiming;
	}
	return NULL;
}
//______________________________________________________________________________________________
const char* AliTestShuttle::GetTriggerDetectorMask()
{
	//returns trigger detector mask
	if (fTriggerDetectorMask.Length()>0){
		return fTriggerDetectorMask;
	}
	return NULL;
}
//______________________________________________________________________________________________
UInt_t AliTestShuttle::GetStartTimeDCSQuery()
{
	// Return Start Time for the DCS query
	//
	// The call is delegated to AliShuttleInterface

	return fTimeCreated-fDCSQueryOffset;
}
//______________________________________________________________________________________________
UInt_t AliTestShuttle::GetEndTimeDCSQuery()
{
	// Return End Time for the DCS query
	//
	// The call is delegated to AliShuttleInterface

	return fEndTime+fDCSQueryOffset;
}
//______________________________________________________________________________________________
void AliTestShuttle::SendMLFromDet(const char* value)
{
	// 
	// Sending an information coming from the current detector to ML
	//
	
	Printf("%s will be sent to monalisa in the $currentdetector_RunCondition tag",value); 
	return;
}
//______________________________________________________________________________________________
void AliTestShuttle::SetLTUConfig(TString* ltuConfig, const char* det){

	// 
	// Setting LTU configuration for detector det
	//

	AliInfo(Form("LTU Config will be added for detector %s",det));
	AliInfo(Form("First element of the array of strings will correspond to LTUFineDelay1 --> %s",ltuConfig[0].Data()));
	AliInfo(Form("Second element of the array of strings will correspond to LTUFineDelay2 --> %s",ltuConfig[1].Data()));
	AliInfo(Form("Third element of the array of strings will correspond to LTUBCDelayAdd --> %s",ltuConfig[2].Data()));
	Float_t ltuFineDelay1 = ltuConfig[0].Atof();
	Float_t ltuFineDelay2 = ltuConfig[1].Atof();
	Float_t ltuBCDelayAdd = ltuConfig[2].Atof();
	AliLTUConfig* ltu = new AliLTUConfig((UChar_t)AliDAQ::DetectorID(det),ltuFineDelay1,ltuFineDelay2,ltuBCDelayAdd);
	Int_t idet = AliDAQ::DetectorID(det);
	fltuConfig->AddAtAndExpand(ltu,idet);
	return;
}
//______________________________________________________________________________________________
TString* AliTestShuttle::GetLTUConfig(const char* det){

	// 
	// Getting LTU configuration for detector det
	//

	TString* ltuConfigString = new TString[3];
	Int_t idet = -1;
	for (Int_t index = 0; index < AliDAQ::kNDetectors; index++){
		AliDebug(3,Form("index = %d, det = %s, CTP name = %s",index,det,AliTriggerInput::fgkCTPDetectorName[index]));
		if (strcmp(det,AliTriggerInput::fgkCTPDetectorName[index]) == 0){
			AliInfo(Form("Getting LTU configuration for detector %s",AliTriggerInput::fgkCTPDetectorName[index]));
			idet = index;
			break;
		}
	}
	if (idet != -1){
		AliLTUConfig* ltu = (AliLTUConfig*)fltuConfig->At(idet);
		if (!ltu){
			AliInfo(Form("ltu for detector %s not added in the simulated run, returning a null pointer",det));
			return 0x0;
		}
		else{
			ltuConfigString[0]=Form("%f",ltu->GetFineDelay1());
			ltuConfigString[1]=Form("%f",ltu->GetFineDelay2());
			ltuConfigString[0]=Form("%f",ltu->GetBCDelayAdd());
		}
	}
	else{
		AliError(Form("Detector %s not found in the list of CTP detector names",det));
	}
	return ltuConfigString;

}


