#include "AliTestPreprocessor.h"

#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliTestDataDCS.h"

#include <TTimeStamp.h>
#include <TObjString.h>
#include <TList.h>

//
// This class is an example for a simple preprocessor.
// It takes data from DCS and passes it to the class AliTestDataDCS, which
// reformats its. This class is then written to the CDB.
//

ClassImp(AliTestPreprocessor)

//______________________________________________________________________________________________
AliTestPreprocessor::AliTestPreprocessor(AliShuttleInterface* shuttle) :
  AliPreprocessor("TPC", shuttle),
  fData(0)
{
  // constructor
  
  AddRunType("PHYSICS");
  AddRunType("CALIBRATION");
}

//______________________________________________________________________________________________
AliTestPreprocessor::~AliTestPreprocessor()
{
  // destructor
}

//______________________________________________________________________________________________
void AliTestPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliTestDataDCS object

  AliPreprocessor::Initialize(run, startTime, endTime);

	Log(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

	fData = new AliTestDataDCS(fRun, fStartTime, fEndTime);
}

//______________________________________________________________________________________________
Bool_t AliTestPreprocessor::ProcessDCS()
{
	//
	// decide here if DCS data is to be processed
	//
	
	// TODO implement a decision, e.g. based on the run type
	// In this example: Skip DCS if run type is CALIB
	if (strcmp(GetRunType(), "CALIB") == 0)
		return kFALSE;
	
	return kTRUE;
}

//______________________________________________________________________________________________
UInt_t AliTestPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into a AliTestDataDCS object

  if (!dcsAliasMap)
  {
  	Log("ERROR: No DCS map provided by SHUTTLE!");
  	return 1;
  }

  // The processing of the DCS input data is forwarded to AliTestDataDCS
  fData->ProcessData(*dcsAliasMap);

  // Example of how to retrieve the run type from the Shuttle, using GetRunType() function
  // TODO Here the run type for the "DET" detector must be set manually with SetInputRunType function,
  // in reality it will be read from the "run type" logbook!
  TString runType = GetRunType();
  Log(Form("Run type for run %d: %s", fRun, runType.Data()));

  // Example of how to retrieve the list of sources that produced the file with id DRIFTVELOCITY
  TList* sourceList = GetFileSources(kDAQ, "DRIFTVELOCITY");
  sourceList = GetFileSources(kHLT, "HLTData");
  if (!sourceList)
  {
  	Log("Error retrieving list of sources from FXS!");
  	return 1;
  }

  if (sourceList->GetEntries() == 0)
  {
  	Log("No sources found for id HLTData!");
  	// TODO What to do now depends on the expected behaviour of the
  	// online DA: if it expected to produce data for the FXS every run, then
  	// if no sources are found it means a problem happened and you should return error;
  	// if DA may or may not send files to FXS, then you shouldn't return error
  	// and go on with the analysis!
  	return 1;
  }

  // TODO We have the list of sources that produced the files with Id DRIFTVELOCITY.
  // Now we will loop on the list and we'll query the files one by one. 
  Log("The following sources produced files with the id DRIFTVELOCITY");
  sourceList->Print();

  TIter iter(sourceList);
  TObjString *source = 0;
  while((source=dynamic_cast<TObjString*> (iter.Next()))){
  	TString fileName = GetFile(kDAQ, "DRIFTVELOCITY", source->GetName());
  	if (fileName.Length() > 0)
    		Log(Form("Got the file %s, now we can extract some values.", fileName.Data()));
  }

  delete sourceList;

  // Example of retrieving files from HLT, including how to query HLT status

  Bool_t hltStatus = GetHLTStatus(); // 1 = HLT ON (=> Query HLT), 0 = HLT OFF (=> skip HLT query)

  if (hltStatus)
  {
  	Log("HLT is ON, let's query our files");
  	sourceList = GetFileSources(kHLT, "HLTData");
  	if (!sourceList)
  	{
  		Log("Error retrieving list of sources from FXS!");
		return 1;
  	}

  	if (sourceList->GetEntries() == 0)
  	{
  		Log("No sources found for id HLTData!");
		// TODO What to do now depends on the expected behaviour of the
		// online DA: if it expected to produce data for the FXS every run, then
		// if no sources are found it means a problem happened and you should return error;
		// if DA may or may not send files to FXS, then you shouldn't return error
		// and go on with the analysis!
		return 1;
  	}
  	Log("The following sources produced files with the id HLTData");

	sourceList->Print();

  	TIter iter(sourceList);
  	TObjString *source = 0;
  	while((source=dynamic_cast<TObjString*> (iter.Next()))){
  		TString fileName = GetFile(kHLT, "HLTData", source->GetName());
  		if (fileName.Length() > 0)
  			Log(Form("Got the file %s from HLT, now we can extract some values.", fileName.Data()));
  	}

  	delete sourceList;

  } else {
	Log("HLT is OFF, skipping query...");
  }

  // Example of how to retrieve the list of sources that produced files
  sourceList = GetFileSources(kDAQ);
  if (!sourceList)
  {
  	Log("Error: No sources found!");
	return 1;
  }
  
  Log("The following sources produced files");
  sourceList->Print();
  delete sourceList;
  
  // Example of how to retrieve the list of ids from a given source
  TList* idList = GetFileIDs(kDAQ, "LDC0");
  if (!idList)
  {
  	Log("Error: No IDs found!");
	return 1;
  }
  
  Log("The following ids are available");
  idList->Print();
  delete idList;
  
  // Example to store a file directly to the reference storage
  // Suppose we have queried the file from the FXS. Now the file is available locally and is called "file1.root".
  const char* refFileName="file1.root";
  if (!StoreReferenceFile(refFileName, "InputData.root"))
  	return 1;
  

  // Example of how to retrieve a run parameter using GetRunParameter function
  // TODO Here the parameter must be set manually with SetInputRunParameter function,
  // in reality it will be read from the run logbook!

  // note that the parameters are returned as character strings!
  const char* nEvents = GetRunParameter("totalEvents");
  if (nEvents) {
  	Log(Form("Number of events for run %d: %s",fRun, nEvents));
  } else {
	Log(Form("Number of events not put in logbook!"));
  }

  // Example of how to retrieve a condition object from OCDB

  AliCDBEntry *entry = GetFromOCDB("Calib", "Data");
  if (!entry)
  {
	Log("No object found in OCDB!");
  } else {
	Log("Got TPC/Calib/Data object from OCDB. The object's metadata is: ");
	entry->PrintMetaData();
  }


  //Now we have to store the final CDB file
  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("TPC expert");
	metaData.SetComment("This preprocessor fills an AliTestDataDCS object.");

	Bool_t result = Store("Calib", "Data", fData, &metaData, 0, 0);
	delete fData;
	fData = 0;

  if (!result)
  	return 1;
  
  return 0;
}

