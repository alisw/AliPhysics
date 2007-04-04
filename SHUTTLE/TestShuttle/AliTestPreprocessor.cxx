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
    return 1;

  // The processing of the DCS input data is forwarded to AliTestDataDCS
  fData->ProcessData(*dcsAliasMap);

  // Example of how to retrieve the run type from the Shuttle, using GetRunType() function
  // TODO Here the run type for the "DET" detector must be set manually with SetInputRunType function,
  // in reality it will be read from the "run type" logbook!
  TString runType = GetRunType();
  Log(Form("Run type for run %d: %s", fRun, runType.Data()));

  TString fileName = GetFile(kDAQ, "PEDESTALS", "GDC");
  if (fileName.Length() > 0)
    Log(Form("Got the file %s, now we can extract some values.", fileName.Data()));
  //TODO here the file could be opened, some values extracted and  written to e.g. fData

  //Example to store a file directly to the reference storage
  if (!StoreReferenceFile(fileName, "InputData.root"))
  	return 1;
  
  TList* list = GetFileSources(kDAQ, "DRIFTVELOCITY");
  if (list)
  {
    Log("The following sources produced files with the id DRIFTVELOCITY");
    list->Print();
    delete list;
  }
  //TODO here the files could be opened, some values extracted and  written to e.g. fData

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
	TObjString *obj = dynamic_cast<TObjString*> (entry->GetObject());
	Log(Form("Got TPC/Calib/Data object from OCDB. The object says: %s",obj->GetName()));
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

