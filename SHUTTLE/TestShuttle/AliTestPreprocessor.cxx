#include "AliTestPreprocessor.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliTestDataDCS.h"

#include <TTimeStamp.h>

//
// This class is an example for a simple preprocessor.
// It takes data from DCS and passes it to the class AliTestDataDCS, which
// reformats its. This class is then written to the CDB.
//

ClassImp(AliTestPreprocessor)

AliTestPreprocessor::AliTestPreprocessor(const char* detector, AliShuttleInterface* shuttle) :
  AliPreprocessor(detector, shuttle),
  fData(0)
{
  // constructor
}

AliTestPreprocessor::~AliTestPreprocessor()
{
  // destructor
}

void AliTestPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliTestDataDCS object

  AliPreprocessor::Initialize(run, startTime, endTime);

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

	fData = new AliTestDataDCS(fRun, fStartTime, fEndTime);
}

Int_t AliTestPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into a AliTestDataDCS object

  if (!dcsAliasMap)
    return -1;

	fData->ProcessData(*dcsAliasMap);

  const char* fileName = GetFile(kDAQ, "PEDESTALS", "GDC");
  if (fileName)
    AliInfo(Form("Got the file %s, now we can extract some values.", fileName));
  // open file, extract some values, write them to fData

  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Alberto Colla");
	metaData.SetComment("This preprocessor fills an AliTestDataDCS object.");

	Store(fData, &metaData);
	delete fData;
	fData = 0;

  return 0;
}

