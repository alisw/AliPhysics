#include "AliZDCPreprocessor.h"

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliZDCDataDCS.h"

#include <TTimeStamp.h>

//
// This class is an example for a simple preprocessor.
// It takes data from DCS and passes it to the class AliZDCDataDCS, which
// reformats its. This class is then written to the CDB.
//

ClassImp(AliZDCPreprocessor)

//______________________________________________________________________________________________
AliZDCPreprocessor::AliZDCPreprocessor(const char* detector, AliShuttleInterface* shuttle) :
  AliPreprocessor(detector, shuttle),
  fData(0)
{
  // constructor
}

//______________________________________________________________________________________________
AliZDCPreprocessor::~AliZDCPreprocessor()
{
  // destructor
}

//______________________________________________________________________________________________
void AliZDCPreprocessor::Initialize(Int_t run, UInt_t startTime,
	UInt_t endTime)
{
  // Creates AliZDCDataDCS object

  AliPreprocessor::Initialize(run, startTime, endTime);

	AliInfo(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

	fData = new AliZDCDataDCS(fRun, fStartTime, fEndTime);
}

//______________________________________________________________________________________________
UInt_t AliZDCPreprocessor::Process(TMap* dcsAliasMap)
{
  // Fills data into a AliZDCDataDCS object

  if (!dcsAliasMap)
    return 0;

  // The processing of the DCS input data is forwarded to AliZDCDataDCS
	fData->ProcessData(*dcsAliasMap);

  const char* PedfileName = GetFile(kDAQ, "PEDESTALS", "LDC0");
  if(PedfileName){
    AliInfo(Form("Got the file %s, now we can extract some values.", PedfileName));
    //TODO here the file could be opened, some values extracted and  written to e.g. fData
  }
  else AliInfo(Form("File %s not found", PedfileName));

  const char* EMDfileName = GetFile(kDAQ, "MUTUALEMD", "GDC");
  if(EMDfileName){
    AliInfo(Form("Got the file %s, now we can extract some values.", EMDfileName));
    //TODO here the file could be opened, some values extracted and  written to e.g. fData
  }
  else AliInfo(Form("File %s not found", EMDfileName));

  //Now we have to store the final CDB file
  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Chiara");
	metaData.SetComment("This preprocessor fills an AliZDCDataDCS object.");

	UInt_t result = Store(fData, &metaData);
	delete fData;
	fData = 0;

  return result;
}

