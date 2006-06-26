#include "AliRICHPreprocessor.h" //header

#include <AliCDBMetaData.h>
#include <AliDCSValue.h>
#include "AliLog.h"
#include "AliTestDataDCS.h"

#include <TTimeStamp.h>

//
// This class is an example for a simple preprocessor.
// It takes data from DCS and passes it to the class AliTestDataDCS, which
// reformats its. This class is then written to the CDB.
//

ClassImp(AliRICHPreprocessor)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHPreprocessor::Initialize(Int_t run, UInt_t startTime,UInt_t endTime)
{

  AliPreprocessor::Initialize(run, startTime, endTime);

  Log(Form("\n\tRun %d \n\tStartTime %s \n\tEndTime %s", run,
		TTimeStamp(startTime).AsString(),
		TTimeStamp(endTime).AsString()));

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
UInt_t AliRICHPreprocessor::Process(TMap* pDcsMap)
{

  if(!pDcsMap)  return 0;

  const char* fileName = GetFile(kDAQ, "PEDESTALS", "GDC");
  if(fileName) AliInfo(Form("Got the file %s, now we can extract some values.", fileName));

  TList* list = GetFileSources(kDAQ, "DRIFTVELOCITY");
  if (list){
    AliInfo("The following sources produced files with the id DRIFTVELOCITY");
    list->Print();
    delete list;
  }

  //Now we have to store the final CDB file
  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("RICH expert");
	metaData.SetComment("This data produced by AliRICHPreprocessor from simulated input.");

  UInt_t result = Store(pTempFreon, &metaData); //use AliPreprocessor::Store(), not allowed to use AliCDBManager directly

  return result;
}

