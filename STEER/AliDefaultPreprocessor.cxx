#include <AliDefaultPreprocessor.h>

#include <TMap.h>

#include "AliCDBMetaData.h"
#include "AliSimpleValue.h"

//
// This class is used for all subdectectors that dont have an own
// preprocessor.
// The DCS data is written into CDB, no files are processed.
//

ClassImp(AliDefaultPreprocessor)

AliDefaultPreprocessor::AliDefaultPreprocessor(const char* detector, AliShuttleInterface* shuttle) :
AliPreprocessor(detector, shuttle)
{
}

AliDefaultPreprocessor::~AliDefaultPreprocessor()
{
}

Int_t AliDefaultPreprocessor::Process(TMap* dcsAliasMap)
{
  // store to default CDB object

  AliCDBMetaData metaData;
  metaData.SetProperty("StartTime",
      new AliSimpleValue(fStartTime));
  metaData.SetProperty("EndTime",
      new AliSimpleValue(fEndTime));
  metaData.SetComment("Automatically stored by AliDefaultPreprocessor!");

  Store(dcsAliasMap, &metaData);

  return 0;
}
