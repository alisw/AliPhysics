#include <AliDefaultPreprocessor.h>

#include <TMap.h>

#include "AliCDBMetaData.h"
#include "AliDCSValue.h"

//
// This class is used for all subdectectors that dont have an own
// preprocessor.
// The DCS data is written into CDB, no files are processed.
//

ClassImp(AliDefaultPreprocessor)

//______________________________________________________________________________________________
AliDefaultPreprocessor::AliDefaultPreprocessor(const char* detector, AliShuttleInterface* shuttle) :
AliPreprocessor(detector, shuttle)
{
}

//______________________________________________________________________________________________
AliDefaultPreprocessor::~AliDefaultPreprocessor()
{
}

//______________________________________________________________________________________________
UInt_t AliDefaultPreprocessor::Process(TMap* dcsAliasMap)
{
  // store to default CDB object

  AliCDBMetaData metaData;
  metaData.SetProperty("StartEndTime",
      new AliDCSValue(fStartTime, fEndTime));
  metaData.SetComment("Automatically stored by AliDefaultPreprocessor!");

  return Store(dcsAliasMap, &metaData);
}
