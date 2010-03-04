#if !defined (__CINT__) || defined (__MAKECINT__)

#include "AliCDBManager.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "TObjString.h"
#endif

TString readCDBentry(TString arg){
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0);
  AliCDBEntry* entry = AliCDBManager::Instance()->Get(arg.Data());
  TObjString* objstr = dynamic_cast<TObjString*>(entry->GetObject());
  return objstr->GetString();
}
