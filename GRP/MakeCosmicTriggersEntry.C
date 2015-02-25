#include "ARVersion.h"
#include <iostream>
#include <fstream>
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include <TROOT.h>
#include <TSystem.h>
#include <THashTable.h>
#include <TString.h>
#include <TError.h>
#include <TObjString.h>
#endif

Bool_t MakeCosmicTriggersEntry(const char *fileName, const char* cdbUri)
{
  const char* macroname = "MakeCosmicTriggersEntry.C";

  if (gSystem->AccessPathName(fileName)) {
    Error(macroname,Form("file (%s) not found", fileName));
    return kFALSE;
  }

  ifstream *file = new ifstream(fileName);
  if (!*file) {
    Error(macroname,Form("Error opening file (%s) !",fileName));
    file->close();
    delete file;
    return kFALSE;
  }

  THashTable *table = new THashTable();
  table->SetName("List of defined cosmic triggers");

  TString strLine;
  while (strLine.ReadLine(*file)) {

    if (strLine.BeginsWith("#")) continue;

    strLine.ReplaceAll(" ","");
    strLine.ReplaceAll("\t","");
    if (strLine.IsNull()) continue;

    TObjString *obj = new TObjString(strLine.Data());
    table->Add(obj);
  }

  file->close();
  delete file;


  // create OCDB storage
  TString Storage(cdbUri);
  if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
    Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
    return kFALSE;
  }
  AliCDBManager* cdb = AliCDBManager::Instance();
  AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
  if(!storage){
    Error(macroname,"Unable to open storage %s\n",Storage.Data());
    return kFALSE;
  }

  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Federico Antinori");
  md->SetComment("List of the defined cosmic triggers. It is used in order to steer the reconstruction, namely in the selection of the proper event specie. It is maintained and updated by the trigger coordinator.");
  // Get root and AliRoot versions and set them in the metadata
  const char* rootv = gROOT->GetVersion();
  TString av(ALIROOT_VERSION);
  TString revnum(ALIROOT_REVISION);
  av+=" - revision: ";
  av+=revnum;
  md->SetAliRootVersion(av.Data());

  AliCDBId id("GRP/Calib/CosmicTriggers",0,AliCDBRunRange::Infinity());
  Info(macroname,"Saving the list of defined cosmic triggers in the OCDB storage \"%s\"",Storage.Data());
  storage->Put(table,id,md);

  table->Delete();
  delete table;

  return kTRUE;
}
