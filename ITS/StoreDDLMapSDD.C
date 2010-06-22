#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliITSDDLModuleMapSDD.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include <TObjArray.h>
#include <TRandom3.h>
#endif

void StoreDDLMapSDD(Int_t firstRun=0, Int_t lastRun=AliCDBRunRange::Infinity()){
  ///////////////////////////////////////////////////////////////////////
  // Macro to generate and store the DDL map for SDD                   //
  // Generates:                                                        //
  //  1 file with 1 AliITSDDLModuleMapSDD object (DDLmapSDD)           //
  ///////////////////////////////////////////////////////////////////////
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://OCDB");
  }
  

  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("AliITSDDLModuleMapSDD");
  md1->SetResponsible("Francesco Prino");
  md1->SetBeamPeriod(0);
  md1->SetAliRootVersion("Head 20 dec. 2007"); //root version
  md1->SetComment("This is a test");

  AliCDBId idDDLSDD("ITS/Calib/DDLMapSDD",firstRun, lastRun);
  AliITSDDLModuleMapSDD *ddlmap=new AliITSDDLModuleMapSDD();
  ddlmap->SetDefaultMap();
  AliCDBManager::Instance()->GetDefaultStorage()->Put(ddlmap, idDDLSDD, md1);  
}
