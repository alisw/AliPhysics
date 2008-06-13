#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliITSresponseSDD.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#endif

void StoreRespSDD(Int_t firstRun=0, Int_t lastRun=999999999 ){
  ///////////////////////////////////////////////////////////////////////
  // Macro to generate and store the calibration files for SDD         //
  // Generates:                                                        //
  //  1 file with the AliITSrespionseSDD object (RespSDD)              //
  ///////////////////////////////////////////////////////////////////////
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://OCDB");
  }
  

  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetObjectClassName("AliITSresponse");
  md->SetResponsible("Francesco Prino");
  md->SetBeamPeriod(0);
  md->SetComment("Simulated data");


  AliCDBId idRespSDD("ITS/Calib/RespSDD",firstRun, lastRun);
  AliITSresponseSDD* rd = new AliITSresponseSDD();
  rd->SetTimeOffset(54.3);
  AliCDBManager::Instance()->GetDefaultStorage()->Put(rd, idRespSDD, md);  
}
