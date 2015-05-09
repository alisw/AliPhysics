#if !defined( __CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <TRandom.h>
#include <TSystem.h>
#include <TDatime.h>
#include <TFile.h>

#include <AliCDBManager.h>
#include <AliCDBStorage.h>
#include <AliCDBEntry.h>
#include <AliCDBMetaData.h>
#include <AliGeometry.h>
#include <AliPID.h>

#include "../TRD/AliTRDgeometry.h"

#include "../TRD/Cal/AliTRDCalROC.h"
#include "../TRD/Cal/AliTRDCalPad.h"
#include "../TRD/Cal/AliTRDCalDet.h"
#include "../TRD/AliTRDcalibDB.h"

#include <AliTRDCalOnlineGainTable.h>

#endif

AliCDBStorage* gStorLoc = 0;


AliCDBMetaData* CreateMetaObject(const char* objectClassName);
void StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData);
void AliTRDCreateLocalGain(Bool_t residual = kFALSE);



//___________________________________________________________________________________________________
AliCDBMetaData* CreateMetaObject(const char* objectClassName)
{
  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName(objectClassName);
  md1->SetResponsible("Annika Passfeld");
  //md1->SetBeamPeriod(1);
  md1->SetAliRootVersion("05-34-18"); //root version
  md1->SetComment("");
  
  return md1;
}
//___________________________________________________________________________________________________
void StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData)
{
  AliCDBId id1(cdbPath, 0, 999999999); 
  //id1.SetVersion(0);
  gStorLoc->Put(object, id1, metaData); 
}
//___________________________________________________________________________________________________
void AliTRDCreateLocalGain(Bool_t residual)
{
 
  
  //*************************************************************************

  
  AliCDBManager *man = AliCDBManager::Instance();
  gStorLoc = man->GetStorage("local://$ALICE_ROOT/OCDB");
  if (!gStorLoc)
    return;

  TObject* obj = 0;
  AliCDBMetaData* metaData = 0;

  //Det//////////////////////////////////////////////////////////////////
  
  metaData = CreateMetaObject("AliTRDCalOnlineGainTable");
  
  file = TFile::Open("$ALICE_ROOT/data/Gaintbl_Uniform_FGAN8_2015-01.root");
  AliTRDCalOnlineGainTable *cal = (AliTRDCalOnlineGainTable *) file->Get("AliTRDCalOnlineGainTable");
  StoreObject("TRD/Calib/Gaintbl_Uniform_FGAN8_2015-01", (TObject *)cal, metaData);


 
}
