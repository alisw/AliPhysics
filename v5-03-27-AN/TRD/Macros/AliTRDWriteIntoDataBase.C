#if !defined( __CINT__) || defined(__MAKECINT__)

#include <iostream>

#include <vector>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TObject.h>
#include <TSystem.h>
#include <TH1F.h>

#include "AliCDBManager.h"
#include "AliCDBId.h"
#include "AliCDBStorage.h"
#include "AliCDBMetaData.h"
#include "AliTRDCalibra.h"
#include "AliTRDCalDet.h"

#endif
const Int_t firstrun = 0;
const Int_t lastrun = 0;
AliCDBStorage* gStorLoc = 0;
AliCDBMetaData* CreateMetaObject(const char* objectClassName);
void StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData);
Bool_t WriteIntoDataBase(TTree *tree, Int_t i);


AliCDBMetaData* CreateMetaObject(const char* objectClassName)
{
  AliCDBMetaData *md1= new AliCDBMetaData(); 
  md1->SetObjectClassName(objectClassName);
  md1->SetResponsible("Raphaelle Bailhache");
  md1->SetBeamPeriod(1);
  md1->SetAliRootVersion("05-06-00"); //root version
  md1->SetComment("The dummy values in this calibration file are for testing only");
  
  return md1;
}

void StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData)
{
  AliCDBId id1(cdbPath, firstrun, lastrun); 
  id1.SetVersion(0);
  gStorLoc->Put(object, id1, metaData); 
}


Bool_t AliTRDWriteIntoDataBase(TTree *tree, Int_t i){
  //
  // To use this macro, you have first to write the coefficients in a file
  // To take the resulted tree from the file and to give it to this macro
  // The macro will write the coefficient into a local database
  // in the directory where you run it
  //


  AliTRDCalDet* obj1 = 0;
  TObject* obj2 = 0;
  AliCDBMetaData *metadata = 0;
 
 
  
  //Single instance of AliCDBManager and AliTRDCalibra 
  AliCDBManager *man = AliCDBManager::Instance();
  if (!man) {
    printf("Could not get CDB Manager");
    return kFALSE;
  }
  gStorLoc = man->GetDefaultStorage();
  if (!gStorLoc) {
    printf("Could not get default Storage");
    return kFALSE;
  }
  AliTRDCalibra *calibra = AliTRDCalibra::Instance();


 
  // Create the new database object Det
  if(i != 2) {
    obj1 = calibra->CreateDetObjectTree(tree, i);
    obj2 = calibra->CreatePadObjectTree(tree, i, obj1);
  }
  if(i == 2) obj2 = calibra->CreatePadObjectTree(tree);
  

  //Store the info for the detector
  if(i != 2) {
    metadata = CreateMetaObject("AliTRDCalDet");
    if(i == 0) StoreObject("TRD/Calib/ChamberGainFactor",(TObject *) obj1, metadata);
    if(i == 1) StoreObject("TRD/Calib/ChamberVdrift",(TObject *) obj1, metadata);
    if(i == 3) StoreObject("TRD/Calib/ChamberT0",(TObject *) obj1, metadata);
    
  }

 
  //Store the info for the pad
  metadata = CreateMetaObject("AliTRDCalPad");
  if(i == 0) StoreObject("TRD/Calib/LocalGainFactor",(TObject *) obj2, metadata); 
  if(i == 1) StoreObject("TRD/Calib/LocalVdrift",(TObject *) obj2, metadata); 
  if(i == 3) StoreObject("TRD/Calib/LocalT0",(TObject *) obj2, metadata); 
  if(i == 2) StoreObject("TRD/Calib/PRFWidth",(TObject *) obj2, metadata); 


  return kTRUE;
}
