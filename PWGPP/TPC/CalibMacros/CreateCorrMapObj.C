#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliTPCChebCorr.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"

#include <TObjArray.h>

#endif

void CreateCorrMapObjRef(TObjArray* mapArr,int firstrun=0,int lastrun=-1, const char* dest="local://")
{
  /*
    Create reference maps OCDB object from array of Cheb.maps for different fields
    Example of creation of the object for ++ and -- polarities:
    .L CreateCorrMapObj.C+
    f0 = TFile::Open("proj244918_row/voxelResTree.root");
    chP = run244918_0_9999999999;
    chP->GetName();
    chP->SetTitle("ref map for ++ field from PbPb2015");
    chP->SetFieldType(AliTPCChebCorr::kFieldPos);
    
    f1 = TFile::Open("proj246392_row/voxelResTree.root");
    chM = run246392_0_9999999999;
    chM->SetTitle("ref map for -- field from PbPb2015");
    chM->SetFieldType(AliTPCChebCorr::kFieldNeg);
    TObjArray ar;
    ar.Add(chP);
    ar.Add(chM);
    CreateCorrMapObjRef(&ar);

  */
  int nmap = mapArr->GetEntriesFast();
  TObjArray* arr = new TObjArray();
  int nAdd = 0;
  for (int i=0;i<nmap;i++) {
    AliTPCChebCorr* map = (AliTPCChebCorr*)mapArr->At(i); 
    if (!map) continue;
    for (int j=0;j<nAdd;j++) { // check if the map of given type is unique
      AliTPCChebCorr* map0 = (AliTPCChebCorr*)arr->At(j); 
      if (map->GetFieldType()==map0->GetFieldType()) {
	printf("Error: the maps of fieldType %d was already added, check input array\n",map->GetFieldType());
	exit(1);
      }
    }
    //
    map->SetTimeStampStart(0);
    map->SetTimeStampEnd(0xffffffff);  
    map->SetTimeDependent(kFALSE);
    arr->Add(map);
    printf("Adding map for field type %d at slot %d\n",map->GetFieldType(),nAdd);
    nAdd++;
  }
  //
  printf("Added references maps for %d field types\n",nAdd);
  AliCDBManager* man = AliCDBManager::Instance();
  //
  man->UnsetDefaultStorage();
  man->SetDefaultStorage("local://");
  man->SetSpecificStorage("TPC/Calib/CorrectionMapsRef",dest);
  //
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Ruben Shahoyan");
  md->SetComment("Test Cheb.Correction map instead of ComposedCorrection");
  AliCDBId id("TPC/Calib/CorrectionMapsRef",firstrun,lastrun<0 ? (AliCDBRunRange::Infinity()) : lastrun);
  //AliCDBStorage* st = man->GetStorage("local//.");
  man->Put(arr,id,md); 
  printf("Created an object with single time independent map\n");
  //
}


void CreateCorrMapObjRef(AliTPCChebCorr* map,int firstrun=0,int lastrun=-1, const char* dest="local://")
{
  map->SetTimeStampStart(0);
  map->SetTimeStampEnd(0xffffffff);  
  map->SetTimeDependent(kFALSE);
  TObjArray* arr = new TObjArray();
  arr->Add(map);
  AliCDBManager* man = AliCDBManager::Instance();
  //
  man->UnsetDefaultStorage();
  man->SetDefaultStorage("local://");
  man->SetSpecificStorage("TPC/Calib/CorrectionMapsRef",dest);
  //
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Ruben Shahoyan");
  md->SetComment("Test Cheb.Correction map instead of ComposedCorrection");
  AliCDBId id("TPC/Calib/CorrectionMapsRef",firstrun,lastrun<0 ? (AliCDBRunRange::Infinity()) : lastrun);
  //AliCDBStorage* st = man->GetStorage("local//.");
  man->Put(arr,id,md); 
  printf("Created an object with single time independent map\n");
  //
}


void CreateCorrMapObj(AliTPCChebCorr* map,int firstrun=0,int lastrun=-1, const char* dest="local://")
{
  map->SetTimeStampStart(0);
  map->SetTimeStampEnd(0xffffffff);  
  map->SetTimeDependent(kFALSE);
  TObjArray* arr = new TObjArray();
  arr->Add(map);
  AliCDBManager* man = AliCDBManager::Instance();
  //
  man->UnsetDefaultStorage();
  man->SetDefaultStorage("local://");
  man->SetSpecificStorage("TPC/Calib/CorrectionMaps",dest);
  //
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Ruben Shahoyan");
  md->SetComment("Test Cheb.Correction map instead of ComposedCorrection");
  AliCDBId id("TPC/Calib/CorrectionMaps",firstrun,lastrun<0 ? (AliCDBRunRange::Infinity()) : lastrun);
  //AliCDBStorage* st = man->GetStorage("local//.");
  man->Put(arr,id,md); 
  printf("Created an object with single time independent map\n");
  //
}

void CreateCorrMapObjTime(TObjArray* arrMap,int firstrun=0,int lastrun=-1, const char* dest="local://")
{
  // emulate time dependent object
  //
  int nmap = arrMap->GetEntriesFast();
  if (nmap<1) return;
  if (nmap==1) {
    CreateCorrMapObj((AliTPCChebCorr*)arrMap->At(0),firstrun,lastrun,dest);
    return;
  }
  //
  // make sure time stamps are ordered
  for (int im0=0;im0<nmap;im0++) {
    AliTPCChebCorr* map0 = (AliTPCChebCorr*)arrMap->At(im0);
    map0->SetTimeDependent(kTRUE);
    for (int im1=im0+1;im1<nmap;im1++) {
      AliTPCChebCorr* map1 = (AliTPCChebCorr*)arrMap->At(im1);
      if (map1->GetTimeStampCenter() < map0->GetTimeStampCenter()) { // swap
	arrMap->AddAt(map1,im0);
	arrMap->AddAt(map0,im1);
	map0 = map1;
      }
    }
  }
  //
  // make sure there are no gaps 
  AliTPCChebCorr* map0 = (AliTPCChebCorr*)arrMap->At(0);  
  for (int im1=1;im1<nmap;im1++) {
    AliTPCChebCorr* map1 = (AliTPCChebCorr*)arrMap->At(im1);
    Long_t dft =  map1->GetTimeStampStart() - map0->GetTimeStampEnd();
    if (dft!=0) {
      printf("Mismatch of %ld s between END[%d]=%ld and START[%d]=%ld\n",
	     dft, im1-1, map0->GetTimeStampEnd(), im1, map1->GetTimeStampStart());
      //
      Long_t dft0 = dft/2;
      Long_t dft1 = dft-dft0;
      map0->SetTimeStampEnd(map0->GetTimeStampEnd() + dft0);
      map1->SetTimeStampStart(map1->GetTimeStampStart() - dft1);
      printf("Set END[%d]=%ld and START[%d]=%ld\n",im1-1,map0->GetTimeStampEnd(), im1,map1->GetTimeStampStart());
    }
    map0 = map1;
  }
  //
  arrMap->Print();
  //
  AliCDBManager* man = AliCDBManager::Instance();
  //
  man->UnsetDefaultStorage();
  man->SetDefaultStorage("local://");
  man->SetSpecificStorage("TPC/Calib/CorrectionMaps",dest);
  //  man->SetRun(run);
  //
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Ruben Shahoyan");
  md->SetComment("Test Cheb.Correction map instead of ComposedCorrection");
  AliCDBId id("TPC/Calib/CorrectionMaps",firstrun,lastrun<0 ? (AliCDBRunRange::Infinity()) : lastrun);
  //AliCDBStorage* st = man->GetStorage("local//.");
  man->Put(arrMap,id,md); 
  printf("Created an object with %d time dependent maps\n",nmap);
  //
}
