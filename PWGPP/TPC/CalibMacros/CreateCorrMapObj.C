#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliTPCChebCorr.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"

#include <TObjArray.h>

#endif


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
    CreateCorrMapObj((AliTPCChebCorr*)arrMap->At(0),firstrun,lastrun);
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
