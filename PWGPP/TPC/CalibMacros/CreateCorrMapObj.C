#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliTPCChebCorr.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"

#include <TObjArray.h>

#endif

const char* GetObjPath(Bool_t correction, Bool_t ref);


//___________________________________________________
const char* GetObjPath(Bool_t correction, Bool_t ref) 
{
  // create object path
  return Form("TPC/Calib/%sMap%s",correction ? "Correction":"Distortion", ref ? "Ref":"");
}

//_________________________________________________________________________________________________
void CreateCorrMapObjRef(TObjArray* mapArr,int firstrun=0,int lastrun=-1, const char* dest="raw://")
{
  // create reference maps 
  int nmap = mapArr->GetEntriesFast();
  TObjArray* arr = new TObjArray();
  int nAdd = 0;
  int run = -1;
  Bool_t isCorrection = kTRUE;
  for (int i=0;i<nmap;i++) {
    AliTPCChebCorr* map = (AliTPCChebCorr*)mapArr->At(i); 
    if (!map) continue;
    isCorrection = map->IsCorrection(); // is this correction or distortion?
    if (run<0) run = map->GetRun();
    for (int j=0;j<nAdd;j++) { // check if the map of given type is unique
      AliTPCChebCorr* map0 = (AliTPCChebCorr*)arr->At(j); 
      if (map->GetFieldType()==map0->GetFieldType()) {
	printf("Error: the map of fieldType %d was already added, check input array\n",map->GetFieldType());
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
  man->SetDefaultStorage(dest);
  man->SetRun(run);
  //
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Ruben Shahoyan");
  md->SetComment(Form("Cheb. %s reference map",isCorrection ? "Correction":"Distortion"));
  AliCDBId id(GetObjPath(isCorrection,kTRUE),firstrun,lastrun<0 ? (AliCDBRunRange::Infinity()) : lastrun);
  man->Put(arr,id,md); 
  printf("Created an object with single time independent map\n");
  //
}

//_________________________________________________________________________________________________
void CreateCorrMapObjDef(TObjArray* mapArr,int firstrun=0,int lastrun=-1, const char* dest="raw://")
{
  // create default maps 
  int nmap = mapArr->GetEntriesFast();
  TObjArray* arr = new TObjArray();
  int nAdd = 0;
  int run = -1;
  Bool_t isCorrection = kTRUE;
  for (int i=0;i<nmap;i++) {
    AliTPCChebCorr* map = (AliTPCChebCorr*)mapArr->At(i); 
    if (!map) continue;
    isCorrection = map->IsCorrection(); // is this correction or distortion?
    if (run<0) run = map->GetRun();
    for (int j=0;j<nAdd;j++) { // check if the map of given type is unique
      AliTPCChebCorr* map0 = (AliTPCChebCorr*)arr->At(j); 
      if (map->GetFieldType()==map0->GetFieldType()) {
	printf("Error: the map of fieldType %d was already added, check input array\n",map->GetFieldType());
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
  printf("Added default maps for %d field types\n",nAdd);
  AliCDBManager* man = AliCDBManager::Instance();
  //
  man->UnsetDefaultStorage();
  man->SetDefaultStorage(dest);
  man->SetRun(run);
  //
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Ruben Shahoyan");
  md->SetComment(Form("Cheb. %s default map",isCorrection ? "Correction":"Distortion"));
  AliCDBId id(GetObjPath(isCorrection,kFALSE), firstrun,lastrun<0 ? (AliCDBRunRange::Infinity()) : lastrun);
  man->Put(arr,id,md); 
  printf("Created an object with single time independent map\n");
  //
}

//_________________________________________________________________________________________________
void CreateCorrMapObjDef(AliTPCChebCorr* map,int firstrun=0,int lastrun=-1, const char* dest="raw://")
{
  map->SetTimeStampStart(0);
  map->SetTimeStampEnd(0xffffffff);  
  map->SetTimeDependent(kFALSE);
  TObjArray* arr = new TObjArray();
  arr->Add(map);
  AliCDBManager* man = AliCDBManager::Instance();
  //
  Bool_t isCorrection = map->IsCorrection(); // is this correction or distortion?
  man->UnsetDefaultStorage();
  man->SetDefaultStorage(dest);
  man->SetRun(map->GetRun());
  //
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Ruben Shahoyan");
  md->SetComment(Form("Cheb. %s default map",isCorrection ? "Correction":"Distortion"));
  AliCDBId id(GetObjPath(isCorrection,kFALSE),firstrun,lastrun<0 ? (AliCDBRunRange::Infinity()) : lastrun);
  man->Put(arr,id,md); 
  printf("Created an object with single time independent map\n");
  //
}

//_________________________________________________________________________________________________
void CreateCorrMapObjRef(AliTPCChebCorr* map,int firstrun=0,int lastrun=-1, const char* dest="raw://")
{
  map->SetTimeStampStart(0);
  map->SetTimeStampEnd(0xffffffff);  
  map->SetTimeDependent(kFALSE);
  TObjArray* arr = new TObjArray();
  arr->Add(map);
  AliCDBManager* man = AliCDBManager::Instance();
  //
  Bool_t isCorrection = map->IsCorrection(); // is this correction or distortion?
  man->UnsetDefaultStorage();
  man->SetDefaultStorage(dest);
  man->SetRun(map->GetRun());
  //
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Ruben Shahoyan");
  md->SetComment(Form("Cheb. %s reference map",isCorrection ? "Correction":"Distortion"));
  AliCDBId id(GetObjPath(isCorrection,kTRUE),firstrun,lastrun<0 ? (AliCDBRunRange::Infinity()) : lastrun);
  man->Put(arr,id,md); 
  printf("Created an object with single time independent map\n");
  //
}

//_________________________________________________________________________________________________
void CreateCorrMapObj(AliTPCChebCorr* map,int firstrun=0,int lastrun=-1, const char* dest="raw://")
{
  map->SetTimeStampStart(0);
  map->SetTimeStampEnd(0xffffffff);  
  map->SetTimeDependent(kFALSE);
  TObjArray* arr = new TObjArray();
  arr->Add(map);
  AliCDBManager* man = AliCDBManager::Instance();
  //
  Bool_t isCorrection = map->IsCorrection(); // is this correction or distortion?
  man->UnsetDefaultStorage();
  man->SetDefaultStorage(dest);
  man->SetRun(map->GetRun());
  //
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Ruben Shahoyan");
  md->SetComment(Form("Cheb. %s map",isCorrection ? "Correction":"Distortion"));
  AliCDBId id(GetObjPath(isCorrection,kFALSE),firstrun,lastrun<0 ? (AliCDBRunRange::Infinity()) : lastrun);
  man->Put(arr,id,md); 
  printf("Created an object with single time independent map\n");
  //
}

//_________________________________________________________________________________________________
void CreateCorrMapObjTime(TObjArray* arrMap,int firstrun=0,int lastrun=-1, const char* dest="raw://")
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
  int run = -1;
  Bool_t isCorrection = kTRUE;
  for (int im0=0;im0<nmap;im0++) {
    AliTPCChebCorr* map0 = (AliTPCChebCorr*)arrMap->At(im0);
    if (run<0) run = map0->GetRun();
    map0->SetTimeDependent(kTRUE);
    isCorrection = map0->IsCorrection(); // is this correction or distortion?
  }
  // make sure time stamps are ordered
  for (int im0=0;im0<nmap;im0++) {
    AliTPCChebCorr* map0 = (AliTPCChebCorr*)arrMap->At(im0);
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
  man->SetDefaultStorage(dest);
  man->SetRun(run);
  //
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Ruben Shahoyan");
  md->SetComment(Form("Cheb. %s maps array",isCorrection ? "Correction":"Distortion"));
  AliCDBId id(GetObjPath(isCorrection,kFALSE),firstrun,lastrun<0 ? (AliCDBRunRange::Infinity()) : lastrun);
  man->Put(arrMap,id,md); 
  printf("Created an object with %d time dependent maps\n",nmap);
  //
}
