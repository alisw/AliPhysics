/*

  Primary goal of the proposal was to provide functionality to browse and compare the content of the OCDB
  specified by different means.

  a.) galice.root               - as currently implemented by Ruben in MC (or any file with cdbMap, and cdbList)
  b.) AliESDs.root              - for the reconstructed data
  c.) ocdb snapshot             - as used for grid productions
  d.) TMap(s)                   - as used internally in galice.root and AliESDs,root  
  e.) log file (if possible)    - looks logs aways used similar syntax, tested and working
  f.) C macro                   - custom macro

  Content comparison should be done:
  a.) on the level of symbolic links 
  b.) on the level of content itself 
  - by by byte comparison dif
  - data member by data member comparison

  Implementation assumption:
  All input formats (a .. f) will  be converted to the TMap storages and TList if AliCDBIds 


  Example usage:
  export mcProd=/hera/alice/miranov/alice-tpc-notes/MCprodutions/code
  .L  $mcProd/AliOCDBToolkit.cxx+
d  MakeDiffExampleUseCase();


*/

// STD
#include <iostream>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <functional>
//
#include "TSystem.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"
//
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

//
//
//


void MakeDiffExampleUseCase();    // example usage
//
//
Bool_t   ParseInfoFromOcdbString(TString ocdbString, TString &ocdbPath, Int_t &run0, Int_t &run1, Int_t &version, Int_t &subVersion);   
Bool_t   ParseInfoFromOcdbString(TString ocdbString, AliCDBId &cdbId);
TList  * ConvertListStringToCDBId(const TList *cdbList0); 

void LoadOCDBFromLog(const char *logName, Int_t verbose);
void LoadOCDBFromMap(const TMap *cdbMap, const TList *cdbList);
void MakeDiff(const TMap *cdbMap0, const TList *cdbList0, const TMap *cdbMap1, const TList *cdbList1, Int_t verbose);



void MakeDiffExampleUseCase(){
  //
  // Example usage for the MC 
  // To run example case, assuming presence of following files in working directory: 
  //    - rec.log        
  //    - galice.root   
  //    - AliESDs.root
  //
  AliCDBManager * man = AliCDBManager::Instance();
  LoadOCDBFromLog("rec.log",0);
  const TMap *cdbMapLog= man->GetStorageMap();        // this is map of 
  const TList *cdbListLog=man->GetRetrievedIds();     // this is list of AliCDBId
  //  TList *cdbListLog0=man->GetRetrievedIds();     // this is list of AliCDBId
  //
  TFile *fmc = TFile::Open("galice.root");
  TMap *cdbMapMC= (TMap*)fmc->Get("cdbMap");          // 
  TList *cdbListMC0= (TList*)fmc->Get("cdbList");     // this is list of TObjStrings
  TList *cdbListMC = ConvertListStringToCDBId(cdbListMC0);        // convert to the TObjArray of AliCDBids
  //
  TFile *fesd = TFile::Open("AliESDs.root");
  TList *listESD = ((TTree*)fesd->Get("esdTree"))->GetUserInfo();
  TMap *cdbMapESD= (TMap*)listESD->FindObject("cdbMap");  
  TList *cdbListESD0= (TList*)listESD->FindObject("cdbList"); // this is list of TObjStrings
  TList *cdbListESD = ConvertListStringToCDBId(cdbListESD0);              // convert to the TObjArray  of AliCDBids
  //
  //
  //
  printf("\n\n");
  printf("Diff log>>>ESD\n\n:");
  MakeDiff(cdbMapLog, cdbListLog, cdbMapESD, cdbListESD,0);
  printf("\n\n");
  printf("Diff ESD>>>log\n\n:");
  MakeDiff(cdbMapESD, cdbListESD,cdbMapLog, cdbListLog,0);
  // 
  printf("\n\n");
  printf("Diff ESD>>>MC\n\n:");
  MakeDiff(cdbMapMC, cdbListMC, cdbMapESD, cdbListESD,0);
}




Bool_t ParseInfoFromOcdbString(TString ocdbString, TString &ocdbPath, Int_t &run0, Int_t &run1, Int_t &version, Int_t &subVersion){
  //
  // Parse OCDB id string and provide basic ocdb information
  //
  //  a.) parse ocdbPath
  Int_t indexBeginPath= ocdbString.Index("path: ")+7;
  if (indexBeginPath<0) return kFALSE;
  Int_t indexEndPath=ocdbString.Index(";",indexBeginPath);
  if (indexEndPath<0) return kFALSE;
  ocdbPath=TString(&(ocdbString.Data()[indexBeginPath]), indexEndPath-indexBeginPath-1);
  // b.) parse runRange
  Int_t indexRun0= ocdbString.Index(": [",indexEndPath)+3;
  if (indexRun0<0) return kFALSE;
  Int_t indexRun1= ocdbString.Index(",",indexRun0)+1;
  if (indexRun1<0) return kFALSE;
  run0=atoi(&(ocdbString.Data()[indexRun0]));
  run1=atoi(&(ocdbString.Data()[indexRun1]));
  AliCDBRunRange runRange(run0,run1);
  //c.) parse version, subversion
  Int_t indexVersion= ocdbString.Index("version: v",indexRun1)+10;
  if (indexVersion<0) return kFALSE;
  Int_t indexSubVersion= ocdbString.Index("_s",indexVersion)+2;
  if (indexSubVersion<0) return kFALSE;
  version=atoi(&(ocdbString.Data()[indexVersion]));
  subVersion=atoi(&(ocdbString.Data()[indexSubVersion]));
  return kTRUE;
}

Bool_t ParseInfoFromOcdbString(TString ocdbString, AliCDBId &cdbId){
  //
  // Parse OCDB id string and provide basic ocdb information and fillcdbID object
  //
  TString ocdbPath;
  Int_t run0=0, run1=0;
  Int_t version=0, subVersion=0;
  Bool_t parseStatus = ParseInfoFromOcdbString(ocdbString, ocdbPath, run0,run1,version,subVersion); 
  if (parseStatus) {
    AliCDBRunRange runRange(run0,run1);
    cdbId=AliCDBId(ocdbPath.Data(),runRange,version,subVersion);
    AliCDBId* id = AliCDBId::MakeFromString(ocdbString);
    cdbId=*id;
    delete id;
  }
  //
  return parseStatus;
}

TList  * ConvertListStringToCDBId(const TList *cdbList0){
  //
  // Convert input  list of the TObjString to list to AliCDBid 
  //
  Int_t entriesList0=cdbList0->GetEntries();
  TList * array0 = new TList();
  AliCDBId tmp0;
  for (Int_t ientry0=0; ientry0<entriesList0; ientry0++){
    if (cdbList0->At(ientry0)==0) continue;
    Bool_t isId =  cdbList0->At(ientry0)->IsA()->InheritsFrom("AliCDBId");
    if (isId){
      array0->AddLast(cdbList0->At(ientry0));
    }else{
      Bool_t isString =  cdbList0->At(ientry0)->IsA()->InheritsFrom("TObjString");
      if (isString){
	TObjString* sid0 = dynamic_cast<TObjString*> (cdbList0->At(ientry0));
	Bool_t status =  ParseInfoFromOcdbString(sid0->String(), tmp0);
	if (!status) continue;
	array0->AddLast(new AliCDBId(tmp0));
      }
    }
  }
  return array0;  
}



void LoadOCDBFromLog(const char *logName, Int_t verbose){
  //
  // Initilaize OCDB
  // Load OCDB setting as specified in log
  // Assuming fixed version of the log 
  // AliCDBManager is initilaized - ocdbMap and ID list can be exported
  //

  // Parsing/loading sequence:
  //    0.) SetDefault storage  *** Default Storage URI:
  //    1.) SetSpecific storage  *** Specific storage
  //    2.) SetRunNumber  Run number:
  //    3.) Set used IDs
  //
  AliCDBManager * man = AliCDBManager::Instance();
  //
  // 0.) SetDefault storage  *** Default Storage URI:
  // 
  TString  defaultOCDB = gSystem->GetFromPipe(TString::Format("cat %s| grep \"Storage URI:\"",logName).Data());
  TObjArray *array = defaultOCDB.Tokenize("\"");
  man->SetDefaultStorage(array->Last()->GetName());
  delete array;
  //
  // 1.) SetSpecific storage  *** Specific storage
  //
  TString  specificStorage  = gSystem->GetFromPipe(TString::Format("cat %s| grep \"Specific storage\"",logName).Data());
  array = specificStorage.Tokenize("\"");
  Int_t entries = array->GetEntries();
  for (Int_t i=1; i<entries-2; i+=4){    
    // add protection here line shuld be in expected format
    if (verbose&2>0) printf("%s\t%s\n",array->At(i)->GetName(),array->At(i+2)->GetName());    
    man->SetSpecificStorage(array->At(i)->GetName(),array->At(i+2)->GetName());
  }
  delete array;
  //
  // 2.) SetRunNumber  Run number:
  //
  TString  runLine  = gSystem->GetFromPipe(TString::Format("cat %s| grep \"I-AliCDBManager::Print: Run number =\"",logName).Data());
  array = runLine.Tokenize("=");
  Int_t run = 0;
  if (array->GetEntries()>1) run=atoi(array->At(1)->GetName()); 
  delete array;
  man->SetRun(run);  
  //
  // 3.) Set used IDs
  //   
  TString  ids =   gSystem->GetFromPipe(TString::Format("cat %s| grep I-AliCDB | grep path| grep range | grep version", logName).Data());
  array= ids.Tokenize("\n");
  entries = array->GetEntries();
  //
  for (Int_t i=0; i<entries; i++){
    //
    TString ocdbString = array->At(i)->GetName();
    TString ocdbEntry;
    TString ocdbPath;
    Int_t run0=0, run1=0;
    Int_t version=0, subVersion=0;
    Bool_t parseStatus = ParseInfoFromOcdbString(ocdbString, ocdbPath, run0,run1,version,subVersion); 
    if (!parseStatus) continue;
    AliCDBRunRange runRange(run0,run1);
    //
    if ((verbose&2)!=0) {
      printf("%s/Run%d_%d_v%d_s%d.root\n",ocdbPath.Data(),run0,run1,version,subVersion); 
    }
    try {
      AliCDBEntry*entry= man->Get(ocdbPath.Data(),runRange,version,subVersion);
    } catch(const exception &e){
      cerr << "OCDB retrieval failed!" << endl;
      cerr << "Detailes: " << e.what() << endl;
    }
  }  
  if ((verbose&1)!=0){
    man->Print();
    man->GetStorageMap()->Print();
    man->GetRetrievedIds()->Print(); 
  }
}


void LoadOCDBFromMap(TMap */*cdbMap*/, TList */*cdbList*/){
  //
  // Initilaize OCDB
  // Load OCDB setting as specified in maps
  // Or Do we have already implementation in AliCDBanager?
}



void MakeDiff(const TMap *cdbMap0, const TList *cdbList0, const TMap *cdbMap1, const TList *cdbList1, Int_t verbose){
  //
  //
  // Print difference between the 2 ocdb maps
  // Input:
  //   maps and list charactireizing OCDB setup  
  // Output:
  //   To be decided, currently it is the teinal output
  //
  Int_t entriesList0=cdbList0->GetEntries();
  Int_t entriesList1=cdbList1->GetEntries();
  //
  for (Int_t ientry0=0; ientry0<entriesList0; ientry0++){
    AliCDBId *id0    = dynamic_cast<AliCDBId*> (cdbList0->At(ientry0));
    AliCDBId *id1=0;
    for (Int_t ientry1=0; ientry1<entriesList1; ientry1++){
      AliCDBId *cid1    = dynamic_cast<AliCDBId*> (cdbList1->At(ientry1));
      //id0.Print();
      //cid1.Print();
      if (cid1->GetPath().Contains(id0->GetPath().Data())==0) continue;
      id1=cid1;
    }
    if (!id1) {
      printf("Missing entry\t");
      id0->Print();
      continue;
    }
    Bool_t isOK=kTRUE;
    if (id0->GetFirstRun()!= id1->GetFirstRun() ||id0->GetLastRun()!= id1->GetLastRun()){
      printf("Differrent run range\n");
      id0->Print();
      id1->Print();
    }    
    if (id0->GetVersion()!= id1->GetVersion() ||id0->GetSubVersion()!= id1->GetSubVersion()){
      printf("Differrent version\n");
      id0->Print();
      id1->Print();
    }    
  }
}

