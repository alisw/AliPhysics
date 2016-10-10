/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

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

  


  //=============================================================================
  // Functionality to dump content of  objects in human readable format in root prompt
  //=============================================================================
  Use case examples: 
  1.) compare content of alignent OCDB files for different yers
  2.) compare ClusterParam for different periods
  =================================================================================================================
  // 1.) Compare alignment example:
  //    Compare TPC alignemnt 2013 and 2010
  AliOCDBtoolkit::DumpOCDBFile("/cvmfs/alice-ocdb.cern.ch/calibration/data/2013/OCDB/TPC/Align/Data/Run0_999999999_v1_s0.root","TPCalign2013.dump",1,"docdb");
  AliOCDBtoolkit::DumpOCDBFile("/cvmfs/alice-ocdb.cern.ch/calibration/data/2010/OCDB/TPC/Align/Data/Run0_999999999_v1_s0.root","TPCalign2010.dump",1,"docdb");
  diff  TPCalign2013.dump TPCalign2010.dump > TPCalign2013_TPCalign2010.diff
  =================================================================================================================
  //  2.) Compare ClusterParam OCDB etry
  AliOCDBtoolkit::DumpOCDBFile("/cvmfs/alice-ocdb.cern.ch/calibration/data/2010/OCDB/TPC/Calib/ClusterParam/Run131541_999999999_v2_s0.root","2010_TPC_Calib_ClusterParam_Run131541_999999999_v2_s0.dump",1,"docdb");
  AliOCDBtoolkit::AliOCDBtoolkit::DumpOCDBFile("/cvmfs/alice-ocdb.cern.ch/calibration/data/2010/OCDB/TPC/Calib/ClusterParam/Run0_999999999_v1_s0.root","2010_TPC_Calib_ClusterParam_Run0_999999999_v1_s0.dump",1,"docdb");
  AliOCDBtoolkit::DumpOCDBFile("/cvmfs/alice-ocdb.cern.ch/calibration/data/2013/OCDB/TPC/Calib/ClusterParam/Run0_999999999_v1_s0.root","2013_TPC_Calib_ClusterParam_Run0_999999999_v1_s0.dump",1,"docdb");
  diff 2010_TPC_Calib_ClusterParam_Run131541_999999999_v2_s0.dump 2010_TPC_Calib_ClusterParam_Run0_999999999_v1_s0.dump


  //=============================================================================
  // Infovke OCDB toolkit functionality from prompt
  //=============================================================================
  More examples : AliOCDBtoolkit.sh in shell prompt (see examples in $ALICE_PHYSCS/PWGPP/CalibMacros/AliOCDBtoolkit.sh) 
  E.g: 
  1.) source  AliOCDBtoolkit.sh
      source $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.sh
  2.) make human readable OCDB tables, CDBids trees, OCDB obejct dumps .... specifying approprate input (ESD,galice.root, log files, snapshot) 
      (object dump only for selected items - configured by $ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.cfg
      (Examples for logs and snapshots to be rechecked)
    
      ocdbMakeTable alien:///alice/data/2010/LHC10d/000126158/pass4/10000126158023.10/AliESDs.root ESD OCDBrec.list   
               ===> Output: OCDBrec.list, OCDBrec.list.root *.xml and *.print files for selected OCDB entries 
      ocdbMakeTable alien:///alice/sim/2016/LHC16d1/244983/301/galice.root MC OCDBsim.list
               ===> Output: OCDBsim.list, OCDBsim.list.root  and *.print files for selected OCDB entries
      dumpObject  $ALICE_ROOT/OCDB/TPC/Calib/PadNoise/Run0_999999999_v1_s0.root  "object" "XML"  TPC_Calib_PadNoise_Run0_999999999_v1_s0.xml


*/



// STD
#include <iostream>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <functional>
#include "TRealData.h"
#include "TDataMember.h"
#include "TClass.h"
#include "TROOT.h"
#include <TVectorD.h>
#include <TObjArray.h>
#include "TRegexp.h"
#include "TPRegexp.h"
//
#include "TSystem.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"
#include "TMessage.h"
#include <TGrid.h>
#include <TEnv.h>
//
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliOCDBtoolkit.h"
#include "AliCDBStorage.h"
#include "TTreeStream.h"
#include "AliProdInfo.h"

using std::cout;
using std::cerr;
using std::endl;
using std::exception;


TObjArray *  AliOCDBtoolkit::fgExcludeList=0;
TObjArray *  AliOCDBtoolkit::fgXmlOCDBDumpList=0;
TObjArray *  AliOCDBtoolkit::fgPrintOCDBDumpList=0;
Int_t        AliOCDBtoolkit::fgVerbose=1;
Int_t        AliOCDBtoolkit::fgRun=1;
TString        AliOCDBtoolkit::fgPath="";
AliProdInfo * AliOCDBtoolkit::fgProdInfo=0;          // production information

void AliOCDBtoolkit::MakeDiffExampleUseCase(){
  //
  // Example usage for the MC 
  // To run example case, assuming presence of following files in working directory: 
  //    - rec.log        
  //    - galice.root   
  //    - AliESDs.root
  //
  AliCDBManager * man = AliCDBManager::Instance();
  AliOCDBtoolkit::LoadOCDBFromLog("rec.log",0);
  const TMap *cdbMapLog= man->GetStorageMap();        // this is map of 
  const TList *cdbListLog=man->GetRetrievedIds();     // this is list of AliCDBId
  //  TList *cdbListLog0=man->GetRetrievedIds();     // this is list of AliCDBId
  //
  TFile *fmc = TFile::Open("galice.root");
  TMap *cdbMapMC= (TMap*)fmc->Get("cdbMap");          // 
  TList *cdbListMC0= (TList*)fmc->Get("cdbList");     // this is list of TObjStrings
  TList *cdbListMC = AliOCDBtoolkit::ConvertListStringToCDBId(cdbListMC0);        // convert to the TObjArray of AliCDBids
  //
  TFile *fesd = TFile::Open("AliESDs.root");
  TList *listESD = ((TTree*)fesd->Get("esdTree"))->GetUserInfo();
  TMap *cdbMapESD= (TMap*)listESD->FindObject("cdbMap");  
  TList *cdbListESD0= (TList*)listESD->FindObject("cdbList"); // this is list of TObjStrings
  TList *cdbListESD = AliOCDBtoolkit::ConvertListStringToCDBId(cdbListESD0);              // convert to the TObjArray  of AliCDBids
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



void AliOCDBtoolkit::DumpOCDBAsTxt(const TString fInput, const TString fType, const TString outfile){
  //
  //
  //
  TFile *file;
  const TMap *cdbMap=0;
  const TList *cdbList=0;
  Bool_t cvmfsExist =  gSystem->AccessPathName("/cvmfs/alice-ocdb.cern.ch/calibration/")==kFALSE; // check presence of the OCDB in cvfms
  if (!fgProdInfo) fgProdInfo=new AliProdInfo;
  //
  //
  AliCDBManager * man = AliCDBManager::Instance();
  if (fInput.Contains("alien://") && gGrid==0){    
    TGrid *myGrid = TGrid::Connect("alien://");            //Oddly this will return also a pointer if connection fails
    if(myGrid->GetPort()==0){                       //if connection fails port 0 is saved, using this to check for successful connection
      cerr << "Cannot connect to grid!" << endl;
      return;
    }
  }
  if(fType.EqualTo("MC",TString::kIgnoreCase)){
        file = TFile::Open(fInput.Data());
	fgPath=fInput;
        cdbMap = (TMap*)file->Get("cdbMap");
	if (!cdbMap){
	  printf("cdbMap does not exist in input file\t%s. Exiting\n",fInput.Data());
	  return;
	}
	// 	
	TString cdbStorage=((TPair*)cdbMap->FindObject("default"))->Value()->GetName();
	AliOCDBtoolkit::CleanCDBPath(cdbStorage,cvmfsExist); // use cvmfs instead of alien: if available
	man->SetDefaultStorage(cdbStorage.Data());	
        TList *cdbListMC0 = (TList*)file->Get("cdbList");     // this is list of TObjStrings
        cdbList = AliOCDBtoolkit::ConvertListStringToCDBId(cdbListMC0);        // convert to the TObjArray of AliCDBids
  } 
    else if(fType.EqualTo("ESD",TString::kIgnoreCase)){
      file = TFile::Open(fInput.Data());
      fgPath=fInput;
      if (!file) {
	printf("Input file  does not exist %s. Exiting\n",fInput.Data());
	return;
      }
      TTree * tree =  ((TTree*)file->Get("esdTree"));
      if (!tree || (tree->Draw("fRunNumber","","goff",1)<=0)){
	::Error("AliOCDBtoolkit::DumpOCDBAsTxt","Input esdTree does not exist in %s",fInput.Data());
	return;
      }
      fgRun = tree->GetV1()[0];
      TList *listESD = tree->GetUserInfo();
      if (listESD)  fgProdInfo->Init(tree->GetUserInfo());
      cdbMap = (TMap*)listESD->FindObject("cdbMap");  
      if (!cdbMap){
	printf("cdbMap does not exist in input file\t%s. Exiting\n",fInput.Data());
	return;
      }
      AliOCDBtoolkit::SetStorage(cdbMap);
      TList *cdbListESD0= (TList*)listESD->FindObject("cdbList"); // this is list of TObjStrings
      cdbList = ConvertListStringToCDBId(cdbListESD0);              // convert to the TObjArray  of AliCDBids
    }
    else if(fType.EqualTo("log",TString::kIgnoreCase)){
        LoadOCDBFromLog(fInput.Data(),0);
        cdbMap = man->GetStorageMap();        // this is map of 
        cdbList =man->GetRetrievedIds();     // this is list of AliCDBId
    }
    else{
        printf("unsupported option: %s",fType.Data());
        return;
    }
  cout <<"BEGINDUMP:" << endl;
  DumpOCDB(cdbMap,cdbList,outfile);
}


Bool_t AliOCDBtoolkit::ParseInfoFromOcdbString(TString ocdbString, TString &ocdbPath, Int_t &run0, Int_t &run1, Int_t &version, Int_t &subVersion){
  // Functionalit
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

Bool_t AliOCDBtoolkit::ParseInfoFromOcdbString(TString ocdbString, AliCDBId &cdbId){
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

TList  * AliOCDBtoolkit::ConvertListStringToCDBId(TList *cdbList0){
  //
  // Convert input  list of the TObjString to list to AliCDBid 
  //
  cdbList0->Sort();
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



void AliOCDBtoolkit::LoadOCDBFromLog(const char *logName, Int_t verbose){
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
    if ((verbose&2)>0) printf("%s\t%s\n",array->At(i)->GetName(),array->At(i+2)->GetName());    
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
      man->Get(ocdbPath.Data(),runRange,version,subVersion);      
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

void  AliOCDBtoolkit::SetStorage(const TMap *cdbMap){
  //
  //  Set storages as specified in the map - TO CHECK.. 
  //  Should go to the AliCDBmanager if not alreadyhhere +++MI
  //   
  //  In case OCDB_PATH local variable is defined
  //  alien storage is replaced by OCDB_PATH prefix: e.g:  local:///cvmfs/alice.gsi.de/
  //
  //  Regexp extensivelly used - see documentation in ????
  //       http://wwwacs.gantep.edu.tr/guides/programming/root/htmldoc/examples/tstring.C.html
  AliCDBManager * man = AliCDBManager::Instance();  
  TIter iter(cdbMap->GetTable());
  TPair* aPair=0;
  Bool_t cvmfsExist =  gSystem->AccessPathName("/cvmfs/alice-ocdb.cern.ch/calibration/")==kFALSE; // check presence of the OCDB in cvfms


  while ((aPair = (TPair*) iter.Next())) {
    //    aPair->Value();
    //aPair->Print();
    TString urlOrig = aPair->Value()->GetName();
    TString url=urlOrig; // e.g TString  url="alien://?User=?DBFolder=/alice/data/2010/OCDB?SE=default?CacheFolder=?OperateDisconnected=1?CacheSize=1073741824?CleanupInterval=0"
    man->ExtractBaseFolder(url); // url==alien://Folder=/alice/data/2010/OCDB"
    TString ocdbPrefix(gSystem->Getenv("OCDB_PATHTEST"));
    if (ocdbPrefix.Length()>0){
      TRegexp alienPrefix("^alien://Folder=");      
      url(alienPrefix)=ocdbPrefix+"";
    }
    AliOCDBtoolkit::CleanCDBPath(url,cvmfsExist); // use cvmfs instead of alien: if available

    printf("%s\t%s\t%s\n", aPair->GetName(), urlOrig.Data(), url.Data());
    if (TString(aPair->GetName())=="default") man->SetDefaultStorage(url);
    else
      man->SetSpecificStorage(aPair->GetName(), url);
  }  
}
 
void AliOCDBtoolkit::LoadOCDBFromMap(const TMap *cdbMap, const TList *cdbList){
  //
  // Initilaize OCDB
  // Load OCDB setting as specified in maps
  // Or Do we have already implementation in AliCDBanager?  TO CHECK.. Should go to the AliCDBmanager if not alreadyhhere
  if (!fgExcludeList)  LoadAliOCDBtoolkitSetup(1);
  AliCDBManager * man = AliCDBManager::Instance();  
  AliOCDBtoolkit::SetStorage(cdbMap);  
  TIter iter(cdbList);
  TObjString *ocdbString=0;
  while (( ocdbString= (TObjString*) iter.Next())) {
    AliCDBId* cdbId = AliCDBId::MakeFromString(ocdbString->String());
    if (fgVerbose>0) cdbId->Print();
    if (IsEntrySelected( cdbId->GetPath(),fgExcludeList)) continue;
    try {
      //      AliCDBEntry * cdbEntry = (AliCDBEntry*) man->Get(*cdbId,kTRUE);
      
      man->Get(*cdbId);
    } catch(const exception &e){
      cerr << "OCDB retrieval failed!" << endl;
      cerr << "Detailes: " << e.what() << endl;
    }   
  }    
}

void AliOCDBtoolkit::LoadOCDBFromESD(const char *fname){
  //
  // Load OCDB setup from the ESD file
  // 
  TFile * fesd = TFile::Open(fname);
  TList *listESD = ((TTree*)fesd->Get("esdTree"))->GetUserInfo();
  TMap *cdbMapESD= (TMap*)listESD->FindObject("cdbMap");  
  TList *cdbListESD0= (TList*)listESD->FindObject("cdbList"); // this is list of TObjStrings
  AliOCDBtoolkit::SetStorage(cdbMapESD); 
  AliOCDBtoolkit::LoadOCDBFromMap(cdbMapESD, cdbListESD0);
}


void AliOCDBtoolkit::MakeDiff(const TMap *cdbMap0, const TList *cdbList0, const TMap */*cdbMap1*/, const TList *cdbList1, Int_t /*verbose*/){
  //
  //
  // Print difference between the 2 ocdb maps
  // Input:
  //   maps and list charactireizing OCDB setup  
  // Output:
  //   To be decided.
  //
  AliOCDBtoolkit::SetStorage(cdbMap0);
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
    //   Bool_t isOK=kTRUE;
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

void AliOCDBtoolkit::DumpOCDB(const TMap *cdbMap0, const TList *cdbList0, const TString outfile){
  //
  // Dump the OCDB configuatation as formated text file 
  // with following collumns
  // cdb name  prefix cdb path
  // OCDB entries are sorted alphabetically
  // e.g:
  // TPC/Calib/RecoParam /hera/alice/jwagner/software/aliroot/AliRoot_TPCdev/OCDB/ TPC/Calib/RecoParam/Run0_999999999_v0_s0.root $SIZE_AliCDBEntry_Object $HASH_AliCDBEntry_Object
  if (!fgExcludeList)  {
    LoadAliOCDBtoolkitSetup(1);
  }
  SetXRDTimeOutAll(-1);
  AliCDBManager * man = AliCDBManager::Instance();
  AliOCDBtoolkit::SetStorage(cdbMap0);  
  TList * cdbList = (TList*) cdbList0;   // sorted array
  AliCDBId *CDBId=0;
  TString cdbName="";
  TString cdbPath="";
  TObjString *ostr;
  AliCDBEntry *cdbEntry=0;
  TGrid *myGrid = NULL;
  UInt_t hash=0;
  TMessage * file;
  Int_t size=0; 
  FILE *ofs = fopen(outfile.Data(),"w");
  TTreeSRedirector *pcstream = new TTreeSRedirector((outfile+".root").Data(),"recreate");
  TObjString spass,speriod,syear,srun;
  Int_t year=0;
  if (fgPath.Length()>0){ // in case fgPath haexist decompose year/period/pass/run information
    TObjArray *mArr=0;
    TPRegexp regPass("/.?ass.+?/");
    TPRegexp regYear("/.?20.+?/");
    TPRegexp regPeriod("/.?LHC.+?/");
    TPRegexp regRun("/.?000.+?/");
    if (regPass.Match(fgPath,"i")>0) {mArr=regPass.MatchS(fgPath,"i"); spass=mArr->At(0)->GetName(); delete mArr;}
    if (regYear.Match(fgPath,"i")>0) {mArr=regYear.MatchS(fgPath,"i"); syear=mArr->At(0)->GetName(); delete mArr;}
    if (regPeriod.Match(fgPath,"i")>0) {mArr=regPeriod.MatchS(fgPath,"i"); speriod=mArr->At(0)->GetName(); delete mArr;}
    if (regRun.Match(fgPath,"i")>0) {mArr=regRun.MatchS(fgPath,"i"); srun=mArr->At(0)->GetName(); delete mArr;}
    spass.String().ReplaceAll("/","");
    syear.String().ReplaceAll("/","");
    speriod.String().ReplaceAll("/","");
    srun.String().ReplaceAll("/000","");
    srun.String().ReplaceAll("/","");
    if (srun.String().Atoi()>0) fgRun=srun.String().Atoi();
    if (syear.String().Atoi()>0) year=syear.String().Atoi();
    
  }

  for (Int_t iter=0;iter<2; iter++){ 
    size=0; 
    hash=0;
    // in first dump the table without detailed infromation
    // in second iteration dump full requested information including size/hash/metadata and xml/prin if specified
    TIter next(cdbList);
    while ((CDBId  =(AliCDBId*) next())){
      cdbName = CDBId->GetPath();
      ostr = (TObjString*)cdbMap0->GetValue(cdbName.Data());
      if(!ostr) ostr = (TObjString*)cdbMap0->GetValue("default");
      cdbPath = ostr->GetString();
      if(cdbPath.Contains("local://"))cdbPath=cdbPath(8,cdbPath.Length()).Data();
      if(!myGrid && cdbPath.Contains("alien://")){        //check if connection to alien is initialized
        myGrid = TGrid::Connect("alien://");            //Oddly this will return also a pointer if connection fails
        if(myGrid->GetPort()==0){                       //if connection fails port 0 is saved, using this to check for successful connection
	  cerr << "Cannot connect to grid!" << endl;
	  continue;
        }
      }
      if (iter > 0 && !(IsEntrySelected( cdbName.Data(),fgExcludeList))) {  //get detailed information
	try {
	  ::Info("AliOCDBtoolkit::DumpOCDB","%s",cdbName.Data());
	  cdbEntry = (AliCDBEntry*) man->Get(*CDBId,kFALSE);
	}catch(const exception &e){
	  cerr << "OCDB retrieval failed!" << endl;
	  cerr << "Detailes: " << e.what() << endl;
	  CDBId->Print();
	  hash=0;
	  size=-1;
	}  
	if (!cdbEntry) {
	  printf("Object not avaliable\n");
	  CDBId->Print();
	  continue;
	}
	TObject *obj = cdbEntry->GetObject();
	file = new TMessage(TBuffer::kWrite);
	file->WriteObject(obj);
	size = file->Length();	
	if(!obj){
	  fprintf(ofs,"object %s empty!\n",cdbName.Data());
	  continue;
	}
	hash = TString::Hash(file->Buffer(),size);
	delete file;
	if (IsEntrySelected( cdbName.Data(),fgXmlOCDBDumpList)){
	  ::Info("AliOCDBtoolkit::DumpOCDB","XML dump %s",cdbName.Data());
	  gSystem->mkdir(cdbName.Data(),1);
	  TString fname=TString::Format("%s/Run%d_%d_v%d_s%d.xml", cdbName.Data(), CDBId->GetFirstRun(), CDBId->GetLastRun(),
				       CDBId->GetVersion(),CDBId->GetSubVersion());
	  TFile * f = TFile::Open(fname.Data(),"recreate");
	  cdbEntry->Write("AliCDBEntry");
	  f->Close();

	  // XML dump
	}
	if (IsEntrySelected( cdbName.Data(),fgPrintOCDBDumpList)){  // !!!! code crash in case stdout redirected to the log file
	  // Print dump
	  gSystem->mkdir(cdbName.Data(),1);
	  TString fname=TString::Format("%s/Run%d_%d_v%d_s%d.print", cdbName.Data(), CDBId->GetFirstRun(), CDBId->GetLastRun(),
				       CDBId->GetVersion(),CDBId->GetSubVersion());
	  gROOT->ProcessLine(TString::Format("((TObject*)%p)->Dump(); >%s",cdbEntry, fname.Data()).Data());
	  gROOT->ProcessLine(TString::Format("((TObject*)%p)->Dump(); >>%s",obj, fname.Data()).Data());
	  gROOT->ProcessLine(TString::Format("((TObject*)%p)->Print(\"all\"); >>%s",obj, fname.Data()).Data());
 
	}

      }
      AliOCDBtoolkit::CleanCDBPath(cdbPath,kFALSE);
      (*pcstream)<<"ocdbTable"<<
	"run="<<fgRun<<	           // run number 
	//                         - this infromation can be stored in friend tree .. in future
	"prodInfo.="<<fgProdInfo<< // production information - can be in another tree
        "path.="<<fgPath<<         // OCDB source path
	"period.="<<&speriod<<     // period
	"pass.="<<&spass<<         // path
	"syear.="<<&syear<<        // year as string
	"year="<<year<<            // year as string
	//
	"cdbID.="<<CDBId<<       // ID of the OCDb object
	"size="<<size<<          // size of the OCDB entry
	"\n";
      fprintf(ofs,"%s\t%s\t%s/Run%d_%d_v%d_s%d.root\t%d\t%u\n",
	      cdbName.Data(),
	      cdbPath.Data(),
	      cdbName.Data(),
	      CDBId->GetFirstRun(),
	      CDBId->GetLastRun(),
	      CDBId->GetVersion(),
	      CDBId->GetSubVersion(),
	      size,
	      hash
	      );
      fflush(ofs);
      //if(!(CDBId->GetPathLevel(0)).Contains("TPC")) continue;
      //cout << CDBId.ToString() << endl;
    }
  }
  delete pcstream;
  fclose(ofs);
}


//====================================================================================================
//  Dump object part
//==================================================================================================== 





void AliOCDBtoolkit::DumpOCDBFile(const char *finput , const char *foutput, Bool_t dumpMetaData, TString  printOption ){
  ///
  ///  Function used in the AliOCDBtoolkit.sh interface ( source $ALICE_PHYSICS/PWGPP/CalibMacro/AliOCDBtoolkit.sh dumpObject <parameters>) 
  ///  DumpOCDB file in the human readable format.
  ///  
  ///  Supported options (not case sensitive)
  ///    "pocdb<x>" - TObjPrint  is used with printOption lovelX - assumed user implmenet approprate print
  ///    "docdb"  - TObject::Dump 
  ///    "XML"     - use standard root OCDB dump (not 100 % relyable sometimes failed e.g in floating poit exception)
  ///    "MI"      - custom recursive dump using the dictionary ()
  ///  Example usage:
  /*
    AliOCDBtoolkit::DumpOCDBFile("$ALICE_ROOT/OCDB/TPC/Calib/PadGainFactor/Run0_999999999_v0_s3.root", "TPC_Calib_PadGainFactor_Run0_999999999_v0_s3.print",1,"pocdb");
    AliOCDBtoolkit::DumpOCDBFile("$ALICE_ROOT/OCDB/TPC/Calib/PadGainFactor/Run0_999999999_v0_s3.root", "TPC_Calib_PadGainFactor_Run0_999999999_v0_s3.dump",1,"docdb")
    AliOCDBtoolkit::DumpOCDBFile("$ALICE_ROOT/OCDB/TPC/Calib/PadGainFactor/Run0_999999999_v0_s3.root", "TPC_Calib_PadGainFactor_Run0_999999999_v0_s3.mi",1,"mi")
    AliOCDBtoolkit::DumpOCDBFile("$ALICE_ROOT/OCDB/TPC/Calib/PadGainFactor/Run0_999999999_v0_s3.root", "TPC_Calib_PadGainFactor_Run0_999999999_v0_s3.print",1,"XML")
  */


  TString optionString = printOption;
  optionString.ToLower();
    
  if (finput==0) return ;
  if (TString(finput).Contains("alien://") && gGrid==0){
    TGrid *myGrid = TGrid::Connect("alien://");            //Oddly this will return also a pointer if connection fails
    if(myGrid->GetPort()==0){                       //if connection fails port 0 is saved, using this to check for successful connection
      cerr << "Cannot connect to grid!" << endl;
      return;
    }
  }
  TFile *falignITS  = TFile::Open(finput);
  AliCDBEntry *entry  = (AliCDBEntry*)falignITS->Get("AliCDBEntry");
  if (!entry) return; 
  TObject *obj = ((AliCDBEntry*)falignITS->Get("AliCDBEntry"))->GetObject();  
  //
  // print option indicated ()
  if (optionString.Contains("pocdb")){
    if (dumpMetaData) gROOT->ProcessLine(TString::Format("((TObject*)%p)->Dump(); >%s",entry, foutput).Data());
    if (!obj) return;
    gROOT->ProcessLine(TString::Format("((TObject*)%p)->Print(\"%s\"); >>%s",obj, printOption.Data(), foutput).Data());
    return;
  }
  if (optionString.Contains("docdb")){
    if (dumpMetaData) gROOT->ProcessLine(TString::Format("((TObject*)%p)->Dump(); >%s",entry, foutput).Data());
    if (!obj) return;
    gROOT->ProcessLine(TString::Format("((TObject*)%p)->Dump(); >>%s",obj, foutput).Data());
    return;
  }
  //
  if (optionString.Contains("mi")){
    if (dumpMetaData) gROOT->ProcessLine(TString::Format("((TObject*)%p)->Dump(); >%s",entry, foutput).Data());
    if (!obj) return;
    gROOT->ProcessLine(TString::Format("AliOCDBtoolkit::DumpObjectRecursive((TObject*)%p); >>%s",obj, foutput).Data());
    return;
  }
  if (optionString.Contains("xml")){
    TFile * f = TFile::Open(TString::Format("%s",foutput).Data(),"recreate");
    if (dumpMetaData) entry->Write("AliCDBEntry");
    else obj->Write("AliCDBEntry");
    f->Close();
    return;
  }
  ::Error("AliOCDBtoolkit::DumpOCDBFile","Not recognized option %s", optionString.Data());
}



void AliOCDBtoolkit::DumpObjectRecursive(TObject *obj){
  //
  //
  //
  Int_t counterRec=0;
  printf("==> Dumping object at: %p, name=%s, class=%s)\n", obj, obj->GetName(), (obj->IsA()->GetName()));
  DumpObjectRecursive(obj, TString(obj->IsA()->GetName())+".",counterRec);
}
 
//
//
//
void AliOCDBtoolkit::DumpObjectRecursive(TObject *obj, TString prefix, Int_t &counterRec){
  //
  // Recursive dump of the TObject
  // Dump all basic types and follow pointers to the objects
  // current limitation:
  //    a.) clases and structures not derived from TObject not followed (to fix)
  //    b.) dynamic arrays not followed
  //    c.) std maps,array ....  not followed
  //    
  //
  if (!obj) return;
  //
  // Special case of Collection classes
  //
  if (obj->IsA()->InheritsFrom(TCollection::Class())) {
    TIter myiter((TCollection*)obj);
    TObject  *arObject=0;
    Int_t counter=0;
    while ((arObject = (TObject*)myiter.Next())) {
      TString prefixArr = TString::Format("%s[%d]",prefix.Data(),counter);
      DumpObjectRecursive(arObject,prefixArr,counterRec);
      counter++;
    } 
    counterRec++;
    return;
  }

  TClass * cl = obj->IsA();
  if (!(cl->GetListOfRealData())) cl->BuildRealData();
  TRealData* rd = 0;
  TIter next(cl->GetListOfRealData());  
  while ((rd = (TRealData*) next())) {
    counterRec++;
    TDataMember* dm = rd->GetDataMember();
    TDataType* dtype = dm->GetDataType();
    Int_t offset = rd->GetThisOffset();
    char* pointer = ((char*) obj) + offset;
    
    if (dm->IsaPointer()) {
      // We have a pointer to an object or a pointer to an array of basic types.
      TClass* clobj = 0;
      if (!dm->IsBasic()) {
	clobj = TClass::GetClass(dm->GetTypeName());
      }
      if (clobj) {
	// We have a pointer to an object.
	//
	if (!clobj->InheritsFrom(TObject::Class())) {
	  // It must be a TObject object.
	  continue; 
	}
	char** apointer = (char**) pointer;
	TObject* robj = (TObject*) *apointer;
	//	
	if(!robj)
	  printf("M:%s%s\n",prefix.Data(),dm->GetName()); // Missing - 0 pointer
	else{
	  printf("T:%s\t%s%s\n", clobj->GetName(),prefix.Data(), dm->GetName());
	  TString prefixNew=prefix;
	  prefixNew+=dm->GetName();
	  prefixNew+=".";
	  if (robj!=obj) DumpObjectRecursive(robj,prefixNew,counterRec);  // trivial check 
	  if (robj==obj){
	    printf("R:%s\t%s%s\n",clobj->GetName(),prefix.Data(), dm->GetName());
	  }
	}
      }
    } else if (dm->IsBasic()) {
      //
      // Basic data type
      //
      const char* index = dm->GetArrayIndex();
      if (dm->GetArrayDim()==0){
	printf("B:\t%s%s\t%s\n", prefix.Data(),rd->GetName(), dtype->AsString(pointer));
      }
      //
      // Basic array - fixed length
      //
      //      if (dm->GetArrayDim()>0 && strlen(index) != 0){
      if (dm->GetArrayDim()>0 ){
	printf("A:\t%s%s\t",prefix.Data(),rd->GetName());
	Int_t counter=0;
	for  (Int_t idim=0; idim<dm->GetArrayDim(); idim++){
	  //printf("A:%d\t%d\n", dm->GetArrayDim(),dm->GetMaxIndex(idim));
	  for (Int_t j=0; j<dm->GetMaxIndex(idim); j++){
	    printf("%s\t",dtype->AsString(pointer+dm->GetUnitSize()*counter));
	    counter++;
	    if (counter%5==0) printf("\nA:\t%s%s\t",prefix.Data(),rd->GetName());
	  }
	}
	printf("\n");
      }
      //
      // Basic array - dynamic length
      //
      if (dm->GetArrayDim()>0 && strlen(index) != 0){
	//
	// Dump first only for the moment
	//  
	printf("B:\t%s%s\t%s\n",prefix.Data(),rd->GetName(), dtype->AsString(pointer));
      }
    } else {
    }
  }
}  

//
// Small checks to test the TRealData and TDataType
//



void DumpDataSimple(){
  //
  // Dump example for elenatr data types 
  //
  TObject *obj = new TVectorD(20);
  TClass * cl = obj->IsA();
  if (!cl->GetListOfRealData()) cl->BuildRealData();
  //
  TRealData* rd = 0;
  rd = (TRealData*)(cl->GetListOfRealData()->FindObject("fNrows"));
  TDataMember* dm = rd->GetDataMember();
  TDataType* dtype = dm->GetDataType();
  //
  Int_t offset = rd->GetThisOffset();
  char* pointer = ((char*) obj) + offset;
  printf("%s\n",dtype->AsString(pointer));
}

void DumpDataArray(){
  //
  // print array example
  // 
  TObject *obj = new TVectorD(20);
  TClass * cl = obj->IsA();
  if (!cl->GetListOfRealData()) cl->BuildRealData();
  TRealData* rd = 0;
  rd = (TRealData*)(cl->GetListOfRealData()->FindObject("*fElements"));
  TDataMember* dm = rd->GetDataMember();
  TDataType* dtype = dm->GetDataType();
  dtype->Print();
  //
  Int_t offset = rd->GetThisOffset();
  char* pointer = ((char*) obj) + offset; 
  printf("%s\n",dtype->AsString(pointer));
}

void DumpTObjectArray(){
  //
  //
  //
  TObjArray *array = new TObjArray(10);
  for (Int_t i=0; i<10; i++) array->AddLast(new TNamed(Form("n%d",i), Form("n%d",i)));  
   AliOCDBtoolkit::DumpObjectRecursive(array);
  //
  //
  TObject *obj = array;
  TClass * cl = obj->IsA();
  if (!cl->GetListOfRealData()) cl->BuildRealData();
  TRealData* rd = 0;
  rd = (TRealData*)(cl->GetListOfRealData()->FindObject("*fCont"));
  TDataMember* dm = rd->GetDataMember();
  TDataType* dtype = dm->GetDataType();
  //
  Int_t offset = rd->GetThisOffset();
  char* pointer = ((char*) obj) + offset;
  char** apointer = (char**) pointer;
  //we have pointer to pointer here
  TObject** ppobj = (TObject**) *apointer;
  (*ppobj)->Print();
  //
  TIter myiter(array);
  TObject  *arObject; 
  dtype->Print();
  while ((arObject = (TObject*)myiter.Next())) {
    AliOCDBtoolkit::DumpObjectRecursive(arObject);
  } 
}


Bool_t AliOCDBtoolkit::AddoptOCDBEntry( const char *finput, const char *output,  Int_t ustartRun, Int_t uendRun){
  //
  // Addopt OCDB entry - keeping all of the CDBentry quantities
  // // Example usage: 
  //  AliOCDBtoolkit::AddoptOCDBEntry("/cvmfs/alice.gsi.de/alice/simulation/2008/v4-15-Release/Residual/TPC/Calib/ClusterParam/Run127712_130850_v4_s0.root",0,0,AliCDBRunRange::Infinity())
  TFile * fin = TFile::Open(finput);
  if (!fin) return kFALSE;
  AliCDBEntry * entry = (AliCDBEntry*) fin->Get("AliCDBEntry");
  if (!entry) return kFALSE;
  
  AliCDBStorage* pocdbStorage = 0;
  if (output!=0) AliCDBManager::Instance()->GetStorage(output);
  else{
    TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
    pocdbStorage = AliCDBManager::Instance()->GetStorage(localStorage.Data());
  }
  //
  AliCDBId  idIn = entry->GetId();
  AliCDBMetaData *metaDataIn = entry->GetMetaData();

  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName(metaDataIn->GetObjectClassName());
  metaData->SetResponsible(TString::Format("%s: copy",metaDataIn->GetResponsible()).Data());
  metaData->SetBeamPeriod(metaDataIn->GetBeamPeriod());
  //
  metaData->SetAliRootVersion(metaDataIn->GetAliRootVersion()); //root version
  metaData->SetComment((TString::Format("%s: copy",metaDataIn->GetComment()).Data()));
  AliCDBId* id1=NULL;
  id1=new AliCDBId(idIn.GetPath(), ustartRun, uendRun);
  pocdbStorage->Put(entry->GetObject(), (*id1), metaData);
  return kTRUE;
}


void AliOCDBtoolkit::MakeSnapshotFromTxt(const TString fInput, const TString outfile, Bool_t singleKeys){
  //
  // Make snasphot form the txt file
  //
  AliCDBManager * man = AliCDBManager::Instance();
  LoadOCDBFromList(fInput.Data());
  man->DumpToSnapshotFile(outfile.Data(), singleKeys);

}


void AliOCDBtoolkit::CleanCDBPath(TString &cdbPath, Bool_t useCVMFS){
  //
  // 1.) clean alien path removing caching flags
  // 2.) repalce alien path wiht the cvmfs in case specified and cvmfs accessible
  // return value - change the content of the cdbPath string
  /* Example:

     Input:          alien://?User=?DBFolder=/alice/data/2010/OCDB?SE=default?CacheFolder=?OperateDisconnected=1?CacheSize=1073741824?CleanupInterval=0
     AliOCDBtoolkit::CleanCDBPath(cdbStorage,kTRUE);   
     AliOCDBtoolkit::CleanCDBPath Output:local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2010/OCDB

   */
  Bool_t isAlienPath = cdbPath.Contains("alien://");
  if (!isAlienPath) return;
  if (useCVMFS){
    useCVMFS=  gSystem->AccessPathName("/cvmfs/alice-ocdb.cern.ch/calibration/")==kFALSE; // is cvfms available ? bizzare function convention    
  }

  if (useCVMFS){
    TPRegexp repAlien0("alien://.*=/alice/");   // replace alien cache prefix
    repAlien0.Substitute(cdbPath, "local:///cvmfs/alice-ocdb.cern.ch/calibration/");
  }else{
    TPRegexp repAlien0("alien://.*=/");   // replace alien cache prefix
    repAlien0.Substitute(cdbPath, "alien:///");
  }
  TPRegexp repAlien1("\\?.*");        // replace ? suffix - caching flags
  repAlien1.Substitute(cdbPath, "");  
}



 
void AliOCDBtoolkit::LoadAliOCDBtoolkitSetup(Int_t verbose){
  //
  // load OCDB toolkit configuration 
  //
  TString config=gSystem->Getenv("AliOCDBtoolkitCFG");
  if (config.Length()==0){
    config="$ALICE_PHYSICS/PWGPP/CalibMacros/AliOCDBtoolkit.cfg";
  }
  if (gSystem->AccessPathName(config.Data())!=0) {
    ::Error("AliOCDBtoolkit::LoadAliOCDBtoolkitSetup","Invalid path to configuration file");
  }
  // Exclude OCDB list - Some entries has to be removed form the list as OCDB is not always back compatible
  // Exclude pattern have to start with Exclude
  TString excludeStr=gSystem->GetFromPipe(TString::Format("cat %s|grep ^Exclude: | sed s_^E.*:__  | sed 's/[[:blank:]]//g' | tr ',' '\n'", config.Data()).Data());
  fgExcludeList=excludeStr.Tokenize("\n,");
  //
  TString xmlStr=gSystem->GetFromPipe(TString::Format("cat %s|grep ^XML: | sed s_^XML.*:__  | sed 's/[[:blank:]]//g' | tr ',' '\n'", config.Data()).Data());
  fgXmlOCDBDumpList=xmlStr.Tokenize("\n,");
  //
  TString printStr=gSystem->GetFromPipe(TString::Format("cat %s|grep ^Print: | sed s_^Print.*:__  | sed 's/[[:blank:]]//g' | tr ',' '\n'", config.Data()).Data());
  fgPrintOCDBDumpList=printStr.Tokenize("\n,");    
  if (verbose>0){
    ::Info("AliOCDBtoolkit::LoadAliOCDBtoolkitSetup","ExludeList");
    fgExcludeList->Print();
    ::Info("AliOCDBtoolkit::LoadAliOCDBtoolkitSetup","XMLlist");
    fgXmlOCDBDumpList->Print();
    ::Info("AliOCDBtoolkit::LoadAliOCDBtoolkitSetup","Printlist");
    fgPrintOCDBDumpList->Print();
  }
}

Bool_t AliOCDBtoolkit::IsEntrySelected(TString entry, TObjArray *selList){
  //
  // check if the entry was selected - means it was specified in the list
  /*
    TString entry="TPC/Calib/RecoParam";
    selList=printOCDBDumpList;
  */
  //
  if (selList==NULL) return 0;
  if (selList->FindObject(entry.Data())!=NULL) return kTRUE;

  Bool_t isSelected=kFALSE;
  for (Int_t ientry=0; ientry<selList->GetEntriesFast(); ientry++){
    TRegexp regexp(selList->At(ientry)->GetName());    
    TPRegexp pregexp(selList->At(ientry)->GetName());
    if (entry.Contains(regexp)>0) isSelected=kTRUE;
    if (entry.Contains(pregexp)>0) isSelected=kTRUE;
  }
  return isSelected;
}


Double_t AliOCDBtoolkit::SetXRDTimeOutAll(Double_t timeOut ){
  //
  // Set timeout for XRD
  //   Problems oberved accessing small UserInfo objects (waiting for minutes and afterwards failing in some T2,T3)
  //   normally access aroud 2 seconds
  //   
  if (timeOut<0 && gSystem->Getenv("XRDCLIENTMAXWAIT")){
    timeOut=TString(gSystem->Getenv("XRDCLIENTMAXWAIT")).Atoi();
  }
  if (timeOut>0){
    gEnv->SetValue("XNet.RequestTimeout", timeOut);
    gEnv->SetValue("XNet.ConnectTimeout", timeOut);
    gEnv->SetValue("XNet.TransactionTimeout", timeOut);
    gEnv->SetValue("XNet.FirstConnectMaxCnt", 5);
    TFile::SetOpenTimeout(timeOut);
  }
  return timeOut;
}

