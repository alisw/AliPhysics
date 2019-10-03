/// \class AliOfflineTrigger
/// \brief This class provides fucntionality for OFFLINE Trigger raw data selection
///        and  consistency checks
///
/// Related task: https://alice.its.cern.ch/jira/browse/PWGPP-6 and PWGPP-134
/// \author Marian Ivanov   - marian.ivanov@cern.ch
/// \author Mesut Arslandok, Mikolaj  - older versions ($ALICE_PHYSICS/../src/PWGPP/rawmerge/rawmerge.C)


/*
  Example usage:
  // 0.) Init offline trigger class
  AliOfflineTrigger trigger("default", 30,100000000);
  //
  // 1.) Dump event header info. Used to make consitency checks (diff, meld, TTree diff)
  //
  trigger.DumpGIDRAWReader("alien:///alice/data/2013/LHC13b/000195483/raw/13000195483000.10.root"); 
  trigger.DumpGIDRAWTree("alien:///alice/data/2013/LHC13b/000195483/raw/13000195483000.10.root"); 
  trigger.DumpGIDESD("alien:///alice/data/2013/LHC13b/000195483/pass4/13000195483000.10/AliESDs.root","1"); 
  //
  trigger.DumpGIDESD("alien:///alice/data/2013/LHC13b/000195483/ESDs/pass3/13000195483000.10/AliESDs.root","1","gidesdpass3.list"); 
  trigger.DumpGIDESD("alien:///alice/data/2013/LHC13b/000195483/pass4/13000195483000.10/AliESDs.root","1","gidesdpass4.list"); 
  //
  //  2.) Make custom esd filter. Example  set sum pt and multiplicity trigger
  //     Define ESD trigger aliases - used later in the AliOfflineTrigger::DumpGIDESD
  //
  trigger.AddESDAlias("sumPt","Sum$(Tracks[].Pt()*(Tracks[].fTPCncls>120&&Tracks[].fITSncls>4))");  
  trigger.AddESDAlias("mult","Sum$(Tracks[].fITSncls>4&&Tracks[].fTPCncls>120&&Tracks[].fTRDncls>50)");
  trigger.AddESDAlias("triggerPt","sumPt>20");  
  trigger.AddESDAlias("triggerMult","mult>50");
  trigger.DumpGIDESD("alien:///alice/data/2013/LHC13b/000195483/pass4/13000195483000.10/AliESDs.root","(sumPt>30&&mult>50)||sumPt>50","gidesdTrigger.list"); 
  //     eventually trigger can be set in tree for cut tuning
  TFile *finput = TFile::Open("alien:///alice/data/2013/LHC13b/000195483/pass4/13000195483000.10/AliESDs.root"); 
  TTree* esdTree = (TTree*)finput->Get("esdTree");
  trigger.SetTriggerAlias(esdTree,"(sumPt>30&&mult>50)||sumPt>50");
  //
  // 4.) Make a diff between event header trees. Example compare ESD pass4, rawreader, rawtree and ESD pass3
  //     Particular case of pass3 there is an event mismatch
  //  
  tree = AliOfflineTrigger::MakeDiffTree("gidesdpass4.list","gidrawReader.list;gidrawTree.list;gidesdpass3.list")
  tree->Scan("eventID:T0.eventID:T1.eventID:T2.eventID","eventID>0"); 
  tree->Scan("timeStamp:T0.timeStamp:T1.timeStamp:T2.timeStamp","eventID>0","colsize=20"); 
  //
  // 5.) Trigger raw data file  example: 
  //
  AliOfflineTrigger trigger("default", 30,100000000);
  trigger.LoadTriggerList("gidesdTrigger.list");
  trigger.ExtractSelected("raw.list", "gidesdTrigger.list", "rawSelected[].root",1000000, 1  );
  //
  // 6.) Example to compare content of the raw filtered data with input flitered list as used in the filterning produnction
  TTree * tree = AliOfflineTrigger::MakeDiffTree("gidrawTree.list","filteredHighPt.list;filteredHighPtV0s.list;filteredMult.list");
  Int_t entriesAll=tree->GetEntries();
  Int_t entriesHighPt = tree->Draw("gid/10.","(gid==T0.gid)","goff");
  Int_t entriesHighPtV0 = tree->Draw("gid/10.","(gid==T1.gid)","goff");
  Int_t entriesMult = tree->Draw("gid/10.","(gid==T2.gid)","goff");
  ::Info("AliOfflineTrigger::MakeDiffTree.KeyValue","All:%d",entriesAll);
  ::Info("AliOfflineTrigger::MakeDiffTree.KeyValue","entriesHighPt+entriesHighPtV0+entriesMult:%d",entriesHighPt+entriesHighPtV0+entriesMult);  

*/


#include <fstream>  
#include <iostream>
#include "TError.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TEntryList.h"
#include "TPRegexp.h"
#include "TGrid.h"
#include "TEnv.h"
#include "TPad.h"
#include "TMath.h"
#include "AliRawReaderRoot.h"
#include "AliXRDPROOFtoolkit.h"
#include "AliRawVEvent.h"
#include "AliRawEventHeaderBase.h"
#include "AliOfflineTrigger.h"
#include "AliSysInfo.h"
#include "TTimeStamp.h"
#ifdef WITHALIEN
#include "TGridCollection.h"
#endif
#include "TPRegexp.h"
using std::cout;
using std::endl;

ClassImp(AliOfflineTrigger)


AliOfflineTrigger::AliOfflineTrigger(const char *triggerName, Int_t timeOut, Int_t cacheSize):
  TNamed(triggerName,triggerName),
  fTrgGIDChunkName(),                         /// GID -> ChunkName  
  fTrgGIDEventNr(),                           /// GID -> EventNumber
  fTrgGIDTimeStamp(),                         /// triger map GID -> TimeStamp map
  fCounterFileInput(0),                       ///  input file counter
  fCounterEventInput(0),                      ///  input event counter
  fCounterFileOutput(0),                      ///  output file counter
  fCounterEventOutput(0),                     ///  output event counter
  fRawName(""),                               ///  name of the output file 
  fRawTriggerFile(0),                         ///! pointer to ouput  raw trigger files
  fRawTriggerTree(0),                         ///! pointer to output raw trigger tree
  fRAWGIDChunkName(),
  fRAWGIDEventNr(),
  fRAWGIDTimeStamp(),
  fDefaultTimeOut(timeOut),
  fDefaultTreeCache(cacheSize),
  fESDTriggerList(NULL)
{
  // set timeouts
  gSystem->Setenv("XRDCLIENTMAXWAIT",Form("%d",timeOut));
  gEnv->SetValue("XNet.RequestTimeout", timeOut);
  gEnv->SetValue("XNet.ConnectTimeout", timeOut);
  gEnv->SetValue("XNet.TransactionTimeout", timeOut);
  gEnv->SetValue("XNet.FirstConnectMaxCnt", 2);  
  gEnv->SetValue("TFile.AsyncPrefetching", 1);  
}

void  AliOfflineTrigger::AddESDAlias(const char *aliasName, const char *aliasValue ){
  //
  // esd aliases to define ESD triggers
  if (fESDTriggerList==NULL) fESDTriggerList = new TObjArray(100);
  fESDTriggerList->AddLast(new TNamed(aliasName,aliasValue));
}

void AliOfflineTrigger::SetTriggerAlias(TTree* tree, const char * trigger){
  if (tree==NULL || fESDTriggerList ==NULL) return;
  Int_t entries=fESDTriggerList->GetEntries();
  for (Int_t i=0; i<entries; i++){
    tree->SetAlias(fESDTriggerList->At(i)->GetName(), fESDTriggerList->At(i)->GetTitle());
  }
  tree->SetAlias("trigger",trigger); // in trigger sets of the aliases can be used - in the list trigger mask to be returned
}


void AliOfflineTrigger::DumpGIDRAWReader(const char *rawFile){
  //
  // DUMP event ids using AliRawReaderRoot
  // Input:
  //    rawFile    - ID of raw file
  // Output: csv ascii file with following columns
  //    fname/C:eventCounter/i:gid/l:timeStamp/i:period/i:orbit/i:bcid/i:eventID/i
  /*
    Example usage:
    AliOfflineTrigger::DumpRAWGIDReader("alien:///alice/data/2013/LHC13b/000195483/raw/13000195483000.10.root"); 
  */  
  if (rawFile==NULL) return;
  if (TPRegexp("^alien").Match(rawFile) && gGrid==0)  TGrid::Connect("alien");
  TString chunkName=rawFile;
  TPRegexp reg0(".*/");
  reg0.Substitute(chunkName,"");
  TPRegexp reg1(".root$");
  reg1.Substitute(chunkName,"");

  std::ofstream file_out("gidrawReader.list");
  //  file_out<<"gid/D:eventCounter/D:period/D:orbit/D:bcid/D:fname/C:eventID/D\n";  // print header
  file_out<<"fname/C:eventCounter/i:gid/l:timeStamp/i:period/i:orbit/i:bcid/i:eventID/i"<<endl;  // print header
  Int_t counter=0;
  AliRawReaderRoot *preader = new AliRawReaderRoot(rawFile);
  if (preader==NULL){
    ::Error("dumpRAWGIDTree","rawFile %s not accessible",rawFile);
    return;
  }
  
  while (preader->NextEvent()){    
    file_out<<
      chunkName<<                        // chunkName
      "\t"<<counter<<                     // event counter
      "\t"<<preader->GetEventIdAsLong()<< // gid
      "\t"<<preader->GetTimestamp()<<     // time stamp
      "\t"<<preader->GetPeriod()<<        // period
      "\t"<<preader->GetOrbitID()<<       // orbit
      "\t"<<preader->GetBCID()<<          // bcid
      "\t"<<counter<<                     // the same as coutner for the moment. Checking the AliRawReaderRoot there is no approprate getter, I will try to add this function
      endl;                           // 
    if (counter%100==0)     printf("EventNumber\t%d\t%llu\n",counter,preader->GetEventIdAsLong());
    counter++;
  }
  file_out.close();
}

void  AliOfflineTrigger::DumpGIDRAWTree(const char *rawFile){
  // 
  // DUMP event ids using AliRawReaderRoot.  
  // Input:
  //    rawFile    - ID of raw file
  // Output: csv ascii file with following columns
  //    fname/C:eventCounter/i:gid/l:timeStamp/i:period/i:orbit/i:bcid/i:eventID/i
  
  //   In the raw data triggering  raw data are acessed using TTree functionality
  //   TTree::GetEntry() and ttree->Fill() for paricular entry numbers
  //   In some circumstancies entry number in tree and entry number if ESD can be different
  //    
  //   Dumping GID and other IDs using the tree->GetEntry()
  /* Example usage:
     dumpRAWGIDTree("alien:///alice/data/2013/LHC13b/000195483/raw/13000195483000.10.root",20,100000000); 
  */

  if (rawFile==NULL) return;
  if (TPRegexp("^alien").Match(rawFile) && gGrid==0)  TGrid::Connect("alien");
  TFile * finput=NULL;
  TTree * tree=NULL;
  TPRegexp expRoot(".root$");
  TPRegexp expList(".list$");
  TObjArray *fileList=NULL;
  
  if (TString(rawFile).Contains(expRoot)){
    fileList = new TObjArray(1);
    fileList->AddLast(new TObjString(rawFile));
  }
  if (TString(rawFile).Contains(expList)){  
    fileList = gSystem->GetFromPipe(TString::Format("cat %s", rawFile)).Tokenize("\n");
  }
  if (fileList==NULL){
    ::Error("dumpRAWGIDTree","Invalid input r %s",rawFile);
  }
  std::ofstream file_out("gidrawTree.list");
  //  file_out<<"gid/D:eventCounter/D:period/D:orbit/D:bcid/D:fname/C:eventID/D\n";  // print headaer
  //file_out<<"fname/C:eventCounter/D:gid/D:timeStamp/D:period/D:orbit/D:bcid/D:eventID/D"<<endl;  // print header
  file_out<<"fname/C:eventCounter/i:gid/l:timeStamp/i:period/i:orbit/i:bcid/i:eventID/i"<<endl;  // print header

  Int_t nFiles=  fileList->GetEntries();
  for (Int_t iFile=0; iFile<nFiles; iFile++){
    if (TPRegexp("^alien").Match(fileList->At(iFile)->GetName()) && gGrid==0)  TGrid::Connect("alien");    
    finput =TFile::Open(fileList->At(iFile)->GetName());
    if (finput==NULL){
      ::Error("dumpRAWGIDTree","rawFile %s not accessible", fileList->At(iFile)->GetName());
      continue;
    }
    tree = (TTree*)finput->Get("RAW");
    if (tree==NULL){
      ::Error("dumpRAWGIDTree","rawFile %s  does not contained requested tree",rawFile);
      continue;
    }
  
    TString chunkName=fileList->At(iFile)->GetName();
    TPRegexp yearExp("./20.*/");
    Int_t indexShort=chunkName.Index(yearExp);
    TString chunkShort(&(chunkName[indexShort+7]));
    //
    TPRegexp reg0(".*/");
    reg0.Substitute(chunkName,"");
    TPRegexp reg1(".root$");
    reg1.Substitute(chunkName,"");
    
    AliRawVEvent*    fEvent=0;              // (super) event
    tree->SetBranchStatus("TPC*",kFALSE);
    tree->SetBranchStatus("TRD*",kFALSE);
    tree->SetBranchStatus("rawevent",kTRUE);
    TBranch *fBranch = tree->GetBranch("rawevent");  // as in AliRawReaderRoot::AliRawReaderRoot  
    //    fBranch->SetAddress(&fEvent);  
    if (fDefaultTreeCache>0) {
      tree->SetCacheSize(fDefaultTreeCache);
      //tree->AddBranchToCache("*",kFALSE); 
    }
    Int_t nevents=tree->GetEntries();
    
    AliSysInfo::AddStamp(TString::Format("%s_0",chunkShort.Data()).Data(), 1, fCounterFileInput,fCounterFileInput);    
    for (Int_t ievent=0; ievent<nevents; ievent++){
      fBranch->SetAddress(&fEvent);  
      fBranch->GetEntry(ievent);
      AliRawEventHeaderBase *header = fEvent->GetHeader();
      Int_t run = header->Get("RunNb");
      //
      const UInt_t *id  = header->GetP("Id");                             // copy of AliRawReaderRoot::GetEventId()
      ULong64_t  period = id ? (((id)[0]>>4)&0x0fffffff): 0;                           // AliRawReader::Get<>
      ULong64_t  orbit  = id ? ((((id)[0]<<20)&0xf00000)|(((id)[1]>>12)&0xfffff)) : 0; // AliRawReader::Get<>
      ULong64_t  bcID   =  id ? ((id)[1]&0x00000fff) : 0;                              // AliRawReader::Get<>
      ULong64_t  gid    = (((ULong64_t)period << 36) | ((ULong64_t)orbit << 12) |(ULong64_t)bcID); // AliRawReader::GetEventIdAsLong() 
      UInt_t     timeStamp=header->Get("Timestamp");  
      delete fEvent;
      fEvent=NULL;
      //
      AliSysInfo::AddStamp(TString::Format("%s",chunkShort.Data()).Data(), 2, iFile, ievent);
      file_out<<    // dump contet
	chunkShort.Data()<<       // chunkName
	"\t"<< ievent<<           // event counter
	"\t"<<gid<<               // gid
	"\t"<<timeStamp<<         // time stamp
	"\t"<< period<<           // period
	"\t"<< orbit<<            // orbit 
	"\t"<< bcID<<             // bcID
	"\t"<<ievent<<            // event ID  - for the tree ir is the same as counter
	endl;      
      if (ievent%100==0)     printf("EventNumber\t%d\t%llu\n",ievent,gid);      
    }
    AliSysInfo::AddStamp(TString::Format("%s_1",chunkShort.Data()).Data(), 3, fCounterFileInput,fCounterFileInput);
    delete fBranch;
    delete fEvent;
    delete tree;
    finput->Close();
    delete finput;
    AliSysInfo::AddStamp(TString::Format("%s_2",chunkShort.Data()).Data(), 4, fCounterFileInput,fCounterFileInput);
  }
}


void  AliOfflineTrigger::DumpGIDESD(const char * chinput, const char *trigger,const char *choutput){
  //
  // DUMP event ids for selected events. Can be used for the fast OFFLINE triggering using ESD files only
  // 
  // Input:
  //    chinput    - esdFile
  //    trigger    - trigger filter (TCut)
  //    chouput    - name of output file
  //
  // Output: csv ascii file with following columns
  //    fname/C:eventCounter/i:gid/l:timeStamp/i:period/i:orbit/i:bcid/i:eventID/i
  /*  Example usage:
      dumpGIDESD("alien:///alice/data/2013/LHC13b/000195483/pass4/13000195483000.10/AliESDs.root");
  */
  if (TPRegexp("^alien").Match(chinput) && gGrid==0)  TGrid::Connect("alien");
  TFile * finput = TFile::Open(chinput);
  if (finput==NULL){
    ::Error("dumpRAWGIDTree","rawFile %s not accessible",chinput);
    return;
  }
  TTree * esdTree= (TTree*)finput->Get("esdTree");
  if (esdTree==NULL){
    ::Error("dumpRAWGIDTree","rawFile %s  does not contained requested tree",chinput);
    return;
  }
  TString chunkName=chinput;
  TPRegexp regE("/AliESDs.root$");
  regE.Substitute(chunkName,"");
  TPRegexp reg0(".*/");
  reg0.Substitute(chunkName,"");
  //
  if (fDefaultTreeCache>0) {
    esdTree->SetCacheSize(fDefaultTreeCache);
    // esdTree->AddBranchToCache("*",kFALSE);  // check which branches needed to evaluate trigger  formula 
    // esdTree->AddBranchToCache("Tracks*",kTRUE); 
    // esdTree->AddBranchToCache("*Header*",kTRUE); 
  }

  Int_t entries=esdTree->GetEntries();
  ::Info("dumpAndCompareGID::DumpGIDESD","Entries=%llu",esdTree->GetEntries());  
  esdTree->SetAlias("gid","((fPeriodNumber << 36) | (fOrbitNumber << 12) | fBunchCrossNumber)"); 
  SetTriggerAlias(esdTree,trigger);
  entries=esdTree->Draw("gid:Entry$:fPeriodNumber:fOrbitNumber:fBunchCrossNumber:fEventNumberInFile:fTimeStamp",trigger,"goffpara");
  //
  std::ofstream file_out(choutput);
  //  file_out<<"fname/C:eventCounter/D:gid/D:timeStamp/D:period/D:orbit/D:bcid/D:eventID/D"<<endl;  // print header
  file_out<<"fname/C:eventCounter/i:gid/l:timeStamp/i:period/i:orbit/i:bcid/i:eventID/i"<<endl;  // print header

  if (entries>0) for (Int_t i=0;i<entries; i++){    
    file_out<<
      chunkName.Data()<<                                 // chunkname
      "\t"<<Long64_t(esdTree->GetVal(1)[i])<<             // event counter
      "\t"<<Long64_t(esdTree->GetVal(0)[i])<<             // gid
      "\t"<<Long64_t(esdTree->GetVal(6)[i])<<             // timeStamp
      "\t"<<Long64_t(esdTree->GetVal(2)[i])<<             // period
      "\t"<<Long64_t(esdTree->GetVal(3)[i])<<             // orbit
      "\t"<<Long64_t(esdTree->GetVal(4)[i])<<             // bcid
      "\t"<<Long64_t(esdTree->GetVal(5)[i])<<             // eventID
      endl;
  }
  file_out.close();
}




void AliOfflineTrigger::TestDiffGIDList(){
  //
  // Routine to test internal consitency of the IDs files
  //
  // Check if the ESD gidlist and rawesd.list are identical
  // This is expected to be the case for the standard ppass reconstruction 
  //
  TTree * treeGIDESD=new TTree;
  TTree * treeGIDRAWReader=new TTree;
  TTree * treeGIDRAWTree=new TTree;
  treeGIDRAWReader->ReadFile("gidrawReeader.list","",'\t');
  treeGIDRAWTree->ReadFile("gidrawTree.list","",'\t');
  treeGIDESD->ReadFile("gidesd.list","",'\t');

  treeGIDESD->BuildIndex("gid");
  treeGIDRAWTree->BuildIndex("gid");
  treeGIDRAWReader->BuildIndex("gid");
  treeGIDESD->AddFriend(treeGIDRAWReader,"RAW");
  treeGIDESD->AddFriend(treeGIDRAWTree,"RAWTree");
  //
  if ( treeGIDRAWReader->GetEntries()!= treeGIDESD->GetEntries()){
    ::Error("dumpAndCompareGID::diffList","TestFailed. Mismatch in number of events RAW=%llu\tESD=%llu", treeGIDRAWReader->GetEntries(), treeGIDESD->GetEntries());
  }else{
    ::Info("dumpAndCompareGID::diffList","TestOK. Match in number of events RAW=%llu\tESD=%llu", treeGIDRAWReader->GetEntries(), treeGIDESD->GetEntries());     
  }
  if ( treeGIDRAWTree->GetEntries()!= treeGIDESD->GetEntries()){
    ::Error("dumpAndCompareGID::diffList","TestFailed. Mismatch in number of events RAWTree=%llu\tESD=%llu", treeGIDRAWTree->GetEntries(), treeGIDESD->GetEntries());
  }else{
    ::Info("dumpAndCompareGID::diffList","TestOK. Match in number of events RAWTree=%llu\tESD=%llu", treeGIDRAWTree->GetEntries(), treeGIDESD->GetEntries());     
  }
  Int_t entriesMatch=treeGIDESD->Draw("event-RAW.event","1");
  Double_t diff= 0;
  if (entriesMatch>0) diff= TMath::Mean(entriesMatch, treeGIDESD->GetV1());
  if (entriesMatch!= treeGIDRAWReader->GetEntries() || diff!=0){
    ::Error("dumpAndCompareGID::diffList","TestFailed. Mismatch in GID");
  }else{
    ::Info("dumpAndCompareGID::diffList","TestOK. Match in GID");     
  }   
}

TTree*    AliOfflineTrigger::MakeDiffTree(const char *refTree, const char *friendTrees){
  //
  // Build trees out of the csv files. 
  // For expert usage and debugging
  // Input: 
  //    refTree      - csv file name for main tree 
  //    friendTrees  - semicolomn separated list of the csv files defining freind trees
  // Output:
  //    tree with friend trees properly indesed attached 
  if (!friendTrees) return NULL;
  TObjArray * array=TString(friendTrees).Tokenize(";");
  Int_t ntrees=array->GetEntries();
  TTree * tree= new TTree();
  tree->ReadFile(refTree,"",'\t');
  tree->BuildIndex("gid");
  for (Int_t jf=0; jf<ntrees; jf++){
    TTree * ftree= new TTree;
    ftree->ReadFile(array->At(jf)->GetName(),"",'\t');
    ftree->BuildIndex("gid");
    tree->AddFriend(ftree,TString::Format("T%d",jf));
  }
  return tree;
}

void AliOfflineTrigger::LoadTriggerList(const char * triggerList){
  //
  // Load triger maps from the ascii trigger list csv file into internal maps
  //
  //  Assuming fname,gid, eventID  and timeStamp branch has to be present
  //           trigger branch is a recomendation  
  //
  //  Format of branches: 
  //  fname/C:eventCounter/i:gid/l:timeStamp/i:period/i:orbit/i:bcid/i:eventID/i
  /*
    const char * triggerList = "gidesd.list"
  */
  TTree * tree= new TTree();
  tree->ReadFile(triggerList,"",'\t');
  Int_t entries=tree->GetEntries();
  if (entries<=0) {
    ::Error("AliOfflineTrigger::LoadTriggerList","Invalid input trigger list\t%s", triggerList);
    return;
  }
  TBranch *chBranch=tree->GetBranch("fname");
  TBranch *gidBranch=tree->GetBranch("gid");
  TBranch *eventBranch=tree->GetBranch("eventID");
  TBranch *timeBranch=tree->GetBranch("timeStamp");
  TBranch *triggerBranch=tree->GetBranch("trigger");
  Bool_t isOK=(chBranch!=NULL && gidBranch!=NULL &&eventBranch && timeBranch!=NULL);
  if (isOK==kFALSE){
    ::Error("AliOfflineTrigger::LoadTriggerList","Missing information - in trigger list\t%s", triggerList);
    return;
  }
  ULong64_t gid;
  UInt_t eventNr,timeStamp;
  char * chbuffer= new char[10000];
  char * chtrigger= new char[10000];
  chBranch->SetAddress(chbuffer);
  gidBranch->SetAddress(&gid);
  timeBranch->SetAddress(&timeStamp);
  eventBranch->SetAddress(&eventNr);
  if (triggerBranch) triggerBranch->SetAddress(chtrigger);
  for (Int_t i=0; i<entries; i++){
    tree->GetEntry(i);
    if (triggerBranch)  fTrgGIDTrigger[gid]=chtrigger;
    fTrgGIDChunkName[gid]=chbuffer;                     
    fTrgGIDEventNr[gid]=i;                               
    fTrgGIDTimeStamp[gid]=timeStamp;                   // for consitency check
  }
  delete tree;
}



Int_t AliOfflineTrigger::LoadMapFromRawData(const char *rawFile, Int_t verbose){
  //
  //  1.) Load RAW headers and fill  event ID in std::maps
  //  2.) Calculate number of selected events
  //  3.) To be added: 
  //         hardware trigger mask 
  //         calibration event without unique gid (e.g calibration event)
  if (rawFile==NULL) return 0;
  if (TPRegexp("^alien").Match(rawFile) && gGrid==0)  TGrid::Connect("alien");  
  TFile * finput =TFile::Open(rawFile);
  if (finput==NULL){
    ::Error("dumpRAWGIDTree","rawFile %s not accessible",rawFile);
    return 0;
  }
  TTree * tree = (TTree*)finput->Get("RAW");
  if (tree==NULL){
    ::Error("dumpRAWGIDTree","rawFile %s  does not contained requested tree",rawFile);
    return 0;
  }
  TString chunkName=rawFile;
  TPRegexp reg0(".*/");
  reg0.Substitute(chunkName,"");
  TPRegexp reg1(".root$");
  reg1.Substitute(chunkName,"");
  //
  AliRawVEvent*    fEvent=0;              // (super) event
  tree->SetBranchStatus("TPC*",kFALSE);
  tree->SetBranchStatus("TRD*",kFALSE);
  tree->SetBranchStatus("rawevent",kTRUE);
  TBranch *fBranch = tree->GetBranch("rawevent");  // as in AliRawReaderRoot::AliRawReaderRoot  
  fBranch->SetAddress(&fEvent);  
  if (fDefaultTreeCache>0) {
    tree->SetCacheSize(fDefaultTreeCache);
    //tree->AddBranchToCache("*",kFALSE); 
  }
  Int_t nevents=tree->GetEntries();
  Int_t nTriggered=0;
  ULong64_t gidOld=0;
  for (Int_t ievent=0; ievent<nevents; ievent++){
    fBranch->GetEntry(ievent);
    Int_t run = fEvent->GetHeader()->Get("RunNb");
    const UInt_t *id  = fEvent->GetHeader()->GetP("Id");                             // copy of AliRawReaderRoot::GetEventId()
    ULong64_t  period = id ? (((id)[0]>>4)&0x0fffffff): 0;                           // AliRawReader::Get<>
    ULong64_t  orbit  = id ? ((((id)[0]<<20)&0xf00000)|(((id)[1]>>12)&0xfffff)) : 0; // AliRawReader::Get<>
    ULong64_t  bcID   =  id ? ((id)[1]&0x00000fff) : 0;                              // AliRawReader::Get<>
    ULong64_t  gid    = (((ULong64_t)period << 36) | ((ULong64_t)orbit << 12) |(ULong64_t)bcID); // AliRawReader::GetEventIdAsLong() 
    UInt_t     timeStamp=fEvent->GetHeader()->Get("Timestamp");  
    if (gid==gidOld && verbose>0){
      ::Error("AliOfflineTrigger::LoadMapFromRawData)","Event %d not unique GID. GID %llu can not be used as unique ID. Time %d", ievent, gid, timeStamp);
      continue;
    }
    fRAWGIDChunkName[gid]=rawFile;
    fRAWGIDTimeStamp[gid]=timeStamp;
    fRAWGIDEventNr[gid]=ievent;    
    fRAWEventNrGID[ievent]=gid;
    if (verbose==1&&ievent%100==0)  ::Info("AliOfflineTrigger::LoadMapFromRawData", "Event # \t%d\t%llu",ievent,gid);
    if ((verbose&2)>0)  ::Info("AliOfflineTrigger::LoadMapFromRawData", "Event # \t%d\t%llu",ievent,gid);
    // check if event triggered
    UInt_t trgTimeStamp=fTrgGIDTimeStamp[gid];
    if (trgTimeStamp>0){
      nTriggered++;
      ::Info("AliOfflineTrigger::LoadMapFromRawData()","Trig.event=%d. GID=%llu\tTime=%d", ievent, gid, timeStamp);
      if (trgTimeStamp!=timeStamp){
	::Error("AliOfflineTrigger::LoadMapFromRawData()","Inconsistent timestamp");
      }else{
	fTrgGIDEventNr[gid]=ievent;
      }
    }
    gidOld=gid;
  }
  tree->SetCacheSize(0);
  delete tree; 
  delete finput;
  return nTriggered;
}



void  AliOfflineTrigger::ExtractSelected(const char *rawFile, Int_t verbose){
  //
  // Extract selected events from raw data file
  //    Trigger list has to be loaded before filtering
  //
  AliRawVEvent*    fEvent=0;              // event
  TBranch *fEventBranch = 0;              // branch for event header   

  if (fRawTriggerFile==NULL){
    ::Error("AliOfflineTrigger::ExtractSelected","Output file not intitialized");
    return;
  }
  Int_t nTriggered=0;
  if (fTrgGIDTimeStamp.size()<=0){
    ::Info("AliOfflineTrigger::LoadMapFromRawData", "Trigger list is empty. Nothing to filter for file %s", rawFile);
    return;
  }
  nTriggered=LoadMapFromRawData(rawFile, verbose);
  ::Info(" AliOfflineTrigger::ExtractSelected", "%s:\t N_{trig}=%d", rawFile,nTriggered);
  if ( nTriggered==0){
    return;
  }
  TFile * finput = TFile::Open(rawFile);
  if (finput==NULL){
    ::Info(" AliOfflineTrigger::ExtractSelected", "File %s not exist or not accessible within TimeOut %d", rawFile, fDefaultTimeOut);
    return;
  }
  if (finput->IsOpen()==kFALSE || finput->IsZombie()){
    ::Info(" AliOfflineTrigger::ExtractSelected", "File %s not accessible within TimeOut %d", rawFile, fDefaultTimeOut);
    return;
  }
  TTree *itree=dynamic_cast<TTree*>(finput->Get("RAW"));
  if (finput==NULL){
    ::Info(" AliOfflineTrigger::ExtractSelected", "Tree in file %s not exist or not accessible within TimeOut %d", rawFile, fDefaultTimeOut);
  }
  
  fEventBranch = itree->GetBranch("rawevent");  // as in AliRawReaderRoot::AliRawReaderRoot  
  fEventBranch->SetAddress(&fEvent);           // access event header


  if (fRawTriggerTree==NULL){
    fRawTriggerFile->cd();
    fRawTriggerTree=itree->CloneTree(0);
  }
  Int_t nEvents=itree->GetEntries();
  Int_t size= itree->GetEntry(0);
  fRawTriggerTree->CopyAddresses(itree);
  for (Int_t iEvent=0; iEvent<nEvents; iEvent++){
    ULong64_t gid=fRAWEventNrGID[iEvent];
    Bool_t isSelected = kFALSE;
    isSelected=(fTrgGIDEventNr[gid]==iEvent && fTrgGIDTimeStamp[gid]>0);   // gid in the list of triggered events
    fCounterFileInput++;
    if (!isSelected) continue;
    if ((verbose&1)>0) {
      ::Info(" AliOfflineTrigger::ExtractSelected", "%s\t%d\t%llu\t%d\t%s",rawFile,iEvent,gid,fTrgGIDTimeStamp[fRAWEventNrGID[iEvent]], TTimeStamp(fTrgGIDTimeStamp[fRAWEventNrGID[iEvent]]).AsString("short"));
    }
    AliSysInfo::AddStamp(TString::Format("%s_BR",rawFile).Data(), 10, fCounterFileInput,fCounterFileInput);    
    size= itree->GetEntry(iEvent);
    AliSysInfo::AddStamp(TString::Format("%s_ER",rawFile).Data(), 11, fCounterFileInput,fCounterFileInput);    
    Int_t readEntry=itree->GetReadEntry();   
    //fRawTriggerTree->CopyAddresses(itree);
    AliSysInfo::AddStamp(TString::Format("%s_BF",rawFile).Data(), 100, fCounterFileInput,fCounterFileInput);
    fRawTriggerTree->Fill();
    AliSysInfo::AddStamp(TString::Format("%s_EF",rawFile).Data(), 101, fCounterFileInput,fCounterFileInput);
    fCounterFileOutput++;    
  }
  itree->ResetBranchAddresses();
  delete fEvent;
  delete itree;
  delete finput;
  fCounterFileInput++;
}


 
void   AliOfflineTrigger::ExtractSelected(const char *rawList, const char * triggerList, const char * outputName, Long_t maxSize,  Int_t verbose){
  //
  /*
   const char *rawList="raw.list";
   const char *triggerList="gidesdTrigger.list";
   const char *outputName="filteredraw_[].root"
   Int_t maxSize=100000000;  // max size of raw trees before moving to next file
   verbose=2
  */
  TObjArray* rawArray = 0;
  if (TPRegexp(".xml$").Match(rawList)){
#ifdef WITHALIEN
    TGridCollection *coll = gGrid->OpenCollection(rawList);
    Int_t nFiles =  coll->GetNofGroups();
    rawArray=new TObjArray(nFiles);
    while( coll->Next()){
      rawArray->AddLast(new TObjString(coll->GetTURL()));
    }
#endif
  }else{
    rawArray = gSystem->GetFromPipe(TString::Format("cat %s", rawList)).Tokenize("\n");  
  }
  if (rawArray == 0) {
    ::Fatal("AliOfflineTrigger::ExtractSelected","Unable to get list of RAW files");
  }
  Int_t nFiles=rawArray->GetEntries();
  if (nFiles<=0){
    ::Error("AliOfflineTrigger::ExtractSelected","Empty input list");
  }
  LoadTriggerList(triggerList);
  for (Int_t iFile=0; iFile<nFiles; iFile++){
    if (fRawTriggerTree!=NULL){
      if (fRawTriggerTree->GetZipBytes()>maxSize){ // close file if content bigger than max Size
	fRawTriggerFile->cd();
	fRawTriggerTree->Write();
	fRawTriggerTree->ResetBranchAddresses();
	delete fRawTriggerTree;
	fRawTriggerFile->Close();
	fRawTriggerFile=NULL;
	fRawTriggerTree=NULL;
	fCounterFileOutput++;
      }
    }
    if (fRawTriggerFile==NULL){
      TString fName=outputName;
      fName.ReplaceAll("[]",TString::Format("%d",fCounterFileOutput));
      fRawTriggerFile= TFile::Open(fName.Data(),"recreate");      
    }
    ExtractSelected(rawArray->At(iFile)->GetName(),verbose); 
  }
  if (fRawTriggerTree!=NULL){  
    fRawTriggerFile->cd();
    fRawTriggerTree->Write();
    fRawTriggerTree->ResetBranchAddresses();
    fRawTriggerFile->Close();
    //delete fRawTriggerTree;
    delete fRawTriggerFile;
    fRawTriggerFile=NULL;
    fRawTriggerTree=NULL;
    fCounterFileOutput++;
  }
}
