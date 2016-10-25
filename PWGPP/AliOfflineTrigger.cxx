/// \class AliOfflineTrigger
/// \brief This class provides fucntionality for OFFLINE Trigger raw data selection
///        and  conistency checks
///
/// Related task: https://alice.its.cern.ch/jira/browse/PWGPP-6 and PWGPP-134
/// \author Marian Ivanov   - marian.ivanov@cern.ch
/// \author Mesut Arslandok, Mikolaj  - older versions (rawmerege.C)


/*
  Example usage:
  //
  // 0.) Init offline trigger class
  //
  .L $ALICE_PHYSICS/../src/PWGPP/AliOfflineTrigger.cxx+
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
  //  2.) Make custom esd filer. Example  set sum pt and multiplicity trigger
  //     Define ESD trigger aliases - used later in the AliOfflineTrigger::DumpGIDESD
  //
  trigger.AddESDAlias("sumPt","Sum$(Tracks[].Pt()*(Tracks[].fTPCncls>120&&Tracks[].fITSncls>4))");  
  trigger.AddESDAlias("mult","Sum$(Tracks[].fITSncls>4&&Tracks[].fTPCncls>120&&Tracks[].fTRDncls>50)");
  trigger.AddESDAlias("triggerPt","sumPt>20");  
  trigger.AddESDAlias("triggerMult","mult>50");
  trigger.DumpGIDESD("alien:///alice/data/2013/LHC13b/000195483/pass4/13000195483000.10/AliESDs.root","(sumPt>30&&mult>50)||sumPt>50"); 
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
using std::cout;
using std::endl;

ClassImp(AliOfflineTrigger)


AliOfflineTrigger::AliOfflineTrigger(const char *triggerName, Int_t timeOut, Int_t cacheSize):
  TNamed(triggerName,triggerName),
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
  //
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
  //    timeOut    - to limit access time (in case fraction of raw file not staged on alien or temporary not available)
  //    cacheSize  - ttree cache size to speedup reading
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

  ofstream file_out("gidrawReader.list");
  //  file_out<<"gid/D:eventCounter/D:period/D:orbit/D:bcid/D:fname/C:eventID/D\n";  // print header
  file_out<<"fname/C:eventCounter/D:gid/D:timeStamp/D:period/D:orbit/D:bcid/D:eventID/D"<<endl;  // print header
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
  //   In the raw dat triggering  we are accessing  raw data only using TTree functionality
  //   TTree::GetEntry() and ttree->Fill() for paricular entry numbers
  //   In some ciscumstancies entry number in tree and entry number if ESD can be different
  //    
  //   Dumping GID and other IDs using the tree->GetEntry()
  //   Not this code is a hack to get map: entry number <--> gid
  /* Example usage:
     dumpRAWGIDTree("alien:///alice/data/2013/LHC13b/000195483/raw/13000195483000.10.root",20,100000000); 
  */
  //   The same varaibles and content can be crosschecked with routine using AliRawReaderRoot

  if (rawFile==NULL) return;
  if (TPRegexp("^alien").Match(rawFile) && gGrid==0)  TGrid::Connect("alien");

  TFile * finput =TFile::Open(rawFile);
  if (finput==NULL){
    ::Error("dumpRAWGIDTree","rawFile %s not accessible",rawFile);
    return;
  }
  TTree * tree = (TTree*)finput->Get("RAW");
  if (tree==NULL){
    ::Error("dumpRAWGIDTree","rawFile %s  does not contained requested tree",rawFile);
    return;
  }
  TString chunkName=rawFile;
  TPRegexp reg0(".*/");
  reg0.Substitute(chunkName,"");
  TPRegexp reg1(".root$");
  reg1.Substitute(chunkName,"");

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

  ofstream file_out("gidrawTree.list");
  //  file_out<<"gid/D:eventCounter/D:period/D:orbit/D:bcid/D:fname/C:eventID/D\n";  // print headaer
  file_out<<"fname/C:eventCounter/D:gid/D:timeStamp/D:period/D:orbit/D:bcid/D:eventID/D"<<endl;  // print header
  
  for (Int_t ievent=0; ievent<nevents; ievent++){
    fBranch->GetEntry(ievent);
    Int_t run = fEvent->GetHeader()->Get("RunNb");
    //
    const UInt_t *id  = fEvent->GetHeader()->GetP("Id");                             // copy of AliRawReaderRoot::GetEventId()
    ULong64_t  period = id ? (((id)[0]>>4)&0x0fffffff): 0;                           // AliRawReader::Get<>
    ULong64_t  orbit  = id ? ((((id)[0]<<20)&0xf00000)|(((id)[1]>>12)&0xfffff)) : 0; // AliRawReader::Get<>
    ULong64_t  bcID   =  id ? ((id)[1]&0x00000fff) : 0;                              // AliRawReader::Get<>
    ULong64_t  gid    = (((ULong64_t)period << 36) | ((ULong64_t)orbit << 12) |(ULong64_t)bcID); // AliRawReader::GetEventIdAsLong() 
    UInt_t     timeStamp=fEvent->GetHeader()->Get("Timestamp");  
    //
    file_out<<    // dump contet
      chunkName.Data()<<       // chunkName
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
}


void  AliOfflineTrigger::DumpGIDESD(const char * chinput, const char *trigger,const char *choutput){
  //
  //  Dump Event ID for selected ESD events
  //    
  //
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
  ofstream file_out(choutput);
  file_out<<"fname/C:eventCounter/D:gid/D:timeStamp/D:period/D:orbit/D:bcid/D:eventID/D"<<endl;  // print header

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
  // Check if the ESD gidlist and rawesd.list are identical
  // This is expected to be the case for the standard ppass reconstruction 
  //
  TTree * treeGIDESD=new TTree;
  TTree * treeGIDRAWReader=new TTree;
  TTree * treeGIDRAWTree=new TTree;
  treeGIDRAWReader->ReadFile("gidrawReeader.list");
  treeGIDRAWTree->ReadFile("gidrawTree.list");
  treeGIDESD->ReadFile("gidesd.list");

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
  // 
  //
  if (!friendTrees) return NULL;
  TObjArray * array=TString(friendTrees).Tokenize(";");
  Int_t ntrees=array->GetEntries();
  TTree * tree= new TTree();
  tree->ReadFile(refTree);
  tree->BuildIndex("gid");
  for (Int_t jf=0; jf<ntrees; jf++){
    TTree * ftree= new TTree;
    ftree->ReadFile(array->At(jf)->GetName());
    ftree->BuildIndex("gid");
    tree->AddFriend(ftree,TString::Format("T%d",jf));
  }
  return tree;
}
