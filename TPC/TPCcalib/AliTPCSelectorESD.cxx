// The class definition in esdClus.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("esdClus.C")
// Root > T->Process("esdClus.C","some options")
// Root > T->Process("esdClus.C+")
//
// Modification log:
// 05/11/2006 HH  Correct for large pads (outer sectors) in amplitude plots

#include "TSystem.h"
#include <TPDGCode.h>
#include <TStyle.h>
#include "TCint.h"
#include "TH1I.h"
#include "TTimeStamp.h"
#include "TProof.h"
#include "TTree.h"
//
#include "AliTracker.h"
#include "AliMagF.h"
// 
#include "AliESDEvent.h"   // new container
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
//
#include "AliSysInfo.h"
#include "AliTPCSelectorESD.h" 


 


AliTPCSelectorESD::AliTPCSelectorESD(TTree *tree) : 
   TSelector(),
   fChain(0),
   fESDevent(0),
   fESD(0),
   fESDfriend(0),
   fFileNo(0),
   fSysWatch(0),      // system info        
   fFileWatch(0),      // file watch - write the status of the analyzed files
   fDebugLevel(0)
 {
   G__SetCatchException(0);     
   //
   //
   //if (somthing){  
   fSysWatch  = new fstream("syswatch2.log", ios_base::out|ios_base::trunc);
   fFileWatch = new fstream("filewatch.log", ios_base::out|ios_base::trunc);
   if (gProof) fDebugLevel = gProof->GetLogLevel();
   if (tree) fChain=tree;
 }   


void AliTPCSelectorESD::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();

}

void AliTPCSelectorESD::SlaveBegin(TTree * tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
  if (tree) fChain = tree;
  Init(tree);
  //
  fNtracks       = new TH1I("ntracks","Number of tracks",100,0,400);
  fNtracksFriend = new TH1I("ntracksF","Number of friend tracks",100,0,400);
  fNClusters     = new TH1I("ncluster","Number of clusters",100,0,200);
  fOutput->AddLast(fNtracks);
  fOutput->AddLast(fNtracksFriend);
  fOutput->AddLast(fNClusters);
  

}

void   AliTPCSelectorESD::CleanESD(){
  //
  Bool_t isNew =  fESDevent!=0;
  if (isNew) return;
  

  if (fESD!=0){
    delete fESD;
    fESD = 0;
  }
  if (fESDevent!=0){
    delete fESDevent;
    fESDevent=0;
  }
  if (fESDfriend){
    delete fESDfriend;
    fESDfriend =0;
  }
}

Bool_t AliTPCSelectorESD::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either AliTPCSelectorESD::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  Int_t status = ProcessIn(entry);
  if (fFileWatch) {
    (*fFileWatch) << "__" << status ;
  }
  return (status==0);
}




Int_t   AliTPCSelectorESD::ReadEvent(Long64_t entry){
  //
  //
  //


  if (!fChain) return -1;  
  if (!fChain->GetTree()) return -1; 
  try {
    fChain->GetTree()->GetEntry(entry);
  } catch (std::bad_alloc) {
    printf("Pica vyjebana pojebany skurveny kokot piciak\n");
    //    fESD =0;
    //fESDfriend = 0;
    //fESDevent=0;
    return -1;
  }
  //
  // Info("Procces","0");
  if (!fESD && !fESDevent) { 
    //fESD = 0;
    //fESDfriend = 0;
    //CleanESD();
    return -2;
  }
  Int_t ntracks = (fESD) ? fESD->GetNumberOfTracks() : fESDevent->GetNumberOfTracks();   

  if (fESDevent){
    fESDevent->SetESDfriend(fESDfriend);   
  }


  fNtracks->Fill(ntracks);
  Info("Procces", Form("entry\t%d: Ntracks = %d",entry, ntracks));
  
  if (!fESDfriend || fESDfriend->GetNumberOfTracks() != ntracks) {
    try {
      delete fESD;
    }
    catch (std::bad_alloc) {
      printf("Pica vyjebana pojebany skurveny kokot piciak\n");
      fESD =0;
      return -1;
    }
    //fESD = 0;
    //fESDfriend = 0;
    //fESD = 0;
    //    CleanESD(); 
    if (fESDfriend) fNtracksFriend->Fill(fESDfriend->GetNumberOfTracks());
    Info("Procces","2: PROBLEM");
    return -3;
  }
  if (fESD) fESD->SetESDfriend(fESDfriend);
  return 0;
}


Int_t AliTPCSelectorESD::ProcessIn(Long64_t entry)
{
  //
  // User part of proccess method
  //
 
  //
  // USER code to go here
  //
  Int_t status = ReadEvent(entry);
  if (status<0) return status; 
  Int_t ntracks = (fESD) ? fESD->GetNumberOfTracks() : fESDevent->GetNumberOfTracks();   
  //
  AliTPCseed *seed;
  AliTPCclusterMI cl;
  Info("Procces", Form("entry\t%d: NtracksF = %d",entry,fESDfriend->GetNumberOfTracks() ));

  for (Int_t tr = 0; tr < ntracks; tr++){ 
    AliESDtrack *esdTrack = fESD ? (AliESDtrack*) fESD->GetTrack(tr): (AliESDtrack*) fESDevent->GetTrack(tr);
    AliESDfriendTrack *friendtrack = (AliESDfriendTrack*) esdTrack->GetFriendTrack();
    seed = 0;
    TObject *cobject=0;
    for (Int_t i=0;;i++){
      cobject = friendtrack->GetCalibObject(i);
      if (!cobject) break;
      seed = dynamic_cast<AliTPCseed*>(cobject);
      if (seed) break;
    }
    //
    //Info("Process",Form("Proccessing track%d\n",tr));
    if (seed) { 
      fNClusters->Fill(seed->GetNumberOfClusters());
      //
      //
    }
  }
  CleanESD();
  return 0;
  
}


void AliTPCSelectorESD::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
    printf ("SlaveTerminate.. \n");
    
}

void AliTPCSelectorESD::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
  printf ("Terminate... \n");
  if (!fOutput) return;
  TFile file("Output.root","recreate");
  printf("fOutput contains the following: \n");
  fOutput->Print();
  printf("Trying to write the file 'Output.root'... \n");
  fOutput->Write();
  file.Close();
  
  
}



void AliTPCSelectorESD::Init(TTree *tree) 
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
  static Int_t counter=0;
  printf(Form("\nAliTPCSelectorESD::Init Accesing%d time\n",counter));
  counter++;
  if (!tree) return;
  fChain = tree;
  //if (counter>1) return;
  tree->SetBranchStatus("*",1);
  //
  // New AliESDevent format
  //
  if (!fChain->GetBranch("ESD")){
    //
    //
    //
    if (fESDevent) delete fESDevent;
     fESDevent = new AliESDEvent();
     fESDevent->ReadFromTree(tree); // Attach the branch with ESD friends
     fESDfriend = (AliESDfriend*)fESDevent->FindListObject("AliESDfriend");
     tree->SetBranchAddress("ESDfriend.",&fESDfriend); 
     return;
  }
  //
  // if old format
  //
  


   //   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("ESD",&fESD);
   Info("Init","Enter");
   Bool_t isOK=kFALSE;
   if (fChain->GetBranch("ESDfriend")) {
     fChain->SetBranchAddress("ESDfriend",&fESDfriend);
     Info("Init","V0-ESDfriend.");
     isOK=kTRUE;
   }
   if (fChain->GetBranch("ESDfriend.")){
     Info("Init","V1-ESDfriend.");
     fChain->SetBranchAddress("ESDfriend.",&fESDfriend);
     isOK=kTRUE;
   }
   if (isOK) return;

   //
   // Try to solve problem
   //

   Info("Init","Problem");
   if (tree->GetBranch("ESD")){
     Info("InitTree",tree->GetBranch("ESD")->GetFile()->GetName());
     char  fname[1000];
     sprintf(fname,"%s/AliESDfriends.root",gSystem->DirName(tree->GetBranch("ESD")->GetFile()->GetName()));
     Info("InitFile",fname);
     if (tree->AddFriend("esdFriendTree",fname)){
       Info("InitFileOK",fname);
       if (fChain->GetBranch("ESDfriend")) {
	 fChain->SetBranchAddress("ESDfriend",&fESDfriend);
	 Info("Init","V0-ESDfriend.");
	 isOK=kTRUE;
       }
       if (fChain->GetBranch("ESDfriend.")){
	 Info("Init","V1-ESDfriend.");
	 fChain->SetBranchAddress("ESDfriend.",&fESDfriend);
	 isOK=kTRUE;
       }       
     }   
   }
}



Bool_t AliTPCSelectorESD::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  ++fFileNo;
  const char * fname = "UNKNOWN";
  const char * hname = gSystem->HostName();
  if (!fChain) return kFALSE;
  if (fChain->GetCurrentFile()){
    fname = fChain->GetCurrentFile()->GetName();
  }
  Info("Notify",Form("Host %s processing file no %d %s\n",hname,fFileNo,fname));

  //
  // Print statistic to log file
  //
//   if (fname) {
//     (*fFileWatch) << endl;
//     (*fFileWatch) << hname   <<"\t"
// 		  << fname   <<"\t";
//   }
//   DumpSysInfo(-1);
  
  return kTRUE;
}


void   AliTPCSelectorESD::DumpSysInfo(Int_t entry){
  //
  // dump system info to log file
  // entry  - entry number in the chain
  //
  const char * fname = "UNKNOWN";
  const char * hname = gSystem->HostName();
  if (fChain->GetCurrentFile()){
    fname = fChain->GetCurrentFile()->GetName();
  }
  //   //
  if (fSysWatch){
    TTimeStamp stamp;
    CpuInfo_t  cpuInfo;
    MemInfo_t  memInfo;
    ProcInfo_t procInfo;
    
    gSystem->GetCpuInfo(&cpuInfo, 1000);
    gSystem->GetMemInfo(&memInfo);
    gSystem->GetProcInfo(&procInfo);
    
    (*fSysWatch) << hname   <<"\t"               // hname - hostname
		 << fname   <<"\t"               // fname - filename
		 << entry   <<"\t"               // entry - entry number
		 << stamp.GetSec()<<"\t"         // time  - time stamp in seconds
 		 << memInfo.fMemUsed<<"\t"       //  
		 << memInfo.fSwapUsed<<"\t"      //
		 << procInfo.fMemResident<<"\t"  //
		 << procInfo.fMemVirtual<<"\t"   //    
		 << cpuInfo.fUser <<"\t"         //
		 << cpuInfo.fSys  <<"\t"         //
		 << procInfo.fCpuUser<<"\t"      //
		 << procInfo.fCpuSys<<"\t"       //
		 << endl;
  }
  AliSysInfo::AddStamp(fname);
}
