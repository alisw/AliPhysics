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
#include "TH2F.h"
#include "TH1F.h"

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
#include "AliComparisonDraw.h" 
#include "AliComparisonSelector.h" 




AliComparisonSelector::AliComparisonSelector(TTree *) : 
   TSelector(),
   //
   fChain(0),
   fFileNo(0),
   fSysWatch(0),      // system info        
   fFileWatch(0),      // file watch - write the status of the analyzed files
   fDebugLevel(0),
   fInfoMC(0),
   fInfoRC(0),
   fComp(0)
{
  G__SetCatchException(0);     
  //
  //
  //if (somthing){  
  fSysWatch  = new fstream("syswatch.log", ios_base::out|ios_base::trunc);
  fFileWatch = new fstream("filewatch.log", ios_base::out|ios_base::trunc);
  if (gProof) fDebugLevel = gProof->GetLogLevel(); 
}   


void AliComparisonSelector::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();

}

void AliComparisonSelector::SlaveBegin(TTree * tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
  fChain = tree;
  fComp = new AliComparisonDraw;
  Init(tree);
  fOutput->AddLast(fComp);
  //
}

void   AliComparisonSelector::Clean(){
  //
}

Bool_t AliComparisonSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either AliComparisonSelector::GetEntry() or TBranch::GetEntry()
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




Int_t   AliComparisonSelector::ReadEvent(Long64_t entry){
  //
  //
  //


  if (!fChain) return -1;  
  if (!fChain->GetTree()) return -1;  
  TTree * tree = fChain->GetTree();
  if (tree->GetBranch("MC") &&  tree->GetBranch("RC")){
     tree->GetBranch("MC")->SetAddress(&fInfoMC);
     tree->GetBranch("RC")->SetAddress(&fInfoRC);
  } 
  try {
    fChain->GetTree()->GetEntry(entry);
  } catch (std::bad_alloc) {
    printf("Pica vyjebana pojebany skurveny kokot piciak\n");
    return -1;
  }
  return 0;
}


Int_t AliComparisonSelector::ProcessIn(Long64_t entry)
{
  //
  // User part of proccess method
  //
 
  //
  // USER code to go here
  //
  Int_t status = ReadEvent(entry);
  if (status<0) return status; 
  fComp->Process(fInfoMC, fInfoRC);
  Clean();
  return 0;  
}


void AliComparisonSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
    printf ("SlaveTerminate.. \n");
    
}

void AliComparisonSelector::Terminate()
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



void AliComparisonSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers

   if (!tree) return;
   fChain = tree;
   tree->SetBranchStatus("*",1);
   
   if (tree->GetBranch("MC") &&  tree->GetBranch("RC")){
     tree->GetBranch("MC")->SetAddress(&fInfoMC);
     tree->GetBranch("RC")->SetAddress(&fInfoRC);
   }   
}



Bool_t AliComparisonSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  ++fFileNo;
  const char * fname = "UNKNOWN";
  const char * hname = gSystem->HostName();
  if (fChain->GetCurrentFile()){
    fname = fChain->GetCurrentFile()->GetName();
  }


  //
  Info("Notify",Form("Host %s processing file no %d %s\n",hname,fFileNo,fname));

  //
  // Print statistic to log file
  //
  if (fname) {
    (*fFileWatch) << endl;
    (*fFileWatch) << hname   <<"\t"
		  << fname   <<"\t";
  }
  DumpSysInfo(-1);
  
  return kTRUE;
}

 
void   AliComparisonSelector::DumpSysInfo(Int_t entry){
  //
  // dump system info to log file
  // entry  - entry number in the chain
  //
  const char * fname = "UNKNOWN";
  const char * hname = gSystem->HostName();
  if (fChain->GetCurrentFile()){
    fname = fChain->GetCurrentFile()->GetName();
  }
  //
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
}




AliComparisonDraw::AliComparisonDraw():
  TObject(),
  fPtResolLPT(0),
  fPtResolHPT(0)
{
  InitHisto();
}

void AliComparisonDraw::InitHisto(){
  //
  //
  //
  fPtResolLPT = new TH2F("Pt resol","pt resol",10, 0.1,3,200,-0.2,0.2);
  fPtResolHPT = new TH2F("Pt resol","pt resol",10, 2,100,200,-0.3,0.3);  
  //
  fPtPoolLPT = new TH2F("Pt pool","pt pool",10, 0.1,3,200,-6,6);
  fPtPoolHPT = new TH2F("Pt pool","pt pool",10, 2,100,200,-6,6);  
}

void AliComparisonDraw::Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC){
  //
  // 
  //
  Float_t mcpt = infoMC->GetParticle().Pt();

  //
  //
  if (infoRC->GetStatus(1)==0) return;
  if (!infoRC->GetESDtrack()) return;  //buggy line


  if (infoRC->GetESDtrack()->GetTPCNcls()<10) return;

  //  printf("Pt\t%f\t%f\n",mcpt, infoRC->GetESDtrack()->Pt());
  
  Float_t deltaPt= (mcpt-infoRC->GetESDtrack()->Pt())/mcpt;  
  Float_t poolPt= (1/mcpt-infoRC->GetESDtrack()->OneOverPt())/
    TMath::Sqrt(infoRC->GetESDtrack()->GetSigma1Pt2());  
  fPtResolLPT->Fill(mcpt,deltaPt);
  fPtResolHPT->Fill(mcpt,deltaPt);
  fPtPoolLPT->Fill(mcpt,poolPt);
  fPtPoolHPT->Fill(mcpt,poolPt); 
}




