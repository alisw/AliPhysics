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
//
#include "AliMagF.h"
#include "AliTracker.h"
//
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
#include "AliClusterMap.h"
//
#include "AliTPCcalibTracks.h"
//#include "AliTPCcalibTracks.cxx"
#include "AliTPCSelectorTracks.h" 





AliTPCSelectorTracks::AliTPCSelectorTracks(TTree *) : 
   TSelector(),
   fChain(0),
   fESD(0),
   fESDfriend(0),
   fFileNo(0)     
 {
   G__SetCatchException(0);     
 }   

void AliTPCSelectorTracks::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();

}

void AliTPCSelectorTracks::SlaveBegin(TTree * tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
  fChain = tree;
  Init(tree);
  //
  fNtracks       = new TH1I("ntracks","Number of tracks",100,0,200);
  fNtracksFriend = new TH1I("ntracksF","Number of friend tracks",100,0,200);
  fNClusters     = new TH1I("ncluster","Number of clusters",100,0,200);
  fOutput->AddLast(fNtracks);
  fOutput->AddLast(fNtracksFriend);
  fOutput->AddLast(fNClusters);
  fCalibTracks = new AliTPCcalibTracks;
    //
  fCalibTracks->ProofSlaveBegin(fOutput);


}

void   AliTPCSelectorTracks::CleanESD(){
  //
  if (fESD!=0){
    delete fESD;
    fESD = 0;
  }
  if (fESDfriend){
    delete fESDfriend;
    fESDfriend =0;
  }
}


Bool_t AliTPCSelectorTracks::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either AliTPCSelectorTracks::GetEntry() or TBranch::GetEntry()
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
  
  if (!fChain) return kFALSE;  
  if (!fChain->GetTree()) return kFALSE; 
  try {
    fChain->GetTree()->GetEntry(entry);
  } catch (std::bad_alloc) {
    printf("Pica vyjebana pojebany skurveny kokot piciak\n");
    fESD =0;
    fESDfriend = 0;
    return 0;
  }
  //
  Info("Procces","0");
  if (!fESD) { 
    fESD =0;
    fESDfriend=0;
    //CleanESD();
    return kFALSE;
  }
  Int_t ntracks = fESD->GetNumberOfTracks();   

  fNtracks->Fill(ntracks);
  Info("Procces",Form("1-Ntracks=%d",ntracks));
  
  if (!fESDfriend || fESDfriend->GetNumberOfTracks()!=ntracks) {
    fESD =0;
    fESDfriend=0;
    //    CleanESD(); 
    if (fESDfriend) fNtracksFriend->Fill(fESDfriend->GetNumberOfTracks());
    Info("Procces","2- PROBLEM");
    return kFALSE;
  }
  fESD->SetESDfriend(fESDfriend);
  //
  // USER code to go here
  //
  AliTPCseed *seed;
  AliTPCclusterMI cl;
  
  for (Int_t tr=0; tr < ntracks; tr++){ 
    AliESDtrack * esdTrack = (AliESDtrack*) fESD->GetTrack(tr);
    AliESDfriendTrack *friendtrack = (AliESDfriendTrack*) fESD->GetTrack(tr)->GetFriendTrack();
    seed = (AliTPCseed*)(friendtrack->GetCalibObject(0));
    if (seed) { 
      if (!fCalibTracks->AcceptTrack(seed)) continue;
      //      FillHistoCluster(seed);
      fCalibTracks->FillResolutionHistoLocal(seed);
      fCalibTracks->AlignUpDown(seed,esdTrack);
      fNClusters->Fill(seed->GetNumberOfClusters());
      //
      //
    }
  }
  CleanESD();
  return kTRUE;
  
}


void AliTPCSelectorTracks::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
    printf ("SlaveTerminate.. \n");
    
}

void AliTPCSelectorTracks::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  if (!fOutput) return;
  TFile file("Output.root","recreate");
  fOutput->Write();
}
void AliTPCSelectorTracks::Init(TTree *tree)
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

Bool_t AliTPCSelectorTracks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

    ++fFileNo;
    printf ("Processing file no %d\n",fFileNo);
 
   return kTRUE;
}


