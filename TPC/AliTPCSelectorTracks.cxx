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
#include "TList.h"
#include "TH1I.h"
#include "TChain.h"
//
#include "AliTracker.h"
#include "AliMagF.h"
//
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
//
//
#include "AliTPCcalibTracks.h"
#include "AliTPCcalibTracksGain.h"

#include "AliTPCSelectorESD.h" 
#include "AliTPCSelectorTracks.h" 
#include "TProof.h"

const char* AliTPCSelectorTracks::fgkOutputFileName = "Output.root";


AliTPCSelectorTracks::AliTPCSelectorTracks(TTree *) : 
  AliTPCSelectorESD(),
  fInit(kFALSE),
  fCalibTracks(0),
  fCalibTracksGain(0)
{
   // 
   // 
   // 
   G__SetCatchException(0);
}   

AliTPCSelectorTracks::~AliTPCSelectorTracks(){
  //
  //
  //
}

void AliTPCSelectorTracks::InitComponent(){
  //
  // Init Components
  //
  //
  // USER -COMPONENT definder part
  // 
  // before calling the process function, two objects have to be added to the chain/tree:
  // clusterParam, the cluster parametrization and cuts, a AliTPCcalibTracksCuts object
  // they have to be added in the following manner:
  //    chain->GetUserInfo()->AddLast(clusterParam);
  //    chain->GetUserInfo()->AddLast(cuts);
  // 
  static Int_t counter=0;
  if (!fChain){
    Error("InitComponent","ERROR - chain not initialized\n");
  }
  Info("InitComponent",Form("Selector initialization No\t%d\n", counter));
  counter++;
  //
  

  AliTPCClusterParam *clusterParam  = (AliTPCClusterParam*)fChain->GetUserInfo()->FindObject("AliTPCClusterParam");
  //
  if (clusterParam == 0) Error("InitComponent","CLUSTER PARAM NOT FOUND IN CHAIN! \n");
  //
  AliTPCcalibTracksCuts *cuts = (AliTPCcalibTracksCuts*)fChain->GetUserInfo()->FindObject("calibTracksCuts");
  if (cuts != 0) Info("InitComponent","cuts found in fChain! \n");
  else{
    Error("InitComponent","CUTS NOT FOUND IN CHAIN\n");
  }
  if (clusterParam==0 && fInput){
    clusterParam = (AliTPCClusterParam*)fInput->FindObject("AliTPCClusterParam");
    Error("InitComponent","CLUSTER PARAM NOT FOUND IN PROOF\n");
  }
  if (cuts==0 &&fInput ){    
    cuts = (AliTPCcalibTracksCuts*)fInput->FindObject("calibTracksCuts");
    Error("InitComponent","CUTS NOT FOUND IN PROOF\n");
  }
  if (!cuts || !clusterParam) {
    if (fInput) fInput->Print();
    return;
  }
   
  fCalibTracks = new AliTPCcalibTracks("calibTracks", "Resolution calibration object for tracks", clusterParam, cuts);
  fOutput->AddLast(fCalibTracks);
   
  AliTPCcalibTracksGain* prevIter = (AliTPCcalibTracksGain*)fChain->GetUserInfo()->FindObject("calibTracksGain");
  if (!prevIter) Info("InitComponent", "Previous iteration of calibTracksGain not found, continuing without.");
  TNamed* debugStreamPrefix = (TNamed*)fChain->GetUserInfo()->FindObject("debugStreamPrefix");
  fCalibTracksGain = new AliTPCcalibTracksGain("calibTracksGain", "Gain calibration object for tracks", cuts, debugStreamPrefix, prevIter);
  fOutput->AddLast(fCalibTracksGain);
  fInit=kTRUE;
}

void AliTPCSelectorTracks::SlaveBegin(TTree * tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
  
  AliTPCSelectorESD::SlaveBegin(tree);
 
  printf(" ***** SlaveBegin ***** \n");
}


void AliTPCSelectorTracks::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
  printf ("SlaveTerminate.. \n");
  printf ("Terminate CalibTracksGain.. \n");
  if (fCalibTracksGain) fCalibTracksGain->Terminate();
}




Int_t AliTPCSelectorTracks::ProcessIn(Long64_t entry)
{
   //
   //
   //
  if (!fInit) InitComponent();
  if (!fInit) return 0;
  Int_t status = ReadEvent(entry);
  if (status<0) return status; 
  Int_t ntracks = (fESD) ? fESD->GetNumberOfTracks() : fESDevent->GetNumberOfTracks();     
   //
   //
   // USER code to go here
   //
   AliTPCseed *seed;
   
   for (Int_t tr = 0; tr < ntracks; tr++){ 
      AliESDtrack *esdTrack = fESD ? (AliESDtrack*) fESD->GetTrack(tr): (AliESDtrack*) fESDevent->GetTrack(tr);
      AliESDfriendTrack *friendtrack = (AliESDfriendTrack*) esdTrack->GetFriendTrack();
      seed = 0;
      TObject *cobject = 0;
      for (Int_t i = 0; ; i++){
         cobject = friendtrack->GetCalibObject(i);
         if (!cobject) break;
         seed = dynamic_cast<AliTPCseed*>(cobject);
         if (seed) break;
      }
      
      if (seed) {
         fNClusters->Fill(seed->GetNumberOfClusters());
         //
         fCalibTracks->Process(seed, esdTrack);   // analysis is done in fCalibTracks
         fCalibTracksGain->Process(seed);
      }
   }
   CleanESD();
   return 0;
}


void AliTPCSelectorTracks::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   

   if (!fOutput) return;
   
   TFile file(fgkOutputFileName, "recreate");
   fCalibTracksGain  = (AliTPCcalibTracksGain*)fOutput->FindObject("calibTracksGain");
   // evaluate all fitters before saving them, because it doesn't seem to be possible
   // to evaluate a TLinearFitter after it has been loaded from a root file
   if (fCalibTracksGain) fCalibTracksGain->Evaluate();
   fOutput->Write();
   file.Close();
   printf("Successfully written file to '%s'.", fgkOutputFileName);
 

   Info("Destructor","Destuctor");
   //delete fCalibTracksGain;
   //delete fCalibTracks;
//   printf ("Terminate... \n");
//   if (!fOutput) return;
//   TFile file("Output.root","recreate");
//   printf("fOutput contains the following: \n");
//   fOutput->Print();
//   printf("Trying to write the file 'Output.root'... \n");
//   fOutput->Write();
//   file.Close();  
  
}


