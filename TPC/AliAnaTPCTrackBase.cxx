//
// This class is ment as a base class for doing analysis of
// reconstructed TPC tracks in the AliAnalysisTask framework.
//
// The method FillTrackHistograms should be overloaded by the users
// class and used to make cuts on tracks
//
// Questions, comments, or suggestions can be send to Peter Christiansen (Lund)
// 

// ROOT includes
#include <TChain.h>
#include <TMath.h>

// ALIROOT includes
#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include <AliESD.h>

// STL includes
#include <iostream>
//
#include "AliAnaTPCTrackBase.h"



using namespace std;

ClassImp(AliAnaTPCTrackBase)

//________________________________________________________________________
AliAnaTPCTrackBase::AliAnaTPCTrackBase() : 
  AliAnalysisTask(),  
  fDebug(0),          //  Debug flag
  fESD(0), 
  fESDfriend(0), 
  fListOfHists(0),
  fMaxTracks(0),      // Max tracks in histogram
  fESDTracks(0),      //! N ESD tracks
  fGoodTracks(0)     //! GOOD tracks
{
  //
  // Default constructor (should not be used)
  //
  SetMaxTracks();
}

AliAnaTPCTrackBase::AliAnaTPCTrackBase(const AliAnaTPCTrackBase & ana):
  AliAnalysisTask(ana),  
  fDebug(ana.fDebug),          //  Debug flag
  fESD(ana.fESD), 
  fESDfriend(ana.fESDfriend), 
  fListOfHists(ana.fListOfHists),
  fMaxTracks(ana.fMaxTracks),      // Max tracks in histogram
  fESDTracks(0),      //! N ESD tracks
  fGoodTracks(0)     //! GOOD tracks
{
  //
  // copy constructor
  //  
  fESDTracks  = (TH1F*)ana.fESDTracks->Clone();      //! N ESD tracks
  fGoodTracks = (TH1F*)ana.fGoodTracks->Clone();     //! GOOD tracks
}


AliAnaTPCTrackBase& AliAnaTPCTrackBase::operator=(const AliAnaTPCTrackBase&ana){
  //
  // assignemnt operator
  //
  if (this != &ana) {
    new (this) AliAnaTPCTrackBase(ana);
  }
  return *this;


}



//________________________________________________________________________
AliAnaTPCTrackBase::AliAnaTPCTrackBase(const char *name) : 
  AliAnalysisTask(name, "AliAnaTPCTrackBase"), 
  fDebug(0),          //  Debug flag
  fESD(0), 
  fESDfriend(0), 
  fListOfHists(0),
  fMaxTracks(0),      // Max tracks in histogram
  fESDTracks(0),      //! N ESD tracks
  fGoodTracks(0)     //! GOOD tracks
{
  //
  // Normal constructor
  //

  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList
  DefineOutput(0, TList::Class());
  
  fDebug = 0;
  SetMaxTracks();
}

//________________________________________________________________________
void AliAnaTPCTrackBase::ConnectInputData(Option_t *) 
{
  //
  // Connect the input data
  //
  if(fDebug>3)
    cout << "AnalysisTaskTPCCluster::ConnectInputData()" << endl;

  AliESDInputHandler* esdH = (AliESDInputHandler*) 
    ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  fESD = esdH->GetEvent();

  if(fESD==0) {
    
    cout << endl << "WARNING: NO ESD event found" << endl << endl;
  } else {

    fESDfriend = 
      (AliESDfriend*)fESD->FindListObject("AliESDfriend");  
    
    if(fESDfriend==0) {
      
      cout << endl << "WARNING: NO ESD friend found" << endl << endl;
    }
  }
}

//________________________________________________________________________
void AliAnaTPCTrackBase::CreateOutputObjects() 
{
  //
  // Connect the output objects
  //
  if(fDebug>3)
    cout << "AnalysisTaskTPCCluster::CreateOutputObjects()" << endl;
  
  OpenFile(0);
  fListOfHists = new TList();
  
  fESDTracks = 
    new TH1F("hESDTracks", 
	     "Number of ESD tracks per event; N ESD tracks; Counts", 
	     TMath::Min(fMaxTracks, 100), 0, fMaxTracks);
  fListOfHists->Add(fESDTracks);

  fGoodTracks = 
    new TH1F("hGoodTracks", 
	     "Number of Good tracks per event; N good tracks; Counts", 
	     TMath::Min(fMaxTracks, 100), 0, fMaxTracks);
  fListOfHists->Add(fGoodTracks);
}


//________________________________________________________________________
void AliAnaTPCTrackBase::Exec(Option_t *) {
  //
  // Execute analysis for current event 
  // For the moment I require that mcTruth is there!  I am afraid to
  // get out of sync if it is missing for some events since I use the
  // number of MC events for normalisation
  //

  if(fDebug>3)
    cout << "AliAnaTPCTrackBase::Exec()" << endl;
  
  if(fESD==0 || fESDfriend==0) {
    
    cout << "AliAnaTPCTrackBase::Exec(): WARNING: fESD=" << fESD 
	 << ", fESDfriend=" << fESDfriend << endl;
    // Post final data. It will be written to a file with option "RECREATE"
    PostData(0, fListOfHists);
    return;
  }
 
  fESD->SetESDfriend(fESDfriend);
  const Int_t nESDTracks = fESD->GetNumberOfTracks();
  
  if(fDebug>0)
    cout << "          Number of ESD tracks: " << nESDTracks << endl;
  
  if ( nESDTracks != fESDfriend->GetNumberOfTracks() ) {
     AliWarning("Number of Tracks differs from Number of Friend-Tracks!");
     printf("Number of tracks: %i, number of friend tracks: %i \n", nESDTracks, fESDfriend->GetNumberOfTracks());
     return;
  }
  
  fESDTracks->Fill(nESDTracks);
  Int_t nGoodTracks = 0;

  for(Int_t i = 0; i < nESDTracks; i++) {
    
    AliESDtrack* track = (AliESDtrack*)fESD->GetTrack(i);
    AliESDfriendTrack* friendTrack = (AliESDfriendTrack*) track->GetFriendTrack();
    AliTPCseed* seed = 0;
    
    if(friendTrack) {
      TObject *cobject = 0;
      for (Int_t i = 0; ; i++){
         cobject = friendTrack->GetCalibObject(i);
         if (!cobject) break;
         seed = dynamic_cast<AliTPCseed*>(cobject);
         if (seed) break;
      }
      if (!seed && fDebug>1) Error("Exec", "No seed found!!!");
    }
    else if (fDebug>1) Error("Exec", "No friend track found!!!");

    Int_t accepted = FillTrackHistograms(nESDTracks, track, friendTrack, seed);

    if(accepted) {
      if(fDebug>1)
	cout << "Track " << i << " was accepted" << endl;
      nGoodTracks++;
    } else {

      if(fDebug>1)
	cout << "Track " << i << " was rejected" << endl;
    }
  }
  
  fGoodTracks->Fill(nGoodTracks);
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(0, fListOfHists);
}      

//________________________________________________________________________
Int_t AliAnaTPCTrackBase::FillTrackHistograms(Int_t nTracks, AliESDtrack* track, AliESDfriendTrack* friendTrack, AliTPCseed* seed) {
  //
  // This method should be overloaded and used to make cuts on tracks
  // and fill histograms. 
  // return 0 if track was rejected, 1 if accepted
  //

  if(nTracks && track && friendTrack && seed)
    return 1;
  else
    return 0;
}

//________________________________________________________________________
void AliAnaTPCTrackBase::Terminate(Option_t *) {
    //
    // Terminate loop
    //
  if(fDebug>3)
    printf("AliAnaTPCTrackBase: Terminate() \n");
}
