//
// This class is the task for connecting together 
// MC information and the RC information 
//
// The task is a wrapper over two components
// AliGenInfoMaker
// AliESDRecInfoMaker.h

// ROOT includes
#include <TChain.h>
#include <TMath.h>

// ALIROOT includes
#include <AliAnalysisManager.h>
#include <AliESDInputHandler.h>
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

#include <AliESD.h>
#include "AliGenInfoTask.h"
#include "AliGenInfoMaker.h"
#include "AliHelix.h"

// STL includes
#include <iostream>

using namespace std;

ClassImp(AliGenInfoTask)

//________________________________________________________________________
AliGenInfoTask::AliGenInfoTask() : 
  AliAnalysisTask(), fGenMaker(0),
  fESD(0), fESDfriend(0), fListOfHists(0)
{
  //
  // Default constructor (should not be used)
  //
  fDebug = 0;
  SetMaxTracks();
}

//________________________________________________________________________
AliGenInfoTask::AliGenInfoTask(const char *name) : 
  AliAnalysisTask(name, "AliGenInfoTask"), 
  fGenMaker(0),
  fESD(0), 
  fESDfriend(0), 
  fListOfHists(0)
{
  //
  // Normal constructor
  //
  fGenMaker = new AliGenInfoMaker;
  fGenMaker->SetIO();
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList
  DefineOutput(0, TList::Class());
  
  fDebug = 0;
  SetMaxTracks();
}

//________________________________________________________________________
void AliGenInfoTask::ConnectInputData(Option_t *) 
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
void AliGenInfoTask::CreateOutputObjects() 
{
  //
  // Connect the output objects
  //
  if(fDebug>3)
    cout << "AnalysisTaskTPCCluster::CreateOutputObjects()" << endl;
  
//    OpenFile(0);
//   fListOfHists = new TList();
  
//   hESDTracks = 
//     new TH1F("hESDTracks", 
// 	     "Number of ESD tracks per event; N ESD tracks; Counts", 
// 	     TMath::Min(fMaxTracks, 100), 0, fMaxTracks);
//   fListOfHists->Add(hESDTracks);

//   hGoodTracks = 
//     new TH1F("hGoodTracks", 
// 	     "Number of Good tracks per event; N good tracks; Counts", 
// 	     TMath::Min(fMaxTracks, 100), 0, fMaxTracks);
//   fListOfHists->Add(hGoodTracks);
}


//________________________________________________________________________
void AliGenInfoTask::Exec(Option_t *) {
  //
  // Execute analysis for current event 
  // For the moment I require that mcTruth is there!  I am afraid to
  // get out of sync if it is missing for some events since I use the
  // number of MC events for normalisation
  //

  if(fDebug>3)
    cout << "AliGenInfoTask::Exec()" << endl;

  if (fESD) {
    AliHelix::SetBz(fESD->GetMagneticField());
  }
  
  // Monte carlo info
  AliMCEventHandler* mcinfo = (AliMCEventHandler*) (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());  
  // If MC has been connected   
  
  if (!mcinfo){
    cout << "Not MC info\n" << endl;
  }
  fGenMaker->ProcessEvent(mcinfo);



  if(fESD==0 || fESDfriend==0) {
    
    cout << "AliGenInfoTask::Exec(): WARNING: fESD=" << fESD 
	 << ", fESDfriend=" << fESDfriend << endl;
    // Post final data. It will be written to a file with option "RECREATE"
    PostData(0, fListOfHists);
    return;
  }
 
  fESD->SetESDfriend(fESDfriend);
  const Int_t nESDTracks = fESD->GetNumberOfTracks();
  
  if(fDebug>0){
    cout << " AliGenIfoTask::Exec() Number of ESD tracks: " << nESDTracks << endl;
  }
  if ( nESDTracks != fESDfriend->GetNumberOfTracks() ) {
     AliWarning("Number of Tracks differs from Number of Friend-Tracks!");
     printf("Number of tracks: %i, number of friend tracks: %i \n", nESDTracks, fESDfriend->GetNumberOfTracks());
     return;
  }
//   // Post final data. It will be written to a file with option "RECREATE"
  PostData(0, fListOfHists);
}      

//________________________________________________________________________
Int_t AliGenInfoTask::FillTrackHistograms(Int_t nTracks, AliESDtrack* track, AliESDfriendTrack* friendTrack, AliTPCseed* seed) {
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
void AliGenInfoTask::Terminate(Option_t *) {
    //
    // Terminate loop
    //
  if(fDebug>3)
    printf("AliGenInfoTask: Terminate() \n");  
  fGenMaker->CloseOutputFile();
}
