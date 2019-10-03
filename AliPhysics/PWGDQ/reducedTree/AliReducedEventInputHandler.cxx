//
//     Event handler for AliReducedEvent information
//     Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no
//

#include <TTree.h>
#include <TFile.h>
#include "AliReducedEventInputHandler.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedEventInfo.h"

ClassImp(AliReducedEventInputHandler)

//______________________________________________________________________________
AliReducedEventInputHandler::AliReducedEventInputHandler() :
    AliInputEventHandler(),
    fEventInputOption(kReducedBaseEvent),
    fReducedEvent(0)
{
  // Default constructor
}

//______________________________________________________________________________
AliReducedEventInputHandler::AliReducedEventInputHandler(const char* name, const char* title):
  AliInputEventHandler(name, title),
  fEventInputOption(kReducedBaseEvent),
  fReducedEvent(0)
 {
    // Constructor
}

//______________________________________________________________________________
AliReducedEventInputHandler::~AliReducedEventInputHandler() 
{
// Destructor
}

//______________________________________________________________________________
Bool_t AliReducedEventInputHandler::Init(TTree* tree, Option_t* opt)
{
    // Initialisation necessary for each new tree
    fTree = tree;
    if (!fTree) return kFALSE;
    fTree->GetEntries();

    SwitchOffBranches();
    SwitchOnBranches();
    
    // Get pointer to the event
    if (!fReducedEvent) {
       switch(fEventInputOption) {
          case kReducedEventInfo:
             fReducedEvent = new AliReducedEventInfo();
             break;
          default:   
             fReducedEvent = new AliReducedBaseEvent();   
       }
    }
    
    tree->SetBranchAddress("Event",&fReducedEvent);
    
    return kTRUE;
}


//______________________________________________________________________________
Bool_t AliReducedEventInputHandler::BeginEvent(Long64_t entry)
{
    // Begin event
    static Int_t prevRunNumber = -1;
    if (prevRunNumber != fReducedEvent->RunNo() ) {
      prevRunNumber = fReducedEvent->RunNo();
    } 
    fTree->GetEvent(entry);
    
    // set transient pointer to event inside tracks
    // fEvent->ConnectTracks();

    return kTRUE;
}


//______________________________________________________________________________
Bool_t AliReducedEventInputHandler::Notify(const char* path)
{
  // Notification of directory change
  //SwitchOffBranches();
  //SwitchOnBranches();
  //fUserInfo=fTree->GetTree()->GetUserInfo();
    
  //TTree *ttree = fTree->GetTree();
  //if (!ttree) ttree = fTree;
  //TString statFname(ttree->GetCurrentFile()->GetName());
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliReducedEventInputHandler::FinishEvent()
{
  // Finish event
  if (fReducedEvent) fReducedEvent->ClearEvent();
  return kTRUE;
}
