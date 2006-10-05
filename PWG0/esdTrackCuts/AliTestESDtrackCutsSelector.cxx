/* $Id$ */

#include "AliTestESDtrackCutsSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

#include <TSelector.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliRunLoader.h>
#include <AliStack.h>

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "AliPWG0Helper.h"


ClassImp(AliTestESDtrackCutsSelector)

AliTestESDtrackCutsSelector::AliTestESDtrackCutsSelector() :
  AliSelectorRL(),
  fEsdTrackCutsAll(0),
  fEsdTrackCutsPri(0),
  fEsdTrackCutsSec(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AliTestESDtrackCutsSelector::~AliTestESDtrackCutsSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AliTestESDtrackCutsSelector::Begin(TTree* tree)
{
  // Begin function

  AliSelectorRL::Begin(tree);

  ReadUserObjects(tree);
}

void AliTestESDtrackCutsSelector::ReadUserObjects(TTree* tree)
{
  // read the user objects, called from slavebegin and begin

  if (!fEsdTrackCutsAll && fInput)
    fEsdTrackCutsAll = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("esdTrackCutsAll")->Clone());
  if (!fEsdTrackCutsPri && fInput)
    fEsdTrackCutsPri = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("esdTrackCutsPri")->Clone());
  if (!fEsdTrackCutsSec && fInput)
    fEsdTrackCutsSec = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("esdTrackCutsSec")->Clone());

  if (!fEsdTrackCutsAll && tree)
    fEsdTrackCutsAll = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("esdTrackCutsAll"));
  if (!fEsdTrackCutsPri && tree)
    fEsdTrackCutsPri = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("esdTrackCutsPri"));
  if (!fEsdTrackCutsSec && tree)
    fEsdTrackCutsSec = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("esdTrackCutsSec"));

  if (!fEsdTrackCutsAll || !fEsdTrackCutsPri || !fEsdTrackCutsSec)
     AliDebug(AliLog::kError, "ERROR: Could not read esdTrackCutsXXX from input list.");
}

void AliTestESDtrackCutsSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelectorRL::SlaveBegin(tree);

  ReadUserObjects(tree);
}

Bool_t AliTestESDtrackCutsSelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fTree is the pointer to the TChain being processed,
  //  use fTree->GetTree()->GetEntry(entry).

  if (AliSelectorRL::Process(entry) == kFALSE)
    return kFALSE;

  // Check prerequisites
  if (!fESD) {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }  

  if (!AliPWG0Helper::IsVertexReconstructed(fESD)) {
    AliDebug(AliLog::kDebug+5, "Vertex is not reconstructed");
    return kFALSE;
  }

  // check if the esd track cut objects are there
  if (!fEsdTrackCutsAll || !fEsdTrackCutsPri || !fEsdTrackCutsSec) {
    AliDebug(AliLog::kError, "fEsdTrackCutsXXX not available");
    return kFALSE;
  }

  // get particle stack
  AliStack* stack = GetStack();
  if (!stack) {
    AliDebug(AliLog::kError, "Stack not available");
    return kFALSE;
  }
  Int_t nPrim  = stack->GetNprimary();
  
  // ########################################################
  // loop over esd tracks
  Int_t nTracks = fESD->GetNumberOfTracks();

  // count the number of "good" tracks as parameter for vertex reconstruction efficiency
  for (Int_t t=0; t<nTracks; t++) {
    AliDebug(AliLog::kDebug+1, Form("ESD Loop: Processing track %d.", t));

    AliESDtrack* esdTrack = fESD->GetTrack(t);

    fEsdTrackCutsAll->AcceptTrack(esdTrack);

    // using the properties of the mc particle
    Int_t label = TMath::Abs(esdTrack->GetLabel());
    if (label == 0) {
      AliDebug(AliLog::kWarning, Form("WARNING: cannot find corresponding mc part for track %d.", t));
      continue;
    }
    TParticle* particle = stack->Particle(label);
    if (!particle) {
      AliDebug(AliLog::kError, Form("UNEXPECTED: part with label %d not found in stack (track loop).", label));
      continue;
    }

    if (label < nPrim)
      fEsdTrackCutsPri->AcceptTrack(esdTrack);
    else
      fEsdTrackCutsSec->AcceptTrack(esdTrack);
  }
  
  return kTRUE;
}

void AliTestESDtrackCutsSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelectorRL::SlaveTerminate();
  
  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, Form("ERROR: Output list not initialized."));
    return;
  }

  fOutput->Add(fEsdTrackCutsAll);
  fOutput->Add(fEsdTrackCutsPri);
  fOutput->Add(fEsdTrackCutsSec);
}

void AliTestESDtrackCutsSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelectorRL::Terminate();

  if (fOutput)
    fOutput->Print();

  fEsdTrackCutsAll = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("esdTrackCutsAll"));
  fEsdTrackCutsPri = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("esdTrackCutsPri"));
  fEsdTrackCutsSec = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("esdTrackCutsSec"));

  // check if the esd track cut objects are there
  if (!fEsdTrackCutsAll || !fEsdTrackCutsPri || !fEsdTrackCutsSec) {
    AliDebug(AliLog::kError, Form("fEsdTrackCutsXXX not available %p %p %p", fEsdTrackCutsAll, fEsdTrackCutsPri, fEsdTrackCutsSec));
    return;
  }

  TFile* file = TFile::Open("trackCuts.root", "RECREATE");
  fEsdTrackCutsAll->SaveHistograms("esdTrackCutsAll");
  fEsdTrackCutsPri->SaveHistograms("esdTrackCutsPri");
  fEsdTrackCutsSec->SaveHistograms("esdTrackCutsSec");

  file->Close();
}
