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
#include <TH3F.h>

#include <TSelector.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliRunLoader.h>
#include <AliStack.h>

#include "AliESDtrackCuts.h"
#include "AliPWG0Helper.h"


ClassImp(AliTestESDtrackCutsSelector)

AliTestESDtrackCutsSelector::AliTestESDtrackCutsSelector() :
  AliSelectorRL(),
  fEsdTrackCutsAll(0),
  fEsdTrackCutsNoVtx(0),
  fEsdTrackCutsPri(0),
  fEsdTrackCutsSec(0),
  fEsdTrackCutsPlusZ(0),
  fEsdTrackCutsMinusZ(0),
  fEsdTrackCutsPos(0),
  fEsdTrackCutsNeg(0),
  fPIDAfterCutNoVtx(0),
  fPIDAfterCutAll(0),
  fVertex(0)
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

  // only do it once
  if (fEsdTrackCutsAll)
    return;

  if (!fEsdTrackCutsAll && fInput)
    fEsdTrackCutsAll = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("esdTrackCutsAll")->Clone());

  if (!fEsdTrackCutsAll && tree)
    fEsdTrackCutsAll = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("esdTrackCutsAll"));

  if (!fEsdTrackCutsAll)
     AliDebug(AliLog::kError, "ERROR: Could not read fEsdTrackCutsAll from input list.");

  fEsdTrackCutsNoVtx =    dynamic_cast<AliESDtrackCuts*> (fEsdTrackCutsAll->Clone("fEsdTrackCutsNoVtx"));
  fEsdTrackCutsNoVtx->SetRequireSigmaToVertex(kFALSE);

  fEsdTrackCutsPri =    dynamic_cast<AliESDtrackCuts*> (fEsdTrackCutsAll->Clone("fEsdTrackCutsPri"));
  fEsdTrackCutsSec =    dynamic_cast<AliESDtrackCuts*> (fEsdTrackCutsAll->Clone("fEsdTrackCutsSec"));
  fEsdTrackCutsPlusZ =  dynamic_cast<AliESDtrackCuts*> (fEsdTrackCutsAll->Clone("fEsdTrackCutsPlusZ"));
  fEsdTrackCutsMinusZ = dynamic_cast<AliESDtrackCuts*> (fEsdTrackCutsAll->Clone("fEsdTrackCutsMinusZ"));
  fEsdTrackCutsPos =    dynamic_cast<AliESDtrackCuts*> (fEsdTrackCutsAll->Clone("fEsdTrackCutsPos"));
  fEsdTrackCutsNeg =    dynamic_cast<AliESDtrackCuts*> (fEsdTrackCutsAll->Clone("fEsdTrackCutsNeg"));
}

void AliTestESDtrackCutsSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelectorRL::SlaveBegin(tree);

  ReadUserObjects(tree);

  fPIDAfterCutNoVtx = new TH1F("fPIDAfterCutNoVtx", "fPIDAfterCutNoVtx", 5001, -2500.5, 2500.5);
  fPIDAfterCutAll   = new TH1F("fPIDAfterCutAll", "fPIDAfterCutAll", 5001, -2500.5, 2500.5);

  fVertex = new TH3F("fVertex", "fVertex", 100, -10, 10, 100, -10, 10, 100, -10, 10);
}

void AliTestESDtrackCutsSelector::Init(TTree* tree)
{
  // read the user objects

  AliSelectorRL::Init(tree);

  // Enable only the needed branches
  if (tree)
  {
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("fTriggerMask", 1);
    tree->SetBranchStatus("fSPDVertex*", 1);
    tree->SetBranchStatus("fTracks.fLabel", 1);

    AliESDtrackCuts::EnableNeededBranches(tree);
  }
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
  if (!fEsdTrackCutsAll || !fEsdTrackCutsPri || !fEsdTrackCutsSec || !fEsdTrackCutsPlusZ || !fEsdTrackCutsMinusZ || !fEsdTrackCutsPos || !fEsdTrackCutsNeg) {
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

    Bool_t passed = fEsdTrackCutsAll->AcceptTrack(esdTrack);

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
    {
      fEsdTrackCutsSec->AcceptTrack(esdTrack);
      if (passed)
      {
        fPIDAfterCutAll->Fill(particle->GetPdgCode());
        fVertex->Fill(particle->Vx(), particle->Vy(), particle->Vz());
      }

      if (fEsdTrackCutsNoVtx->AcceptTrack(esdTrack))
        fPIDAfterCutNoVtx->Fill(particle->GetPdgCode());
    }

    TParticlePDG* pdgPart = particle->GetPDG();
    if (pdgPart)
    {
      if (pdgPart->Charge() > 0)
        fEsdTrackCutsPos->AcceptTrack(esdTrack);
      else if (pdgPart->Charge() < 0)
        fEsdTrackCutsNeg->AcceptTrack(esdTrack);
    }

    if (particle->Eta() < 0)
      fEsdTrackCutsPlusZ->AcceptTrack(esdTrack);
    else
      fEsdTrackCutsMinusZ->AcceptTrack(esdTrack);
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
  fOutput->Add(fEsdTrackCutsNoVtx);
  fOutput->Add(fEsdTrackCutsPri);
  fOutput->Add(fEsdTrackCutsSec);
  fOutput->Add(fEsdTrackCutsPlusZ);
  fOutput->Add(fEsdTrackCutsMinusZ);
  fOutput->Add(fEsdTrackCutsPos);
  fOutput->Add(fEsdTrackCutsNeg);
  fOutput->Add(fPIDAfterCutNoVtx);
  fOutput->Add(fPIDAfterCutAll);
  fOutput->Add(fVertex);
}

void AliTestESDtrackCutsSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelectorRL::Terminate();

  fEsdTrackCutsAll = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("esdTrackCutsAll"));
  fEsdTrackCutsNoVtx = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCutsNoVtx"));
  fEsdTrackCutsPri = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCutsPri"));
  fEsdTrackCutsSec = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCutsSec"));
  fEsdTrackCutsPlusZ = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCutsPlusZ"));
  fEsdTrackCutsMinusZ = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCutsMinusZ"));
  fEsdTrackCutsPos = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCutsPos"));
  fEsdTrackCutsNeg = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fEsdTrackCutsNeg"));
  fPIDAfterCutNoVtx = dynamic_cast<TH1F*> (fOutput->FindObject("fPIDAfterCutNoVtx"));
  fPIDAfterCutAll = dynamic_cast<TH1F*> (fOutput->FindObject("fPIDAfterCutAll"));
  fVertex = dynamic_cast<TH3F*> (fOutput->FindObject("fVertex"));

  // check if the esd track cut objects are there
  if (!fEsdTrackCutsAll || !fEsdTrackCutsPri || !fEsdTrackCutsSec || !fEsdTrackCutsPlusZ || !fEsdTrackCutsMinusZ || !fEsdTrackCutsPos || !fEsdTrackCutsNeg) {
    AliDebug(AliLog::kError, Form("fEsdTrackCutsXXX not available %p %p %p %p %p %p %p", fEsdTrackCutsAll, fEsdTrackCutsPri, fEsdTrackCutsSec, fEsdTrackCutsPlusZ, fEsdTrackCutsMinusZ, fEsdTrackCutsPos, fEsdTrackCutsNeg));
    return;
  }

  TFile* file = TFile::Open("trackCuts.root", "RECREATE");

  fEsdTrackCutsAll->SaveHistograms();
  fEsdTrackCutsNoVtx->SaveHistograms();
  fEsdTrackCutsPri->SaveHistograms();
  fEsdTrackCutsSec->SaveHistograms();
  fEsdTrackCutsPlusZ->SaveHistograms();
  fEsdTrackCutsMinusZ->SaveHistograms();
  fEsdTrackCutsPos->SaveHistograms();
  fEsdTrackCutsNeg->SaveHistograms();
  fPIDAfterCutNoVtx->Write();
  fPIDAfterCutAll->Write();
  fVertex->Write();

  file->Close();

	fEsdTrackCutsAll->DrawHistograms();
}
