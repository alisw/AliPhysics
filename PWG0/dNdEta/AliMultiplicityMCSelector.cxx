/* $Id$ */

#include "AliMultiplicityMCSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TParticle.h>

#include <AliLog.h>
#include <AliESD.h>
#include <AliStack.h>

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "AliPWG0Helper.h"

ClassImp(AliMultiplicityMCSelector)

AliMultiplicityMCSelector::AliMultiplicityMCSelector() :
  AliSelectorRL(),
  fMultiplicityESD(0),
  fMultiplicityMC(0),
  fCorrelation(0),
  fEsdTrackCuts(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AliMultiplicityMCSelector::~AliMultiplicityMCSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AliMultiplicityMCSelector::Begin(TTree* tree)
{
  // Begin function

  AliSelectorRL::Begin(tree);

  ReadUserObjects(tree);
}

void AliMultiplicityMCSelector::ReadUserObjects(TTree* tree)
{
  // read the user objects, called from slavebegin and begin

  if (!fEsdTrackCuts && fInput)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts && tree)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from input list.");
}

void AliMultiplicityMCSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  ReadUserObjects(tree);

  fMultiplicityESD = new TH1F("fMultiplicityESD", "fMultiplicityESD;Ntracks;Count", 201, -0.5, 200.5);
  fMultiplicityMC = new TH1F("fMultiplicityMC", "fMultiplicityMC;Npart;Count", 201, -0.5, 200.5);

  fCorrelation = new TH2F("fCorrelation", "fCorrelation;Npart;Ntracks", 201, -0.5, 200.5, 201, -0.5, 200.5);
}

Bool_t AliMultiplicityMCSelector::Process(Long64_t entry)
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
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  if (!fEsdTrackCuts)
  {
    AliDebug(AliLog::kError, "fESDTrackCuts not available");
    return kFALSE;
  }

  AliStack* stack = GetStack();
  if (!stack)
  {
    AliDebug(AliLog::kError, "Stack not available");
    return kFALSE;
  }

  if (AliPWG0Helper::IsEventTriggered(fESD) == kFALSE)
    return kTRUE;

  if (AliPWG0Helper::IsVertexReconstructed(fESD) == kFALSE)
    return kTRUE;

  // TODO put pt cut

  // get number of "good" tracks from ESD
  Int_t nESDTracks = fEsdTrackCuts->CountAcceptedTracks(fESD);
  fMultiplicityESD->Fill(nESDTracks);

  // get number of tracks from MC
  // loop over mc particles
  Int_t nPrim  = stack->GetNprimary();
  Int_t nMCTracks = 0;
  for (Int_t iMc = 0; iMc < nPrim; ++iMc)
  {
    AliDebug(AliLog::kDebug+1, Form("MC Loop: Processing particle %d.", iMc));

    TParticle* particle = stack->Particle(iMc);

    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (mc loop).", iMc));
      continue;
    }

    if (AliPWG0Helper::IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    Float_t eta = particle->Eta();

    if (TMath::Abs(eta) > 0.9)
      continue;

    nMCTracks++;
  }// end of mc particle

  fMultiplicityMC->Fill(nMCTracks);

  fCorrelation->Fill(nMCTracks, nESDTracks);

  return kTRUE;
}

void AliMultiplicityMCSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelector::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, Form("ERROR: Output list not initialized."));
    return;
  }

  fOutput->Add(fMultiplicityESD);
  fOutput->Add(fMultiplicityMC);
  fOutput->Add(fCorrelation);
}

void AliMultiplicityMCSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  fOutput->Print();

  fMultiplicityESD = dynamic_cast<TH1F*> (fOutput->FindObject("fMultiplicityESD"));
  fMultiplicityMC = dynamic_cast<TH1F*> (fOutput->FindObject("fMultiplicityMC"));
  fCorrelation = dynamic_cast<TH2F*> (fOutput->FindObject("fCorrelation"));

  if (!fMultiplicityESD || !fMultiplicityMC || !fCorrelation)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p %p %p", (void*) fMultiplicityESD, (void*) fMultiplicityMC, (void*) fCorrelation));
    return;
  }

  TFile* file = TFile::Open("multiplicityMC.root", "RECREATE");

  fMultiplicityMC->Write();
  fMultiplicityESD->Write();
  fCorrelation->Write();

  file->Close();
}
