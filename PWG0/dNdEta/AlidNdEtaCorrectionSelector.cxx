/* $Id$ */

#include "AlidNdEtaCorrectionSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TParticlePDG.h>

#include <TChain.h>
#include <TSelector.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliGenEventHeader.h>
#include <AliTracker.h>
#include <AliHeader.h>
#include <AliESDVertex.h>
#include <AliESD.h>
#include <AliESDtrack.h>
#include <AliRunLoader.h>
#include <AliStack.h>

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "dNdEta/AlidNdEtaCorrection.h"
#include "AliPWG0Helper.h"

ClassImp(AlidNdEtaCorrectionSelector)

AlidNdEtaCorrectionSelector::AlidNdEtaCorrectionSelector() :
  AliSelectorRL(),
  fEsdTrackCuts(0),
  fdNdEtaCorrection(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AlidNdEtaCorrectionSelector::~AlidNdEtaCorrectionSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AlidNdEtaCorrectionSelector::Begin(TTree * tree)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelectorRL::Begin(tree);
}

void AlidNdEtaCorrectionSelector::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelectorRL::SlaveBegin(tree);

  fdNdEtaCorrection = new AlidNdEtaCorrection();

  if (fTree)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fTree->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
    AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from user info");
}

Bool_t AlidNdEtaCorrectionSelector::Process(Long64_t entry)
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

  // check prerequesites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  AliHeader* header = GetHeader();
  if (!header)
  {
    AliDebug(AliLog::kError, "Header not available");
    return kFALSE;
  }

  AliRunLoader* runLoader = GetRunLoader();
  if (!runLoader)
  {
    AliDebug(AliLog::kError, "RunLoader not available");
    return kFALSE;
  }

  runLoader->LoadKinematics();
  AliStack* stack = runLoader->Stack();
  if (!stack)
  {
    AliDebug(AliLog::kError, "Stack not available");
    return kFALSE;
  }

  if (!fEsdTrackCuts)
  {
    AliDebug(AliLog::kError, "fESDTrackCuts not available");
    return kFALSE;
  }

  Bool_t vertexReconstructed = AliPWG0Helper::IsVertexReconstructed(fESD);

  Bool_t eventTriggered = AliPWG0Helper::IsEventTriggered(fESD);

  fdNdEtaCorrection->IncreaseEventCount();
  if (eventTriggered)
    fdNdEtaCorrection->IncreaseTriggeredEventCount();

  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  // loop over mc particles
  Int_t nPrim  = stack->GetNprimary();

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
    Float_t pt = particle->Pt();

    if (vertexReconstructed)
      fdNdEtaCorrection->FillParticle(vtxMC[2], eta, pt);

    fdNdEtaCorrection->FillParticleAllEvents(eta, pt);
    if (eventTriggered)
      fdNdEtaCorrection->FillParticleWhenEventTriggered(eta, pt);
  }// end of mc particle

  // ########################################################
  // loop over esd tracks
  Int_t nTracks = fESD->GetNumberOfTracks();

  // count the number of "good" tracks for vertex reconstruction efficiency
  // TODO change to number of ITS clusters or similar
  Int_t nGoodTracks = 0;
  for (Int_t t=0; t<nTracks; t++)
  {
    AliDebug(AliLog::kDebug+1, Form("ESD Loop: Processing track %d.", t));

    AliESDtrack* esdTrack = fESD->GetTrack(t);

    // cut the esd track?
    if (!fEsdTrackCuts->AcceptTrack(esdTrack))
      continue;

    nGoodTracks++;

    // using the properties of the mc particle
    Int_t label = TMath::Abs(esdTrack->GetLabel());
    if (label == 0)
    {
      AliDebug(AliLog::kWarning, Form("WARNING: cannot find corresponding mc part for track %d.", t));
      continue;
    }

    TParticle* particle = stack->Particle(label);
    if (!particle)
    {
      AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (track loop).", label));
      continue;
    }

    if (vertexReconstructed)
      fdNdEtaCorrection->FillParticleWhenMeasuredTrack(vtxMC[2], particle->Eta(), particle->Pt());
  } // end of track loop


  fdNdEtaCorrection->FillEvent(vtxMC[2], nGoodTracks);
  if (eventTriggered)
  {
    fdNdEtaCorrection->FillEventWithTrigger(vtxMC[2], nGoodTracks);
    if (vertexReconstructed)
      fdNdEtaCorrection->FillEventWithTriggerWithReconstructedVertex(vtxMC[2], nGoodTracks);
  }

  return kTRUE;
}

void AlidNdEtaCorrectionSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelectorRL::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, "ERROR: Output list not initialized");
    return;
  }

  fOutput->Add(fdNdEtaCorrection);
}

void AlidNdEtaCorrectionSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelectorRL::Terminate();

  fdNdEtaCorrection = dynamic_cast<AlidNdEtaCorrection*> (fOutput->FindObject("dndeta_correction"));

  fdNdEtaCorrection->Finish();

  TFile* fout = new TFile("correction_map.root","RECREATE");

  fEsdTrackCuts->SaveHistograms("esd_track_cuts");
  fdNdEtaCorrection->SaveHistograms();

  fout->Write();
  fout->Close();

  fdNdEtaCorrection->DrawHistograms();
}
