/* $Id$ */

#include "AlidNdEtaCorrectionSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TParticlePDG.h>

#include <TChain.h>
#include <TSelector.h>

#include <AliLog.h>
#include <AliGenEventHeader.h>
#include <AliTracker.h>
#include <AliHeader.h>
#include <AliESDVertex.h>
#include <AliESD.h>
#include <AliESDtrack.h>

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "dNdEtaCorrection.h"

ClassImp(AlidNdEtaCorrectionSelector)

AlidNdEtaCorrectionSelector::AlidNdEtaCorrectionSelector() :
  AliSelectorRL(),
  fEsdTrackCuts(0),
  fdNdEtaCorrection(0),
  fdNdEtaCorrectionFinal(0)
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

  fdNdEtaCorrection = new dNdEtaCorrection();

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

  if (!fEsdTrackCuts)
  {
    AliDebug(AliLog::kError, "fESDTrackCuts not available");
    return kFALSE;
  }

  // ########################################################
  // get the EDS vertex
  const AliESDVertex* vtxESD = fESD->GetVertex();

  Bool_t goodEvent = kTRUE;

  // the vertex should be reconstructed
  if (strcmp(vtxESD->GetName(),"default")==0) 
    goodEvent = kFALSE;

  Double_t vtx_res[3];
  vtx_res[0] = vtxESD->GetXRes();
  vtx_res[1] = vtxESD->GetYRes();
  vtx_res[2] = vtxESD->GetZRes();
  
  // the resolution should be reasonable???
  if (vtx_res[2]==0 || vtx_res[2]>0.1)
    goodEvent = kFALSE;
    

  // ########################################################
  // get the MC vertex
  AliGenEventHeader* genHeader = header->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  fdNdEtaCorrection->FillEvent(vtxMC[2]);

  if (goodEvent)
    fdNdEtaCorrection->FillUsedEvent(vtxMC[2]);
    


  // ########################################################
  // loop over mc particles
  TTree* particleTree = GetKinematics();
  TParticle* particle = 0;
  particleTree->SetBranchAddress("Particles", &particle);

  Int_t nPrim  = header->GetNprimary();
  Int_t nTotal = header->GetNtrack();

  for (Int_t i_mc = nTotal - nPrim; i_mc < nTotal; ++i_mc)
  {
    particleTree->GetEntry(i_mc);

    if (!particle)
      continue;

    if (IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    Float_t eta = particle->Eta();
    
    fdNdEtaCorrection->FillParticleAllEvents(vtxMC[2], eta);	
    
    if (goodEvent)
      fdNdEtaCorrection->FillParticleWhenUsedEvent(vtxMC[2], eta);	
    
  }// end of mc particle

  if (!goodEvent)
    return kTRUE;  

  // ########################################################
  // loop over esd tracks
  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t t=0; t<nTracks; t++)
  {
    AliESDtrack* esdTrack = fESD->GetTrack(t);

    // cut the esd track?
    if (!fEsdTrackCuts->AcceptTrack(esdTrack))
      continue;

    // using the eta of the mc particle
    Int_t label = TMath::Abs(esdTrack->GetLabel());
    if (label == 0)
    {
      AliDebug(AliLog::kWarning, Form("WARNING: cannot find corresponding mc part for track %d.", t));
      continue;
    }
    particleTree->GetEntry(nTotal - nPrim + label);

    fdNdEtaCorrection->FillParticleWhenMeasuredTrack(vtxMC[2], particle->Eta());

  } // end of track loop

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

  fdNdEtaCorrectionFinal = dynamic_cast<dNdEtaCorrection*> (fOutput->FindObject("dndeta_correction"));

  fdNdEtaCorrectionFinal->Finish();

  TFile* fout = new TFile("correction_map.root","RECREATE");

  fEsdTrackCuts->SaveHistograms("esd_track_cuts");
  fdNdEtaCorrectionFinal->SaveHistograms();

  fout->Write();
  fout->Close();

  fdNdEtaCorrectionFinal->DrawHistograms();
}
