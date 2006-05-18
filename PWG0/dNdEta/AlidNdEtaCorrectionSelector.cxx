#include "AlidNdEtaCorrectionSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TParticlePDG.h>

#include <AliLog.h>
#include <AliGenEventHeader.h>
#include <AliTracker.h>

#include "../esdTrackCuts/AliESDtrackCuts.h"
#include "dNdEtaCorrection.h"

ClassImp(AlidNdEtaCorrectionSelector)

AlidNdEtaCorrectionSelector::AlidNdEtaCorrectionSelector(TTree *) :
  AliSelector(),
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

  AliSelector::Begin(tree);
}

void AlidNdEtaCorrectionSelector::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  fdNdEtaCorrection = new dNdEtaCorrection();

  if (fChain)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fChain->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
    AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from user info");
}

Bool_t AlidNdEtaCorrectionSelector::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. Typically here the branch pointers
  // will be retrieved. It is normaly not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed.

  if (AliSelector::Notify() == kFALSE)
    return kFALSE;

  return kTRUE;
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
  //  Assuming that fChain is the pointer to the TChain being processed,
  //  use fChain->GetTree()->GetEntry(entry).

  if (AliSelector::Process(entry) == kFALSE)
    return kFALSE;

  // check prerequesites
  if (!fESD)
  {
    AliDebug(AliLog::kError, "ESD branch not available");
    return kFALSE;
  }

  if (!fHeader)
  {
    AliDebug(AliLog::kError, "Header branch not available");
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

  // the vertex should be reconstructed
  if (strcmp(vtxESD->GetName(),"default")==0)
    return kTRUE;

  Double_t vtx_res[3];
  vtx_res[0] = vtxESD->GetXRes();
  vtx_res[1] = vtxESD->GetYRes();
  vtx_res[2] = vtxESD->GetZRes();

  // the resolution should be reasonable???
  if (vtx_res[2]==0 || vtx_res[2]>0.1)
    return kTRUE;

  // ########################################################
  // get the MC vertex
  AliGenEventHeader* genHeader = fHeader->GenEventHeader();

  TArrayF vtxMC(3);
  genHeader->PrimaryVertex(vtxMC);

  // ########################################################
  // loop over mc particles
  TTree* particleTree = GetKinematics();
  TParticle* particle = 0;
  particleTree->SetBranchAddress("Particles", &particle);

  Int_t nPrim  = fHeader->GetNprimary();
  Int_t nTotal = fHeader->GetNtrack();

  for (Int_t i_mc = nTotal - nPrim; i_mc < nTotal; ++i_mc)
  {
    particleTree->GetEntry(i_mc);

    if (!particle)
      continue;

    if (IsPrimaryCharged(particle, nPrim) == kFALSE)
      continue;

    fdNdEtaCorrection->FillGene(vtxMC[2], particle->Eta());
  }// end of mc particle

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

    fdNdEtaCorrection->FillMeas(vtxMC[2], particle->Eta());

  } // end of track loop

  return kTRUE;
}

void AlidNdEtaCorrectionSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  AliSelector::SlaveTerminate();

  // Add the histograms to the output on each slave server
  if (!fOutput)
  {
    AliDebug(AliLog::kError, "ERROR: Output list not initialized");
    return;
  }

  fOutput->Add(fdNdEtaCorrection->GetGeneratedHistogram());
  fOutput->Add(fdNdEtaCorrection->GetMeasuredHistogram());
}

void AlidNdEtaCorrectionSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  fdNdEtaCorrectionFinal = new dNdEtaCorrection();
  TH2F* measuredHistogram = dynamic_cast<TH2F*> (fOutput->FindObject("etaVsVtx_meas"));
  TH2F* generatedHistogram = dynamic_cast<TH2F*> (fOutput->FindObject("etaVsVtx_gene"));
  if (!measuredHistogram || !generatedHistogram)
  {
     AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p %p", (void*) generatedHistogram, (void*) measuredHistogram));
    return;
  }
  fdNdEtaCorrectionFinal->SetGeneratedHistogram(generatedHistogram);
  fdNdEtaCorrectionFinal->SetMeasuredHistogram(measuredHistogram);

  fdNdEtaCorrectionFinal->Finish();

  TFile* fout = new TFile("correction_map.root","RECREATE");

  fEsdTrackCuts->SaveHistograms("esd_track_cuts");
  fdNdEtaCorrectionFinal->SaveHistograms();

  fout->Write();
  fout->Close();

  fdNdEtaCorrectionFinal->DrawHistograms();
}
