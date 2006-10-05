/* $Id$ */

#include "AlidNdEtaAnalysisESDSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliESDVertex.h>
#include <AliESD.h>

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "dNdEta/dNdEtaAnalysis.h"
#include "dNdEta/AlidNdEtaCorrection.h"
#include "AliPWG0Helper.h"

ClassImp(AlidNdEtaAnalysisESDSelector)

AlidNdEtaAnalysisESDSelector::AlidNdEtaAnalysisESDSelector() :
  AliSelector(),
  fdNdEtaAnalysisMBVtx(0),
  fdNdEtaAnalysisMB(0),
  fdNdEtaAnalysis(0),
  fEsdTrackCuts(0),
  fdNdEtaCorrection(0)
  {
  //
  // Constructor. Initialization of pointers
  //

  AliLog::SetClassDebugLevel("AlidNdEtaAnalysisESDSelector", AliLog::kDebug);
}

AlidNdEtaAnalysisESDSelector::~AlidNdEtaAnalysisESDSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AlidNdEtaAnalysisESDSelector::Begin(TTree* tree)
{
  // Begin function

  ReadUserObjects(tree);
}

void AlidNdEtaAnalysisESDSelector::ReadUserObjects(TTree* tree)
{
  // read the user objects, called from slavebegin and begin

  if (!fEsdTrackCuts && fInput)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts && tree)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (tree->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from input list.");


  if (!fdNdEtaCorrection && fInput)
    fdNdEtaCorrection = dynamic_cast<AlidNdEtaCorrection*> (fInput->FindObject("dndeta_correction"));

  if (!fdNdEtaCorrection && tree)
    fdNdEtaCorrection = dynamic_cast<AlidNdEtaCorrection*> (tree->GetUserInfo()->FindObject("dndeta_correction"));

  if (!fdNdEtaCorrection)
     AliDebug(AliLog::kError, "ERROR: Could not read dndeta_correction from input list.");
}

void AlidNdEtaAnalysisESDSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  ReadUserObjects(tree);

  fdNdEtaAnalysisMBVtx = new dNdEtaAnalysis("dndeta_mbvtx", "dndeta_mbvtx");
  fdNdEtaAnalysisMB = new dNdEtaAnalysis("dndeta_mb", "dndeta_mb");
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
}

void AlidNdEtaAnalysisESDSelector::Init(TTree* tree)
{
  // read the user objects

  AliSelector::Init(tree);

  // Enable only the needed branches
  if (tree)
  {
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("fTriggerMask", 1);
    tree->SetBranchStatus("fSPDVertex*", 1);

    AliESDtrackCuts::EnableNeededBranches(tree);
  }
}

Bool_t AlidNdEtaAnalysisESDSelector::Process(Long64_t entry)
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

  if (AliSelector::Process(entry) == kFALSE)
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

  if (!fdNdEtaCorrection)
  {
    AliDebug(AliLog::kError, "fdNdEtaCorrection not available");
    return kFALSE;
  }

  if (AliPWG0Helper::IsEventTriggered(fESD) == kFALSE)
  {
    AliDebug(AliLog::kDebug+1, Form("Skipping event %d because it was not triggered", (Int_t) entry));
    return kTRUE;
  }

  if (AliPWG0Helper::IsVertexReconstructed(fESD) == kFALSE)
  {
    AliDebug(AliLog::kDebug+1, Form("Skipping event %d because its vertex was not reconstructed", (Int_t) entry));
    return kTRUE;
  }

  // ########################################################
  // get the EDS vertex
  const AliESDVertex* vtxESD = fESD->GetVertex();
  Double_t vtx[3];
  vtxESD->GetXYZ(vtx);

  // get number of "good" tracks
  TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD);
  Int_t nGoodTracks = list->GetEntries();

  Float_t vertexRecoCorr = fdNdEtaCorrection->GetVertexRecoCorrection(vtx[2], nGoodTracks);
  if (vertexRecoCorr <= 0)
  {
    AliDebug(AliLog::kError, Form("INFO: Skipping event because vertexRecoCorr is <= 0 (%f)", vertexRecoCorr));
    delete list;
    return kTRUE;
  }

  Float_t triggerCorr = fdNdEtaCorrection->GetTriggerBiasCorrection(vtx[2], nGoodTracks);
  if (triggerCorr <= 0)
  {
    AliDebug(AliLog::kError, Form("INFO: Skipping event because triggerCorr is <= 0 (%f)", triggerCorr));
    delete list;
    return kTRUE;
  }

  // loop over esd tracks
  for (Int_t t=0; t<nGoodTracks; t++)
  {
    AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (list->At(t));
    if (!esdTrack)
    {
      AliDebug(AliLog::kError, Form("ERROR: Could not retrieve track %d.", t));
      continue;
    }

    Double_t p[3];
    esdTrack->GetConstrainedPxPyPz(p); // ### TODO should be okay because we have a vertex, however GetInnerPxPyPy / GetOuterPxPyPy also exist
    TVector3 vector(p);

    Float_t theta = vector.Theta();
    Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));
    Float_t pt = vector.Pt();

    Float_t track2particleCorr = fdNdEtaCorrection->GetTrack2ParticleCorrection(vtx[2], eta, pt);

    Float_t weight = track2particleCorr * vertexRecoCorr * triggerCorr;
    if (weight <= 0)
    {
      AliDebug(AliLog::kError, Form("INFO: Skipping track because weight is <= 0 (track %d, weight %f) (vtx %f, eta %f, pt %f)", t, weight, vtx[2], eta, pt));
      continue;
    }

    fdNdEtaAnalysisMBVtx->FillTrack(vtx[2], eta, pt, track2particleCorr);
    fdNdEtaAnalysisMB->FillTrack(vtx[2], eta, pt, track2particleCorr * vertexRecoCorr);
    fdNdEtaAnalysis->FillTrack(vtx[2], eta, pt, weight);
  } // end of track loop

  delete list;
  list = 0;

  // for event count per vertex
  fdNdEtaAnalysisMBVtx->FillEvent(vtx[2], 1);
  fdNdEtaAnalysisMB->FillEvent(vtx[2], vertexRecoCorr);
  fdNdEtaAnalysis->FillEvent(vtx[2], vertexRecoCorr * triggerCorr);

  return kTRUE;
}

void AlidNdEtaAnalysisESDSelector::SlaveTerminate()
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

  // Add the objects to the output list and set them to 0, so that the destructor does not delete them.

  fOutput->Add(fdNdEtaAnalysis);
  fdNdEtaAnalysis = 0;

  fOutput->Add(fdNdEtaAnalysisMB);
  fdNdEtaAnalysisMB = 0;

  fOutput->Add(fdNdEtaAnalysisMBVtx);
  fdNdEtaAnalysisMBVtx = 0;
}

void AlidNdEtaAnalysisESDSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  fdNdEtaAnalysis = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta"));
  fdNdEtaAnalysisMB = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta_mb"));
  fdNdEtaAnalysisMBVtx = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta_mbvtx"));

  if (!fdNdEtaAnalysis || !fdNdEtaAnalysisMB || !fdNdEtaAnalysisMBVtx)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p %p %p", (void*) fdNdEtaAnalysis, (void*) fdNdEtaAnalysisMB, (void*) fdNdEtaAnalysisMBVtx));
    return;
  }

  if (fdNdEtaAnalysis)
    fdNdEtaAnalysis->Finish(fdNdEtaCorrection, 0.3);

  if (fdNdEtaAnalysisMB)
    fdNdEtaAnalysisMB->Finish(fdNdEtaCorrection, 0.3);

  if (fdNdEtaAnalysisMBVtx)
    fdNdEtaAnalysisMBVtx->Finish(fdNdEtaCorrection, 0.3);

  TFile* fout = new TFile("analysis_esd.root","RECREATE");

  if (fdNdEtaAnalysis)
    fdNdEtaAnalysis->SaveHistograms();

  if (fdNdEtaAnalysisMB)
    fdNdEtaAnalysisMB->SaveHistograms();

  if (fdNdEtaAnalysisMBVtx)
    fdNdEtaAnalysisMBVtx->SaveHistograms();

  if (fEsdTrackCuts)
    fEsdTrackCuts->SaveHistograms("esd_tracks_cuts");

  if (fdNdEtaCorrection)
    fdNdEtaCorrection->SaveHistograms();

  fout->Write();
  fout->Close();
}
