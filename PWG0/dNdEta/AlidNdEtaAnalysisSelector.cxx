#include "AlidNdEtaAnalysisSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>

#include <AliLog.h>
#include <AliGenEventHeader.h>
#include <AliTracker.h>

#include "../esdTrackCuts/AliESDtrackCuts.h"
#include "dNdEtaCorrection.h"
#include "dNdEtaAnalysis.h"

ClassImp(AlidNdEtaAnalysisSelector)

AlidNdEtaAnalysisSelector::AlidNdEtaAnalysisSelector(TTree *) :
  AliSelector(),
  fEsdTrackCuts(0),
  fdNdEtaAnalysis(0),
  fdNdEtaCorrection(0),
  fdNdEtaAnalysisFinal(0)
{
  //
  // Constructor. Initialization of pointers
  //
}

AlidNdEtaAnalysisSelector::~AlidNdEtaAnalysisSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AlidNdEtaAnalysisSelector::SlaveBegin(TTree * tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta");

  if (fChain)
  {
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fChain->GetUserInfo()->FindObject("AliESDtrackCuts"));
    fdNdEtaCorrection = dynamic_cast<dNdEtaCorrection*> (fChain->GetUserInfo()->FindObject("dNdEtaCorrection"));
  }

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from user info.");

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kWarning, "ERROR: Could not read dNdEtaCorrection from user info.");

  AliLog::SetClassDebugLevel("AliESDtrackCuts", 1);
}

Bool_t AlidNdEtaAnalysisSelector::Process(Long64_t entry)
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

  // Check prerequisites
  if (!fESD || !fEsdTrackCuts)
    return kFALSE;

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

  Double_t vtx[3];
  vtxESD->GetXYZ(vtx);

  // ########################################################
  // loop over esd tracks
  Int_t nTracks = fESD->GetNumberOfTracks();
  for (Int_t t=0; t<nTracks; t++)
  {
    AliESDtrack* esdTrack = fESD->GetTrack(t);
    if (!esdTrack)
    {
      AliDebug(AliLog::kError, Form("ERROR: Could not retrieve track %d.", t));
      continue;
    }

    // cut the esd track?
    if (!fEsdTrackCuts->AcceptTrack(esdTrack))
      continue;

    Double_t p[3];
    esdTrack->GetConstrainedPxPyPz(p); // ### TODO or GetInnerPxPyPy / GetOuterPxPyPy
    TVector3 vector(p);

    Float_t theta = vector.Theta();
    Float_t eta   = -TMath::Log(TMath::Tan(theta/2.));

    Float_t correction = fdNdEtaCorrection->GetCorrection(vtx[2], eta);

    fdNdEtaAnalysis->FillTrack(vtx[2], eta, correction);

  } // end of track loop

  // for event count per vertex
  fdNdEtaAnalysis->FillEvent(vtx[2]);

  return kTRUE;
}

void AlidNdEtaAnalysisSelector::SlaveTerminate()
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

  fOutput->Add(fdNdEtaAnalysis->GetEtaVsVtxHistogram());
  fOutput->Add(fdNdEtaAnalysis->GetEtaVsVtxUncorrectedHistogram());
  fOutput->Add(fdNdEtaAnalysis->GetVtxHistogram());

  fdNdEtaAnalysis->GetVtxHistogram()->Print();
  fOutput->Print();
}

void AlidNdEtaAnalysisSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  TH2F* etaVsVtxHistogram = dynamic_cast<TH2F*> (fOutput->FindObject("eta_vs_vtx"));
  TH2F* etaVsVtxUncorrectedHistogram = dynamic_cast<TH2F*> (fOutput->FindObject("eta_vs_vtx_uncorrected"));
  TH1D* vtxHistogram = dynamic_cast<TH1D*> (fOutput->FindObject("vtx"));

  if (!etaVsVtxHistogram || !vtxHistogram || !etaVsVtxUncorrectedHistogram)
  {
     AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p %p %p", (void*) etaVsVtxHistogram, (void*) etaVsVtxUncorrectedHistogram, (void*) vtxHistogram));
    return;
  }

  fdNdEtaAnalysisFinal = new dNdEtaAnalysis("dNdEtaResult");

  fdNdEtaAnalysisFinal->SetEtaVsVtxHistogram(etaVsVtxHistogram);
  fdNdEtaAnalysisFinal->SetEtaVsVtxUncorrectedHistogram(etaVsVtxUncorrectedHistogram);
  fdNdEtaAnalysisFinal->SetVtxHistogram(vtxHistogram);

  fdNdEtaAnalysisFinal->Finish();

  TFile* fout = new TFile("out.root","RECREATE");

  if (fEsdTrackCuts)
    fEsdTrackCuts->SaveHistograms("esd_tracks_cuts");

  if (fdNdEtaCorrection)
    fdNdEtaCorrection->SaveHistograms();

  fdNdEtaAnalysisFinal->SaveHistograms();

  fout->Write();
  fout->Close();

  fdNdEtaAnalysisFinal->DrawHistograms();
}
