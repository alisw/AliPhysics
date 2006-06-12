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

ClassImp(AlidNdEtaAnalysisESDSelector)

AlidNdEtaAnalysisESDSelector::AlidNdEtaAnalysisESDSelector() :
  AliSelector(),
  fEsdTrackCuts(0)
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

void AlidNdEtaAnalysisESDSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  if (fInput)
  {
    printf("Printing input list:\n");
    fInput->Print();
  }

  if (!fEsdTrackCuts && fInput)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fInput->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from input list.");

  fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
}

void AlidNdEtaAnalysisESDSelector::Init(TTree* tree)
{
  // read the user objects

  AliSelector::Init(tree);

  if (!fEsdTrackCuts && fTree)
    fEsdTrackCuts = dynamic_cast<AliESDtrackCuts*> (fTree->GetUserInfo()->FindObject("AliESDtrackCuts"));

  if (!fEsdTrackCuts)
     AliDebug(AliLog::kError, "ERROR: Could not read EsdTrackCuts from user info.");
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

    fdNdEtaAnalysis->FillTrack(vtx[2], eta);

  } // end of track loop

  // for event count per vertex
  fdNdEtaAnalysis->FillEvent(vtx[2]);

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

  fOutput->Add(fdNdEtaAnalysis);
}

void AlidNdEtaAnalysisESDSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  fdNdEtaAnalysis = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta"));

  if (!fdNdEtaAnalysis)
  {
    AliDebug(AliLog::kError, Form("ERROR: Histograms not available %p", (void*) fdNdEtaAnalysis));
    return;
  }

  TFile* fout = new TFile("analysis_esd.root","RECREATE");

  if (fdNdEtaAnalysis)
    fdNdEtaAnalysis->SaveHistograms();

  if (fEsdTrackCuts)
    fEsdTrackCuts->SaveHistograms("esd_tracks_cuts");

  fout->Write();
  fout->Close();
}
