/* $Id$ */

#include "AlidNdEtaAnalysisESDSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

#include <AliLog.h>
#include <AliESDVertex.h>
#include <AliESD.h>

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "dNdEta/dNdEtaAnalysis.h"
#include "AliPWG0Helper.h"

#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "TParticle.h"

ClassImp(AlidNdEtaAnalysisESDSelector)

AlidNdEtaAnalysisESDSelector::AlidNdEtaAnalysisESDSelector() :
  AliSelector(),
  fdNdEtaAnalysis(0),
  fMult(0),
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
}

void AlidNdEtaAnalysisESDSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliSelector::SlaveBegin(tree);

  ReadUserObjects(tree);

  fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fMult = new TH1F("fMult", "fMult;Ntracks;Count", 201, -0.5, 200.5);
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
    tree->SetBranchStatus("fTracks.fLabel", 1);

    AliESDtrackCuts::EnableNeededBranches(tree);
  }
}

Bool_t AlidNdEtaAnalysisESDSelector::Process(Long64_t entry)
{
  // loop over all events

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

  // get the ESD vertex
  const AliESDVertex* vtxESD = fESD->GetVertex();
  Double_t vtx[3];
  vtxESD->GetXYZ(vtx);

  // get number of "good" tracks
  TObjArray* list = fEsdTrackCuts->GetAcceptedTracks(fESD);
  Int_t nGoodTracks = list->GetEntries();

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

    fdNdEtaAnalysis->FillTrack(vtx[2], eta, pt);
  } // end of track loop

  delete list;
  list = 0;

  // for event count per vertex
  fdNdEtaAnalysis->FillEvent(vtx[2], nGoodTracks);
  fMult->Fill(nGoodTracks);

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
  fOutput->Add(fMult);

  fdNdEtaAnalysis = 0;
}

void AlidNdEtaAnalysisESDSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliSelector::Terminate();

  fdNdEtaAnalysis = dynamic_cast<dNdEtaAnalysis*> (fOutput->FindObject("dndeta"));
  fMult = dynamic_cast<TH1F*> (fOutput->FindObject("fMult"));

  if (!fdNdEtaAnalysis)
  {
    AliDebug(AliLog::kError, "ERROR: Histograms not available");
    return;
  }

  TFile* fout = new TFile("analysis_esd_raw.root", "RECREATE");

  if (fdNdEtaAnalysis)
    fdNdEtaAnalysis->SaveHistograms();

  if (fEsdTrackCuts)
    fEsdTrackCuts->SaveHistograms("esd_tracks_cuts");

  if (fMult)
    fMult->Write();

  fout->Write();
  fout->Close();
}
