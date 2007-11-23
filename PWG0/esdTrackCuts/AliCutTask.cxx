#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TFile.h"

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "AliPWG0Helper.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"

#include "AliCutTask.h"

// simple task that runs the esd track cuts to evaluate the basic plots created during the cuts

ClassImp(AliCutTask)

//________________________________________________________________________
AliCutTask::AliCutTask(const char *name) 
  : AliAnalysisTask(name, ""), fESD(0), fTrackCuts(0), fVertex(0), fOutput(0)
{
  // Constructor

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());
}

//________________________________________________________________________
void AliCutTask::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    //tree->SetBranchStatus("*", kFALSE);

    tree->SetBranchStatus("fTracks.*", kTRUE);
    tree->SetBranchStatus("Tracks.*", kTRUE);

    tree->SetBranchStatus("fTriggerMask", kTRUE);
    tree->SetBranchStatus("AliESDHeader", kTRUE);

    tree->SetBranchStatus("fSPDVertex*", kTRUE);
    tree->SetBranchStatus("SPDVertex", kTRUE);
    //tree->SetBranchStatus("fPosition[3]", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliCutTask::CreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutput = new TList;
  fOutput->SetOwner();

  fOutput->Add(fTrackCuts);

  fVertex = new TH1F("fVertex", "fVertex;z vtx (cm);Count", 201, -20, 20);
  fOutput->Add(fVertex);
}

//________________________________________________________________________
void AliCutTask::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  // Post output data.
  PostData(0, fOutput);

  if (!AliPWG0Helper::IsVertexReconstructed(fESD->GetVertex()))
    return;

  Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
  fTrackCuts->CountAcceptedTracks(fESD);

  // get the ESD vertex
  fVertex->Fill(fESD->GetVertex()->GetZv());
}

//________________________________________________________________________
void AliCutTask::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {
    Printf("ERROR: fOutput not available");
    return;
  }

  fTrackCuts = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("AliESDtrackCuts"));
  if (!fTrackCuts) {
    Printf("ERROR: fTrackCuts not available");
    return;
  }

  fVertex = dynamic_cast<TH1F*> (fOutput->FindObject("fVertex"));
  if (!fVertex) {
    Printf("ERROR: fVertex not available");
    return;
  }

  TFile* file = TFile::Open("trackCuts.root", "RECREATE");

  fTrackCuts->SaveHistograms();
  fVertex->Write();

  file->Close();

	fTrackCuts->DrawHistograms();

  new TCanvas;
  fVertex->Draw();
}
