#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TFile.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#include "esdTrackCuts/AliESDtrackCuts.h"
#include "AliPWG0Helper.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"

#include "AliCutTask.h"

// simple task that runs the esd track cuts to evaluate the basic plots created during the cuts

ClassImp(AliCutTask)

//________________________________________________________________________
AliCutTask::AliCutTask(const char *name) : 
  AliAnalysisTask(name, ""), fESD(0), fTrackCuts(0), 
  fAnalysisMode(AliPWG0Helper::kTPCITS),
  fTrackCutsPrimaries(0),
  fTrackCutsSecondaries(0), fVertex(0), fTriggerStats(0), fOutput(0),
  fPrimStats(0)
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

    // old esd
    tree->SetBranchStatus("fTriggerMask", kTRUE);
    tree->SetBranchStatus("fTracks.*", kTRUE);
    tree->SetBranchStatus("fSPDVertex*", kTRUE);

    // new esd
    tree->SetBranchStatus("TriggerMask", kTRUE);
    tree->SetBranchStatus("AliESDHeader", kTRUE);
    //AliPWG0Helper::SetBranchStatusRecursive(tree, "Tracks", kTRUE);
    //AliPWG0Helper::SetBranchStatusRecursive(tree, "SPDVertex", kTRUE);

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
  if (fTrackCutsPrimaries)
    fOutput->Add(fTrackCutsPrimaries);
  if (fTrackCutsSecondaries)
    fOutput->Add(fTrackCutsSecondaries);

  fVertex = new TH1F("fVertex", "fVertex;z vtx (cm);Count", 201, -20, 20);
  fOutput->Add(fVertex);

  fTriggerStats = new TH1F("fTriggerStats", "fTriggerStats;trigger;Count", 5, -0.5, 4.5);
  fTriggerStats->GetXaxis()->SetBinLabel(1, "!MB1 & !MB2");
  fTriggerStats->GetXaxis()->SetBinLabel(2, "MB1");
  fTriggerStats->GetXaxis()->SetBinLabel(3, "MB2");
  fTriggerStats->GetXaxis()->SetBinLabel(4, "ITS_SPD_GFO_L0");
  fTriggerStats->GetXaxis()->SetBinLabel(5, "VZERO_OR_LEFT | VZERO_OR_RIGHT");
  fOutput->Add(fTriggerStats);

  fPrimStats = new TH1F("fPrimStats", "fPrimStats", 8, 0.5, 8.5);
  fPrimStats->GetXaxis()->SetBinLabel(1, "Primary accepted once");
  fPrimStats->GetXaxis()->SetBinLabel(2, "Primary accepted more than once");
  fPrimStats->GetXaxis()->SetBinLabel(3, "Primary accepted more than once (count)");
  fPrimStats->GetXaxis()->SetBinLabel(4, "Primary track rejected");
  fPrimStats->GetXaxis()->SetBinLabel(5, "Primary track rejected (count)");
  fPrimStats->GetXaxis()->SetBinLabel(6, "Primary track rejected, but other track accepted");
  fPrimStats->GetXaxis()->SetBinLabel(7, "Primary track rejected, but other track accepted (count)");
  fOutput->Add(fPrimStats);

  // disable info messages of AliMCEvent (per event)
  AliLog::SetClassDebugLevel("AliMCEvent", AliLog::kWarning - AliLog::kDebug + 1);
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

  //fESD->GetVertex()->Print();

  if (!AliPWG0Helper::IsEventTriggered(fESD->GetTriggerMask(), AliPWG0Helper::kMB1) && !AliPWG0Helper::IsEventTriggered(fESD->GetTriggerMask(), AliPWG0Helper::kMB2))
    fTriggerStats->Fill(0);

  if (AliPWG0Helper::IsEventTriggered(fESD->GetTriggerMask(), AliPWG0Helper::kMB1))
    fTriggerStats->Fill(1);

  if (AliPWG0Helper::IsEventTriggered(fESD->GetTriggerMask(), AliPWG0Helper::kMB2))
    fTriggerStats->Fill(2);

  if (fESD->GetTriggerMask() & (0x1 << 14))
    fTriggerStats->Fill(3);

  if (fESD->GetTriggerMask() & 1 || fESD->GetTriggerMask() & 2)
    fTriggerStats->Fill(4);

  //if (fESD->GetTriggerMask() & (0x1 << 14) == 0)
  if (!AliPWG0Helper::IsEventTriggered(fESD->GetTriggerMask(), AliPWG0Helper::kMB2))
    return;

  // get the ESD vertex
  const AliESDVertex* vtxESD = AliPWG0Helper::GetVertex(fESD, fAnalysisMode);
  
  if (!vtxESD)
  {
    Printf("No vertex. Skipping event");
    return;
  }

  Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
  Int_t acceptedTracks = fTrackCuts->CountAcceptedTracks(fESD);
  Printf("Accepted %d tracks", acceptedTracks);
  
  if (fTrackCutsPrimaries && fTrackCutsSecondaries)
  {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    AliStack* stack = mcEvent->Stack();
    if (!stack)
    {
      Printf("ERROR: Stack not available");
      return;
    }

    // count if primaries are counted several times
    Int_t max = stack->GetNprimary();
    Int_t* primAcc = new Int_t[max];
    Int_t* primRej = new Int_t[max];
    for (Int_t i=0; i<max; i++)
    {
      primAcc[i] = 0;
      primRej[i] = 0;
    }

    for (Int_t trackId = 0; trackId < fESD->GetNumberOfTracks(); trackId++)
    {
      AliESDtrack* track = fESD->GetTrack(trackId);
      if (!track)
      {
        Printf("UNEXPECTED: Could not get track with id %d", trackId);
        continue;
      }
      
      if (track->GetLabel() == 0)
      {
        Printf("Track %d has no label. Skipping", trackId);
        continue;
      }
      
      Int_t label = TMath::Abs(track->GetLabel());
      if (stack->IsPhysicalPrimary(label) == kTRUE)
      {
	if (label >= max)
	{
	  Printf("Warning label %d is higher than number of primaries %d", label, max);
	  continue;
	}
	    
        if (fTrackCutsPrimaries->AcceptTrack(track)) 
	{
	  primAcc[label]++;
	} 
	else
	  primRej[label]++;

      }
      else
      {
        if (!stack->Particle(label)) {
          Printf("ERROR: Could not get particle with label %d", label);
          continue;
        }
        
        // Deuteron is ok, but missing in Pdg particles in root
        if (stack->Particle(label)->GetPdgCode() != 10010020) {
          if (!stack->Particle(label)->GetPDG()) {
            Printf("ERROR: Could not get PDG particle of particle with label %d (pdg code is %d)", label, stack->Particle(label)->GetPdgCode());
            stack->Particle(label)->Print();
            continue;
          }
          
          if (stack->Particle(label)->GetPDG()->Charge() == 0) {
            Printf("Particle %d has 0 charge. Skipping.", label);
            continue;
          }
        }
          
        fTrackCutsSecondaries->AcceptTrack(track);
      }
    }

    for (Int_t i=0; i<max; i++) {
      if (primAcc[i] == 1) {
	fPrimStats->Fill(1);
      } else if (primAcc[i] > 1) {
	fPrimStats->Fill(2);
	fPrimStats->Fill(3, primAcc[i]);
      }

      if (primRej[i] > 0) {
	if (primAcc[i] == 0) {
	  fPrimStats->Fill(4);
	  fPrimStats->Fill(5, primRej[i]);
	} else if (primAcc[i] > 0) {
	  fPrimStats->Fill(6);
	  fPrimStats->Fill(7, primRej[i]);
	}
      }
    }

    delete[] primAcc;
    delete[] primRej;
  }

  // get the ESD vertex
  fVertex->Fill(vtxESD->GetZv());
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

  fTriggerStats = dynamic_cast<TH1F*> (fOutput->FindObject("fTriggerStats"));
  if (!fTriggerStats) {
    Printf("ERROR: fTriggerStats not available");
    return;
  }

  fTrackCutsPrimaries = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fTrackCutsPrimaries"));
  fTrackCutsSecondaries = dynamic_cast<AliESDtrackCuts*> (fOutput->FindObject("fTrackCutsSecondaries"));
  
  TFile* file = TFile::Open("trackCuts.root", "RECREATE");

  fTrackCuts->SaveHistograms();
  fVertex->Write();
  fTriggerStats->Write();
  
  if (fTrackCutsPrimaries)
    fTrackCutsPrimaries->SaveHistograms();

  if (fTrackCutsSecondaries)
    fTrackCutsSecondaries->SaveHistograms();

  if (fPrimStats)
    fPrimStats->Write();
  
  file->Close();
  
  Printf("Writting results to trackCuts.root.");
  
  fTrackCuts->DrawHistograms();

  new TCanvas;
  fVertex->Draw();

  new TCanvas;
  fTriggerStats->Draw();
}

void AliCutTask::EnableSecondaryStudy()
{ 
  // 
  // clones the fTrackCuts
  //
  
  if (!fTrackCuts) {
    Printf("ERROR: fTrackCuts 0. Do not call this function before SetTrackCuts");
    return;
  }
  
  fTrackCutsPrimaries = dynamic_cast<AliESDtrackCuts*> (fTrackCuts->Clone("fTrackCutsPrimaries"));
  fTrackCutsSecondaries = dynamic_cast<AliESDtrackCuts*> (fTrackCuts->Clone("fTrackCutsSecondaries"));
}
