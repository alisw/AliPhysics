/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/* AliAnalysisTaskExampleV.cxx
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Uses only the virtual interfaces so it is general enought to process ESD, AOD
 * **AND** flatESD (online)
 *
 * origin: by Mikolaj Krzewicki, mkrzewic@cern.ch
 * (vaguely based on some other example from ANALYSIS/examples)
 */
#include "AliAnalysisTaskExampleV.h"

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"

#include "AliVEventHandler.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliLog.h"
#include "AliExternalTrackParam.h"
#include "AliESDVertex.h"

ClassImp(AliAnalysisTaskExampleV)

//________________________________________________________________________
AliAnalysisTaskExampleV::AliAnalysisTaskExampleV() // All data members should be initialised here
   :AliAnalysisTask(),
    fOutput(0),
    fHistPt(0), 
    fHistEta(0),
    fSkipExec(kFALSE)
{
    //default ctor, never init memory here
}

//________________________________________________________________________
AliAnalysisTaskExampleV::AliAnalysisTaskExampleV(const char *name) // All data members should be initialised here
   :AliAnalysisTask(name, ""),
    fOutput(0),
    fHistPt(0), 
    fHistEta(0),
    fSkipExec(kFALSE)
{
    // ctor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container
    DefineInput(0, TChain::Class());                                            // for input
    DefineOutput(0, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskExampleV::~AliAnalysisTaskExampleV()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
    if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fOutput;
    }
}

//________________________________________________________________________
void AliAnalysisTaskExampleV::CreateOutputObjects()
{
    // Create histograms
    // Called once
        
    fOutput = new TList();
    fOutput->SetOwner();  // ****IMPORTANT!*******
    
    // Create output objects
    Int_t ptbins = 15;
    Float_t ptlow = 0.1, ptup = 3.1;
    fHistPt = new TH1F("fHistPt", "P_{T} distribution for reconstructed", ptbins, ptlow, ptup);
    fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPt->SetMarkerStyle(kFullCircle);
        
    Int_t etabins = 40;
    Float_t etalow = -2.0, etaup = 2.0;
    fHistEta = new TH1F("fHistEta","#eta distribution for reconstructed",etabins, etalow, etaup);
    fHistEta->GetXaxis()->SetTitle("#eta");
    fHistEta->GetYaxis()->SetTitle("counts");
        
    fOutput->Add(fHistPt);
    fOutput->Add(fHistEta);
    PostData(0, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskExampleV::ConnectInputData(Option_t*)
{
  AliVEventHandler *vH = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!vH) {
    Printf("ERROR: Could not get VEventHandler");
  }
  else {
    fV = vH->GetEvent();
  }
}

//________________________________________________________________________
void AliAnalysisTaskExampleV::Exec(Option_t *) 
{
    // Main loop
    // Called for each event
        
    // get pointer to reconstructed event
    AliVEvent *event = fV;
    if (!event) { AliInfo("ERROR: Could not retrieve event"); return; }

    // access the vertex:
    AliESDVertex vertex;
    event->GetPrimaryVertex(vertex);
    if(!(vertex.GetStatus())) return;

    //maybe do something with the vertex
    //

    // Track loop
    Int_t ntracks = event->GetNumberOfTracks();
    for(Int_t i = 0; i < ntracks; i++) {
        AliVTrack* track = event->GetVTrack(i); // pointer to reconstructed to track          
        if(!track) { 
            AliError(Form("ERROR: Could not retrieve track %d",i)); 
            continue; 
        }

        //pt comes from the track parameters, so we need to get them out first:
        AliExternalTrackParam trackParams;
        track->GetTrackParam(trackParams);
                
        fHistPt->Fill(trackParams.Pt());
        fHistEta->Fill(trackParams.Eta());
    }
}


//________________________________________________________________________
void AliAnalysisTaskExampleV::Terminate(Option_t *) 
{
    // Draw result to screen, or perform fitting, normalizations
    // Called once at the end of the query
    fOutput = dynamic_cast<TList*> (GetOutputData(1));
    if(!fOutput) { AliInfo("ERROR: could not retrieve TList fOutput"); return; }
        
    fHistPt = dynamic_cast<TH1F*> (fOutput->FindObject("fHistPt"));
    if (!fHistPt) { AliInfo("ERROR: could not retrieve fHistPt"); return;}
    fHistEta = dynamic_cast<TH1F*> (fOutput->FindObject("fHistEta"));
    if (!fHistEta) { AliInfo("ERROR: could not retrieve fHistEta"); return;}
        
    TCanvas *c = new TCanvas("AliAnalysisTaskExampleV","P_{T} & #eta",10,10,1020,510);
    c->Divide(2,1);
    c->cd(1)->SetLogy();
    fHistPt->DrawCopy("E");
    c->cd(2);
    fHistEta->DrawCopy("E");
}
