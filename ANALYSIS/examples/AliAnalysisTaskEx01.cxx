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

/* AliAnalysisTaskEx01.cxx
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Instructions for adding histograms can be found below, starting with NEW HISTO
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 */
#include "AliAnalysisTaskEx01.h"

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"

ClassImp(AliAnalysisTaskEx01)

//________________________________________________________________________
AliAnalysisTaskEx01::AliAnalysisTaskEx01() // All data members should be initialised here
   :AliAnalysisTaskSE(),
    fOutput(0),
    fTrackCuts(0),
    fHistPt(0), 
    fHistEta(0)  // The last in the above list should not have a comma after it
{
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskEx01::AliAnalysisTaskEx01(const char *name) // All data members should be initialised here
   :AliAnalysisTaskSE(name),
    fOutput(0),
    fTrackCuts(0),
    fHistPt(0), 
    fHistEta(0)  // The last in the above list should not have a comma after it
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskEx01::~AliAnalysisTaskEx01()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
    if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fOutput;
    }
    if (fTrackCuts) delete fTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskEx01::UserCreateOutputObjects()
{
    // Create histograms
    // Called once (on the worker node)
        
    fOutput = new TList();
    fOutput->SetOwner();  // IMPORTANT!
    
    fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
    // === Primary Track Selection ===
    //
    // The definition of a primary track is taken from the ALICE Twiki
    // page https://twiki.cern.ch/twiki/bin/view/ALICE/SelectionOfPrimaryTracksForPpDataAnalysis
    // using the following parameters for a standard dN/dPt analysis:
    //  track quality cuts:
    //          esdTrackCuts->SetMinNClustersTPC(70);
    //          esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    //          esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //          esdTrackCuts->SetRequireTPCRefit(kTRUE);
    //          esdTrackCuts->SetRequireITSRefit(kTRUE);
    //          esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
    //                                                                                     AliESDtrackCuts::kAny);
    //  dca:
    //          if(selPrimaries) {
    //                  // 7*(0.0026+0.0050/pt^1.01)
    //                  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    //          }
    //          esdTrackCuts->SetMaxDCAToVertexZ(2);
    //          esdTrackCuts->SetDCAToVertex2D(kFALSE);
    //          esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    //
    // The Primary Track Selection is implemented here by creating an
    // AliESDtrackCuts object, with kTRUE argument to choose primary tracks.
    //
    // By default, it is set to the above conditions which are suitable for
    // a standard inclusive dN/dPt analysis. For others, such as identified
    // dN/dPt or strangeness as well as others, follow the above link for
    // the specific changes to include in the selection.
        
    // To change cuts after selecting some default set, one can use 
    // esdtrackcuts->SetMinNClustersTPC(70) for example

    // Create histograms
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
        
    // NEW HISTO should be defined here, with a sensible name,
        
    fOutput->Add(fHistPt);
    fOutput->Add(fHistEta);
    // NEW HISTO added to fOutput here
    PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskEx01::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
        
        
    // Create pointer to reconstructed event
    return;
    AliVEvent *event = InputEvent();
    if (!event) { Printf("ERROR: Could not retrieve event"); return; }

    // If the task accesses MC info, this can be done as in the commented block below:
    /*
    // Create pointer to reconstructed event
    AliMCEvent *mcEvent = MCEvent();
    if (!mcEvent) { Printf("ERROR: Could not retrieve MC event"); return; }
    Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
         
    // set up a stack for use in check for primary/stable particles
    AliStack* stack = mcEvent->Stack();
    if( !stack ) { Printf( "Stack not available"); return; }
    */  
        
    // create pointer to event
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    if (!esd) {
        AliError("Cannot get the ESD event");
        return;
    }  
//    AliESDHeader* esdheader = (AliESDHeader*)esd->GetHeader();
        
    // === Physics Selection Task ===
    // 
    // To perform a physics selection here, a bitwise operation is used against
    // the UInt_t mask which is extracted in the following way:
    //
    //  UInt_t mask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();   
    //
    // This can be tested to produce the following
    //
    //  Bool_t bMinBias = (mask == AliVEvent::kMB) ? 1 : 0; // check if minimum bias trigger class fired
    //  Bool_t bHighMult = (mask == AliVEvent::kHighMult) ? 1 : 0; // check if high multiplicity trigger class fired
    //
    // For more complicated trigger selections, one can directly test both
    // trigger classes and fired trigger inputs for a particular event, for e.g.
    //
    //  Bool_t bCSH1 = (esd->IsTriggerClassFired("CSH1-B-NOPF-ALLNOTRD")) ? 1 : 0;
    //  Bool_t b0SH1 = (esdheader->IsTriggerInputFired("0SH1")) ? 1 : 0;
    //
    // These booleans can then be used to fill different histograms for specific
    // conditions, or summed to make one cut for all events that fill histos.

    // Do some fast cuts first
    // check for good reconstructed vertex
    if(!(esd->GetPrimaryVertex()->GetStatus())) return;
    // if vertex is from spd vertexZ, require more stringent cut
    if (esd->GetPrimaryVertex()->IsFromVertexerZ()) {
        if (esd->GetPrimaryVertex()->GetDispersion()>0.02 ||  esd->GetPrimaryVertex()->GetZRes()>0.25 ) return; // bad vertex from VertexerZ
    }
        
    // Track loop for reconstructed event
    Int_t ntracks = esd->GetNumberOfTracks();
    for(Int_t i = 0; i < ntracks; i++) {
        AliESDtrack* esdtrack = esd->GetTrack(i); // pointer to reconstructed to track          
        if(!esdtrack) { 
            AliError(Form("ERROR: Could not retrieve esdtrack %d",i)); 
            continue; 
        }
                
        // Some MC checks, if MC is used
        //if(esdtrack->GetLabel() < 0) continue; // get rid of "ghost" tracks
                
        // ... and the thorough checking of ESD cuts after.
        // if this is not a primary track, skip to the next one
        if(!fTrackCuts->AcceptTrack(esdtrack)) continue;
                
        fHistPt->Fill(esdtrack->Pt());
        fHistEta->Fill(esdtrack->Eta());
    }
    // NEW HISTO should be filled before this point, as PostData puts the
    // information for this iteration of the UserExec in the container
    PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskEx01::Terminate(Option_t *) 
{
    // Draw result to screen, or perform fitting, normalizations
    // Called once at the end of the query
    return;    
    fOutput = dynamic_cast<TList*> (GetOutputData(1));
    if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
        
    fHistPt = dynamic_cast<TH1F*> (fOutput->FindObject("fHistPt"));
    if (!fHistPt) { Printf("ERROR: could not retrieve fHistPt"); return;}
    fHistEta = dynamic_cast<TH1F*> (fOutput->FindObject("fHistEta"));
    if (!fHistEta) { Printf("ERROR: could not retrieve fHistEta"); return;}
        
    // Get the physics selection histograms with the selection statistics
    //AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    //AliESDInputHandler *inputH = dynamic_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
    //TH2F *histStat = (TH2F*)inputH->GetStatistics();
   
   
    // NEW HISTO should be retrieved from the TList container in the above way,
    // so it is available to draw on a canvas such as below

    TCanvas *c = new TCanvas("AliAnalysisTaskEx01","P_{T} & #eta",10,10,1020,510);
    c->Divide(2,1);
    c->cd(1)->SetLogy();
    fHistPt->DrawCopy("E");
    c->cd(2);
    fHistEta->DrawCopy("E");
}
