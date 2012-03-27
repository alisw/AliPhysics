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

/* $Id: AliAnalysisTaskEx02.cxx 46301 2011-01-06 14:25:27Z agheata $ */

/* AliAnalysisTaskEx02.cxx
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Instructions for adding histograms can be found below, starting with NEW HISTO
 *
 * Based on tutorial example from offline pages using AliMultiInputHandler
 * Edited by Martin Vala
 */
#include <TCanvas.h>
#include <TList.h>
#include <TH1.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"

#include "AliAnalysisTaskEx02.h"
#include <AliMultiInputEventHandler.h>
#include <AliMixInputEventHandler.h>

ClassImp(AliAnalysisTaskEx02)

//________________________________________________________________________
AliAnalysisTaskEx02::AliAnalysisTaskEx02() // All data members should be initialised here
   : AliAnalysisTaskSE(),
     fOutput(0),
     fHistPt(0),
     fHistEta(0),
     fHistMultiDiff(0),
     fHistZVertexDiff(0) // The last in the above list should not have a comma after it
{
   // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskEx02::AliAnalysisTaskEx02(const char *name) // All data members should be initialised here
   : AliAnalysisTaskSE(name),
     fOutput(0),
     fHistPt(0),
     fHistEta(0),
     fHistMultiDiff(0),
     fHistZVertexDiff(0) // The last in the above list should not have a comma after it
{
   // Constructor
   // Define input and output slots here (never in the dummy constructor)
   // Input slot #0 works with a TChain - it is connected to the default input container
   // Output slot #1 writes into a TH1 container
   DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskEx02::~AliAnalysisTaskEx02()
{
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
   if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
   }
}

//________________________________________________________________________
void AliAnalysisTaskEx02::UserCreateOutputObjects()
{
   // Create histograms
   // Called once (on the worker node)

   fOutput = new TList();
   fOutput->SetOwner();  // IMPORTANT!

   // Create histograms
   Int_t ptbins = 15;
   Float_t ptlow = 0.1, ptup = 3.1;
   fHistPt = new TH1F("fHistPt", "P_{T} distribution for reconstructed", ptbins, ptlow, ptup);
   fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
   fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
   fHistPt->SetMarkerStyle(kFullCircle);

   Int_t etabins = 40;
   Float_t etalow = -2.0, etaup = 2.0;
   fHistEta = new TH1F("fHistEta", "#eta distribution for reconstructed", etabins, etalow, etaup);
   fHistEta->GetXaxis()->SetTitle("#eta");
   fHistEta->GetYaxis()->SetTitle("counts");

   fHistMultiDiff = new TH1F("fHistMultiDiff", "Multiplicity Difference", 100, 0, 100);
   fHistMultiDiff->GetXaxis()->SetTitle("MultiDiff");
   fHistMultiDiff->GetYaxis()->SetTitle("counts");

   fHistZVertexDiff = new TH1F("fHistZVertexDiff", "Multiplicity Difference", 100, 0, 10);
   fHistZVertexDiff->GetXaxis()->SetTitle("ZVertexDiff");
   fHistZVertexDiff->GetYaxis()->SetTitle("counts");

   // NEW HISTO should be defined here, with a sensible name,

   fOutput->Add(fHistPt);
   fOutput->Add(fHistEta);
   fOutput->Add(fHistMultiDiff);
   fOutput->Add(fHistZVertexDiff);
   // NEW HISTO added to fOutput here
   PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskEx02::UserExec(Option_t *)
{
   // Main loop
   // Called for each event


   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *inEvHMainMulti = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
   if (inEvHMainMulti) {
      AliInputEventHandler *inEvMain = dynamic_cast<AliInputEventHandler *>(inEvHMainMulti->GetFirstInputEventHandler());

      AliESDEvent *esd = dynamic_cast<AliESDEvent *>(inEvMain->GetEvent());
      if (esd) {
         Loop(esd);
      } else {
         AliAODEvent *aod = dynamic_cast<AliAODEvent *>(inEvMain->GetEvent());
         if (aod) Loop(aod);
      }
   }

   PostData(1, fOutput);
}

void AliAnalysisTaskEx02::UserExecMix(Option_t *)
{

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *inEvHMain = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());
   if (inEvHMain) {

      AliMixInputEventHandler *mixEH = dynamic_cast<AliMixInputEventHandler *>(inEvHMain->GetFirstMultiInputHandler());
      if (!mixEH) return;
      if (mixEH->CurrentBinIndex() < 0) {
         AliDebug(AliLog::kDebug + 1, "Current event mixEH->CurrentEntry() == -1");
         return ;
      }
      AliDebug(AliLog::kDebug, Form("Mixing %lld %d [%lld,%lld] %d", mixEH->CurrentEntry(), mixEH->CurrentBinIndex(), mixEH->CurrentEntryMain(), mixEH->CurrentEntryMix(), mixEH->NumberMixed()));
      AliInputEventHandler *ihMainCurrent = inEvHMain->GetFirstInputEventHandler();
      AliMultiInputEventHandler *inEvHMixedCurrent = mixEH->GetFirstMultiInputHandler();
      AliInputEventHandler *ihMixedCurrent = inEvHMixedCurrent->GetFirstInputEventHandler();
      AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(ihMainCurrent->GetEvent());
      if (esdEvent) {
         AliESDEvent *esdEventMix = dynamic_cast<AliESDEvent *>(ihMixedCurrent->GetEvent());
         AliDebug(AliLog::kDebug, Form("Multi=%d MultiMix=%d", esdEvent->GetNumberOfTracks(), esdEventMix->GetNumberOfTracks()));
//          AliDebug(AliLog::kDebug, Form("Mask=%lld MaskMix=%lld", esdEvent->GetTriggerMask(), esdEvent->GetTriggerMask()));
         fHistMultiDiff->Fill(TMath::Abs(esdEvent->GetNumberOfTracks() - esdEventMix->GetNumberOfTracks()));

      } else {
         AliAODEvent *aodEvent = dynamic_cast<AliAODEvent *>(ihMainCurrent->GetEvent());
         if (aodEvent) {
            AliAODEvent *aodEventMix = dynamic_cast<AliAODEvent *>(ihMixedCurrent->GetEvent());
            AliDebug(AliLog::kDebug, Form("Multi=%d MultiMix=%d", aodEvent->GetNumberOfTracks(), aodEventMix->GetNumberOfTracks()));
//             AliDebug(AliLog::kDebug, Form("Mask=%lld MaskMix=%lld", aodEvent->GetTriggerMask(), aodEventMix->GetTriggerMask()));
            fHistMultiDiff->Fill(TMath::Abs(aodEvent->GetNumberOfTracks() - aodEventMix->GetNumberOfTracks()));

            AliAODVertex *aodVertex = aodEvent->GetPrimaryVertex();
            AliAODVertex *aodVertexMix = aodEventMix->GetPrimaryVertex();
            if ( aodVertex && aodVertexMix) {
               AliDebug(AliLog::kDebug, Form("Vz=%f VzMix=%f", aodVertex->GetZ(), aodVertexMix->GetZ()));
               fHistZVertexDiff->Fill(TMath::Abs( aodVertex->GetZ()-aodVertexMix->GetZ()));
            } else {
               AliWarning("Wrong format: input is not ESD or AOD !!!");
            }

         }
      }
   }

   PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskEx02::Terminate(Option_t *)
{
   // Draw result to screen, or perform fitting, normalizations
   // Called once at the end of the query

   fOutput = dynamic_cast<TList *>(GetOutputData(1));
   if (!fOutput) { AliError("Could not retrieve TList fOutput"); return; }

   fHistPt = dynamic_cast<TH1F *>(fOutput->FindObject("fHistPt"));
   if (!fHistPt) { AliError("Could not retrieve fHistPt"); return;}

   fHistEta = dynamic_cast<TH1F *>(fOutput->FindObject("fHistEta"));
   if (!fHistEta) { AliError("Could not retrieve fHistEta"); return;}

   fHistMultiDiff = dynamic_cast<TH1F *>(fOutput->FindObject("fHistMultiDiff"));
   if (!fHistMultiDiff) { AliError("Could not retrieve fHistMultiDiff"); return;}

   fHistZVertexDiff = dynamic_cast<TH1F *>(fOutput->FindObject("fHistZVertexDiff"));
   if (!fHistZVertexDiff) { AliError("Could not retrieve fHistZVertexDiff"); return;}


   // NEW HISTO should be retrieved from the TList container in the above way,
   // so it is available to draw on a canvas such as below

   TCanvas *c = new TCanvas("AliAnalysisTaskEx02", "P_{T} & #eta", 10, 10, 820, 410);
   c->Divide(2, 2);
   c->cd(1)->SetLogy();
   fHistPt->DrawCopy("E");
   c->cd(2);
   fHistEta->DrawCopy("E");
   c->cd(3)->SetLogy();
   fHistMultiDiff->DrawCopy();
   c->cd(4)->SetLogy();
   fHistZVertexDiff->DrawCopy();


}

void AliAnalysisTaskEx02::Loop(AliESDEvent *esd)
{
   // Track loop for reconstructed event
   Int_t ntracks = esd->GetNumberOfTracks();
   for (Int_t i = 0; i < ntracks; i++) {
      AliESDtrack *esdtrack = esd->GetTrack(i); // pointer to reconstructed to track
      if (!esdtrack) {
         AliError(Form("ERROR: Could not retrieve esdtrack %d", i));
         continue;
      }
      fHistPt->Fill(esdtrack->Pt());
      fHistEta->Fill(esdtrack->Eta());
   }
}
void AliAnalysisTaskEx02::LoopESDMC()
{
   // TODO
}

void AliAnalysisTaskEx02::Loop(AliAODEvent *aod)
{

   // Track loop for reconstructed event
   Int_t ntracks = aod->GetNumberOfTracks();
   for (Int_t i = 0; i < ntracks; i++) {
      AliAODTrack *aodTrack = aod->GetTrack(i); // pointer to reconstructed to track
      if (!aodTrack) {
         AliError(Form("ERROR: Could not retrieve esdtrack %d", i));
         continue;
      }

      fHistPt->Fill(aodTrack->Pt());
      fHistEta->Fill(aodTrack->Eta());
   }
}
void AliAnalysisTaskEx02::LoopAODMC()
{
   // TODO
}
