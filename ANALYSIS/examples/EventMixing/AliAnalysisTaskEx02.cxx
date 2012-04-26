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
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
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
     fHistZVertexDiff(0),
     fUseLoopInUserExecMix(kFALSE),
     fUseLoopMixedEvent(kFALSE),
     fUseLoopV0(kFALSE),
     fMainInputHandler(0),
     fMixingInputHandler(0) // The last in the above list should not have a comma after it
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
     fHistZVertexDiff(0),
     fUseLoopInUserExecMix(kFALSE),
     fUseLoopMixedEvent(kFALSE),
     fUseLoopV0(kFALSE),
     fMainInputHandler(0),
     fMixingInputHandler(0) // The last in the above list should not have a comma after it
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


   // sets helper pointers for Mixing
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   SetMainInputHandler(mgr);
   if (fMainInputHandler) SetMixingInputHandler(fMainInputHandler);

   PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskEx02::UserExec(Option_t *)
{
   // Main loop
   // Called for each eventn (Simply standard UserExec() like you would do without mixing)
   // If you want to process events which are used for mixing,
   // then you have to check for (multiplicity, Vz, ... ) ranges set in event pool.
   // This action is up to user
   //


   AliESDEvent *esd = dynamic_cast<AliESDEvent *>(GetMainEvent());
   if (esd) {
      if (!fUseLoopInUserExecMix) {
         if (fUseLoopV0)
            LoopV0(esd);
         else
            Loop(esd);
      }
   } else {
      AliAODEvent *aod = dynamic_cast<AliAODEvent *>(GetMainEvent());
      if (aod) {
         if (!fUseLoopInUserExecMix) {
            if (fUseLoopV0)
               LoopV0(aod);
            else
               Loop(aod);
         }
      }
   }

   PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskEx02::UserExecMix(Option_t *)
{
   //
   // Running mixing function
   //

   if (!fMixingInputHandler) return;

   Int_t bufferSize = fMixingInputHandler->BufferSize();
   Int_t numberMixed = fMixingInputHandler->NumberMixed();
   if (numberMixed == 1) {
      // you may process Main Event here instead of UserExec()
      // Here you will have events which are in range of bins in event pool
      // but only in case when other N (bufferSize) events are found for mix
   }

   AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(GetMainEvent());
   if (esdEvent) {

      for(Int_t iBuffer=0; iBuffer<bufferSize; iBuffer++) {
         AliESDEvent *esdEventMix = dynamic_cast<AliESDEvent *>(GetMixedEvent(iBuffer));
         AliDebug(AliLog::kDebug, Form("Multi=%d MultiMix=%d", esdEvent->GetNumberOfTracks(), esdEventMix->GetNumberOfTracks()));
//          AliDebug(AliLog::kDebug, Form("Mask=%lld MaskMix=%lld", esdEvent->GetTriggerMask(), esdEvent->GetTriggerMask()));
         fHistMultiDiff->Fill(TMath::Abs(esdEvent->GetNumberOfTracks() - esdEventMix->GetNumberOfTracks()));

         // For tesing purpose only (no physics)
         if (fUseLoopInUserExecMix) {
            if (fUseLoopV0) {
               if (fUseLoopMixedEvent) LoopV0(esdEventMix);
               else LoopV0(esdEvent);
            }
            else {
               if (fUseLoopMixedEvent) Loop(esdEventMix);
               else Loop(esdEvent);
            }
         }

      }
   } else {
      AliAODEvent *aodEvent = dynamic_cast<AliAODEvent *>(GetMainEvent());
      if (aodEvent) {
         for(Int_t iBuffer=0; iBuffer<bufferSize; iBuffer++) {
            AliAODEvent *aodEventMix = dynamic_cast<AliAODEvent *>(GetMixedEvent(iBuffer));
            AliDebug(AliLog::kDebug, Form("Multi=%d MultiMix=%d", aodEvent->GetNumberOfTracks(), aodEventMix->GetNumberOfTracks()));
//             AliDebug(AliLog::kDebug, Form("Mask=%lld MaskMix=%lld", aodEvent->GetTriggerMask(), aodEventMix->GetTriggerMask()));
            fHistMultiDiff->Fill(TMath::Abs(aodEvent->GetNumberOfTracks() - aodEventMix->GetNumberOfTracks()));

            // For tesing purpose only (no physics)
            if (fUseLoopInUserExecMix) {
               if (fUseLoopV0) {
                  if (fUseLoopMixedEvent) LoopV0(aodEventMix);
                  else LoopV0(aodEvent);
               }
               else {
                  if (fUseLoopMixedEvent) Loop(aodEventMix);
                  else Loop(aodEvent);
               }
            }
            AliAODVertex *aodVertex = aodEvent->GetPrimaryVertex();
            AliAODVertex *aodVertexMix = aodEventMix->GetPrimaryVertex();
            if ( aodVertex && aodVertexMix) {
               AliDebug(AliLog::kDebug, Form("Vz=%f VzMix=%f", aodVertex->GetZ(), aodVertexMix->GetZ()));
               fHistZVertexDiff->Fill(TMath::Abs( aodVertex->GetZ()-aodVertexMix->GetZ()));
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

//________________________________________________________________________
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

//________________________________________________________________________
void AliAnalysisTaskEx02::LoopV0(AliESDEvent *esd)
{
   //
   // V0 loop for reconstructed event (ESD)
   //

   Int_t nv0 = esd->GetNumberOfV0s();
   Int_t lIndexTrackPos;
   AliESDtrack *myTrackPosTest;
   for (Int_t i = 0; i < nv0; i++)
   {
      AliESDv0 *esdv0 = esd->GetV0(i);
      if (!esdv0) {
         AliError(Form("ERROR: Could not retrieve esdv0 %d", i));
         continue;
      }
      lIndexTrackPos = TMath::Abs(esdv0->GetPindex());
      myTrackPosTest = esd->GetTrack(lIndexTrackPos);
      fHistPt->Fill(myTrackPosTest->Pt());
      fHistEta->Fill(myTrackPosTest->Eta());
   }

}

//________________________________________________________________________
void AliAnalysisTaskEx02::LoopESDMC()
{
   //
   // TODO
   //
}

//________________________________________________________________________
void AliAnalysisTaskEx02::Loop(AliAODEvent *aod)
{
   //
   // Loops over AOD event
   //

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



//________________________________________________________________________
void AliAnalysisTaskEx02::LoopV0(AliAODEvent *aod)
{
   //
   // V0 loop for reconstructed event (AOD)
   //

   Int_t nv0 = aod->GetNumberOfV0s();
   AliAODTrack *postrackmix;
   for (Int_t i = 0; i < nv0; i++)
   {
      AliAODv0 *aodv0 = dynamic_cast<AliAODv0 *>(aod->GetV0(i));
      if (!aodv0) {
         AliError(Form("ERROR: Could not retrieve aodv0 %d", i));
         continue;
      }
      postrackmix = (AliAODTrack *)aodv0->GetDaughter(0);
      fHistPt->Fill(postrackmix->Pt());
      fHistEta->Fill(postrackmix->Eta());
   }

}

//________________________________________________________________________
void AliAnalysisTaskEx02::LoopAODMC()
{
   //
   // TODO
   //
}

//________________________________________________________________________
AliVEvent *AliAnalysisTaskEx02::GetMainEvent()
{
   //
   // Access to MainEvent
   //

   AliMultiInputEventHandler *inEvHMainMulti = fMainInputHandler;
   if (inEvHMainMulti) {
      AliInputEventHandler *inEvMain = dynamic_cast<AliInputEventHandler *>(inEvHMainMulti->GetFirstInputEventHandler());
      if (inEvMain) return inEvMain->GetEvent();
   }

   return 0;
}

//________________________________________________________________________
AliVEvent *AliAnalysisTaskEx02::GetMixedEvent(Int_t buffId)
{
   //
   // Access to Mixed event with buffer id
   //

   AliMultiInputEventHandler *inEvHMain = fMainInputHandler;
   if (inEvHMain) {

      AliMixInputEventHandler *mixIH = fMixingInputHandler;
      if (!mixIH) return 0;
      if (mixIH->CurrentBinIndex() < 0) {
         AliDebug(AliLog::kDebug + 1, "Current event mixEH->CurrentEntry() == -1");
         return 0;
      }
//       AliMultiInputEventHandler *inEvHMixedCurrent = mixEH->GetFirstMultiInputHandler();
      AliMultiInputEventHandler *inEvHMixedCurrent = dynamic_cast<AliMultiInputEventHandler *>(mixIH->InputEventHandler(buffId));
      if (!inEvHMixedCurrent) return 0;
      AliInputEventHandler *ihMixedCurrent = inEvHMixedCurrent->GetFirstInputEventHandler();
      if (ihMixedCurrent) return ihMixedCurrent->GetEvent();
   }

   return 0;
}

//________________________________________________________________________
AliMultiInputEventHandler *AliAnalysisTaskEx02::SetMainInputHandler(AliAnalysisManager *mgr)
{
   //
   // Sets main input handler
   //

   if (!fMainInputHandler) fMainInputHandler = dynamic_cast<AliMultiInputEventHandler *>(mgr->GetInputEventHandler());

   return fMainInputHandler;
}

//________________________________________________________________________
AliMixInputEventHandler *AliAnalysisTaskEx02::SetMixingInputHandler(AliMultiInputEventHandler *mainIH)
{
   //
   // Sets mixing input handler
   //

   if (!fMixingInputHandler) fMixingInputHandler = dynamic_cast<AliMixInputEventHandler *>(mainIH->GetFirstMultiInputHandler());

   return fMixingInputHandler;
}


