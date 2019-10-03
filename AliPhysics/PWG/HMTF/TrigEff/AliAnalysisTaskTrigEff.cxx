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

/* AliAnalysisTaskTrigEff.cxx
 * Simple task to compute INT7 trigger efficiency in the INEL>0 sample, based on a more open trigger
 */
#include "AliAnalysisTaskTrigEff.h"

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
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliMultSelection.h"
#include "AliTriggerAnalysis.h"

ClassImp(AliAnalysisTaskTrigEff)


//________________________________________________________________________
AliAnalysisTaskTrigEff::AliAnalysisTaskTrigEff() // All data members should be initialised here
:AliAnalysisTaskSE(),
  fOutput(0),
  fHistEvCount(0),
  fTriggerAnalysis(0),
  fNtrkCut(1)// The last in the above list should not have a comma after it
{
  // Dummy constructor ALWAYS needed for I/O.
  for(Int_t ihist = 0; ihist < kNHist; ihist++){
    fHistV0M[ihist]       = 0;
    fHistNtrk70100[ihist] = 0;
    
  }
  
}

//________________________________________________________________________
AliAnalysisTaskTrigEff::AliAnalysisTaskTrigEff(const char *name) // All data members should be initialised here
  :AliAnalysisTaskSE(name),
   fOutput(0),
   fHistEvCount(0),
   fTriggerAnalysis(0),
   fNtrkCut(1)// The last in the above list should not have a comma after it
{
  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  for(Int_t ihist = 0; ihist < kNHist; ihist++){
    fHistV0M[ihist]       = 0;
    fHistNtrk70100[ihist] = 0;    
  }

  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());                                            // for output list
  
}

//________________________________________________________________________
AliAnalysisTaskTrigEff::~AliAnalysisTaskTrigEff()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
  }
}

//________________________________________________________________________
void AliAnalysisTaskTrigEff::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)
        
  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!

  fTriggerAnalysis = new AliTriggerAnalysis();

  // Create histograms
  const TString histoNames[] = {"kHistoINT11", "kHistoINT7", "kHistoINT7Offline", "kHistoINT11EvSel", "kHistoINT7EvSel", "kHistoINT7OfflineEvSel"};
  for(Int_t ihist = 0; ihist < kNHist; ihist++){
    //    fHistV0M[ihist] = new TH1F(histoNames[ihist], histoNames[ihist], 110, -0.5, 109.5);
    fHistV0M[ihist] = new TH1F(histoNames[ihist], histoNames[ihist], 110, 0, 110);
    fHistV0M[ihist] -> Sumw2();
    fOutput->Add(fHistV0M[ihist]);
    // The efficiency is very different from 1 in thre 70-100% V0M multiplicity class: we look at the details in the trigger class using the tracklets distribution
    fHistNtrk70100[ihist] = new TH1F(histoNames[ihist]+"_ntrk70100", histoNames[ihist]+"_ntrk70100", 50, -0.5, 49.5);
    fHistNtrk70100[ihist]->Sumw2();
    fOutput->Add(fHistNtrk70100[ihist]);
    
  }
  
        
  fHistEvCount  = new TH1I ("fHistEvCount", "event counter", kNEvSel, -0.5, kNEvSel-0.5);

  fHistEvCount->GetXaxis()->SetBinLabel(fHistEvCount->GetXaxis()->FindBin(kEvAll),      "All");
  fHistEvCount->GetXaxis()->SetBinLabel(fHistEvCount->GetXaxis()->FindBin(kEvINT11),    "INT11");
  fHistEvCount->GetXaxis()->SetBinLabel(fHistEvCount->GetXaxis()->FindBin(kEvFOOnline), "FO Online");
  fHistEvCount->GetXaxis()->SetBinLabel(fHistEvCount->GetXaxis()->FindBin(kEvInelGT0),  Form("INEL > %d (trk)", fNtrkCut));

  fOutput->Add(fHistEvCount);

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskTrigEff::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  // Reference is INT11 + FO offline. On those events we want to check how many have int7, as a function of V0 centrality.
  
        
  // Create pointer to reconstructed event
  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }
  // create pointer to event (change to ESD if needed)
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  if (!esd) {
    AliFatal("Cannot get the ESD event");
    return;
  }  
  // WARNING: INT11 needs to be hooked to INT1 with a custom OADB!!
  Bool_t isINT11 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT1);
  Bool_t isINT7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
  

  // Get Helper Classes
  AliMultSelection *multSelection = (AliMultSelection*) event->FindListObject("MultSelection");
  Bool_t isEvSelCentralityOk = (multSelection->GetEvSelCode()) == 0;
  const AliMultiplicity* mult = esd->GetMultiplicity();
  if (!mult){
    AliFatal("No multiplicity object"); 
  }

  // Histos:
  // - Centrality distribution (int5/int7 with/without isEvSelCentralityOk)

  Float_t v0mCentr = multSelection->GetMultiplicityPercentile("V0M");


  // FASTOR Online (only layer 1)
  TBits fastorOnlBit  = mult->GetFastOrFiredChips();
  TBits fastorOfflBit =  mult->GetFiredChipMap();
  UInt_t fastorOnl1 = fastorOnlBit.CountBits(400);

  // FASTOR offline (only layer 1)
  UInt_t fastorOffl1 = mult->GetNumberOfFiredChips(1);


  // ================= TRIGGER SELECTION ==================
  fHistEvCount->Fill(kEvAll);
  // Basic trigger selection
  if(!isINT11) return;  
  fHistEvCount->Fill(kEvINT11);
  
  // 1. Validate trigger requesting a SPD fastor (make sure the event was triggered also because of SPD)
  //    If we don't do this, we may be putting a bias on the trigger efficiency (an event triggered by V0OR
  //    is more likely to fire also V0AND than a generic INEL event
  if (fastorOnl1  < 1) return;
  fHistEvCount->Fill(kEvFOOnline);
  // 2. Request a reconstructed tracklet in |eta| < 1 (INEL>0 definition)
  //    For this, we have to loop over all tracklets to check eta
  UInt_t ntracklet = mult->GetNumberOfTracklets();		// tracklets
  Bool_t isOfflineTrig = kFALSE;
  Int_t ntrkEtaLess1 = 0;
  for(UInt_t itracklet = 0; itracklet < ntracklet; itracklet++){
    if (mult->GetEta(itracklet) < 1) {
      ntrkEtaLess1++;
    }
  }
  if (ntrkEtaLess1 >= fNtrkCut) { // needs to be out of the loop, becasue we use the total number of tracklets below.
    isOfflineTrig = kTRUE;
  }    
  if(!isOfflineTrig) return; //I return here because the assumption is that the efficiency of V0AND does not depend on the fact that a tracklet is reconstructed.
  fHistEvCount->Fill(kEvInelGT0);
  
  // ================== INT7 SELECTION =====================

  // 3. Check if int7 is also present (do this twice, both using isINT7 and filling fraction of triggers)
  //    The online one is already defined above as "isINT7"
  //    The offline is defined here:
  Bool_t isV0ANDOffline = fTriggerAnalysis -> IsOfflineTriggerFired(event, AliTriggerAnalysis::kV0AND);

  // ================= FILL HISTOGRAMS =====================
  // Selection Types
  // kAll, kINT11, kINT11EvSel, kINT11EVSelandINT7Online,
  fHistV0M                   [kHistoINT11]      ->Fill(v0mCentr);
  if(isINT7)         fHistV0M[kHistoINT7]       ->Fill(v0mCentr);
  if(isV0ANDOffline) fHistV0M[kHistoINT7Offline]->Fill(v0mCentr);

  
  if(isEvSelCentralityOk) {
    fHistV0M                   [kHistoINT11EvSel]       ->Fill(v0mCentr);
    if(isINT7)         fHistV0M[kHistoINT7EvSel]       ->Fill(v0mCentr);
    if(isV0ANDOffline) fHistV0M[kHistoINT7OfflineEvSel]->Fill(v0mCentr);    
  }

  // Ntracklets
  if (v0mCentr > 70 && v0mCentr < 100) {
    fHistNtrk70100                   [kHistoINT11]      ->Fill(ntrkEtaLess1);
    if(isINT7)         fHistNtrk70100[kHistoINT7]       ->Fill(ntrkEtaLess1);
    if(isV0ANDOffline) fHistNtrk70100[kHistoINT7Offline]->Fill(ntrkEtaLess1);

  
    if(isEvSelCentralityOk) {
      fHistNtrk70100                   [kHistoINT11EvSel]       ->Fill(ntrkEtaLess1);
      if(isINT7)         fHistNtrk70100[kHistoINT7EvSel]       ->Fill(ntrkEtaLess1);
      if(isV0ANDOffline) fHistNtrk70100[kHistoINT7OfflineEvSel]->Fill(ntrkEtaLess1);    
    }
    
  }
  
  



  // NEW HISTO should be filled before this point, as PostData puts the
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskTrigEff::Terminate(Option_t *) 
{
  // Draw result to screen, or perform fitting, normalizations
  // Called once at the end of the query
  return;
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
                
  // Get the physics selection histograms with the selection statistics
  //AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  //AliESDInputHandler *inputH = dynamic_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
  //TH2F *histStat = (TH2F*)inputH->GetStatistics();
   
   
  // NEW HISTO should be retrieved from the TList container in the above way,
  // so it is available to draw on a canvas such as below

}
