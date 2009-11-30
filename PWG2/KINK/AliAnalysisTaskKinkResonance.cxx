/**************************************************************************
 * Author: Paraskevi Ganoti, University of Athens (pganoti@phys.uoa.gr)   *
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
 
//----------------------------------------------------------------------------------------------------------------
//                        class AliAnalysisTaskKinkResonance
//        Example of an analysis task for reconstructing resonances having at least one kaon-kink in their decay 
//        products. 
//-----------------------------------------------------------------------------------------------------------------
#include "TCanvas.h"
#include "TH1.h"

#include "AliESDEvent.h"
#include "AliAnalysisTaskKinkResonance.h"
#include "AliResonanceKink.h"

ClassImp(AliAnalysisTaskKinkResonance)

//______________________________________________________________________________
AliAnalysisTaskKinkResonance::AliAnalysisTaskKinkResonance(const char *dname) 
  : AliAnalysisTaskSE(dname), fList(0), fKinkResonance(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskKinkResonance::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  
  fList= fKinkResonance->GetHistogramList();
   
}

//________________________________________________________________________
void AliAnalysisTaskKinkResonance::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
   
  AliVEvent *event = InputEvent();
  if (!event) {
     Printf("ERROR: Could not retrieve event");
     return;
  }

  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  
  if (!esd) {
     Printf("ERROR: Could not retrieve esd");
     return;
  }
  
  if((fKinkResonance->GetAnalysisType() == "MC")||(fKinkResonance->GetAnalysisType() == "ESD")) {
    AliMCEvent* mcEvent = MCEvent();
    if (!mcEvent) {
       Printf("ERROR: Could not retrieve MC event");
       return;
    }
    fKinkResonance->SetDebugLevel(fDebug);
    fKinkResonance->Analyse(esd, mcEvent);
  }
  
  if(fKinkResonance->GetAnalysisType() == "DATA") {
    fKinkResonance->SetDebugLevel(fDebug);
    fKinkResonance->Analyse(esd);
  }
   
  PostData(1, fList);
   
}      

//________________________________________________________________________
void AliAnalysisTaskKinkResonance::Terminate(Option_t *) 
{
  // Draw result to the screen 
  // Called once at the end of the query

  fList = dynamic_cast<TList*> (GetOutputData(1));

  if (!fList) {
    Printf("ERROR: fList not available");
    return;
  }
  
  TH1D* h=(TH1D *) fList->At(1);
  TCanvas *c1 = new TCanvas("AliAnalysisTaskKinkResonance","Inv mass",10,10,510,510);
  c1->cd(1);
  if (h) h->DrawCopy("E");
  
}

