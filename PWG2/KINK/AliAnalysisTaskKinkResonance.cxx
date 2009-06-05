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
//        products. It provides basic plots as well as plots helping to calculate the corrections.
//        Usage: To analyse a resonance having a kaon in its decay products, one has to modify the integer 
//        variables resonancePDG, daughter1 and daughter2 accordingly as well as daughter1Mass  and daughter2Mass.
//        Also, depending on the analysis mode (ESD or MC), fAnalysisType in the constructor must also be changed 
//-----------------------------------------------------------------------------------------------------------------
#include "TCanvas.h"
#include "TChain.h"
#include "TTree.h"
#include "TH2D.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliResonanceKink.h"
#include "AliAnalysisTaskKinkResonance.h"

ClassImp(AliAnalysisTaskKinkResonance)

//_____________________________________________________________________________
AliAnalysisTaskKinkResonance::AliAnalysisTaskKinkResonance() 
  : AliAnalysisTaskSE(), fESD(0), fmcEventH(0), fList(0), fKinkResonance(0)

{
  // Constructor
}

//______________________________________________________________________________
AliAnalysisTaskKinkResonance::AliAnalysisTaskKinkResonance(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fmcEventH(0), fList(0), fKinkResonance(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskKinkResonance::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
      
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    
    if (!eventHandler) {
       Printf("ERROR: Could not retrieve MC event handler");
       return;
     }

     fmcEventH = eventHandler->MCEvent(); 
  }
}

//________________________________________________________________________
void AliAnalysisTaskKinkResonance::CreateOutputObjects() 
{
  // Create histograms
  // Called once
  
  fList=new TList();
  fList= fKinkResonance->GetHistogramList();
   
}

//________________________________________________________________________
void AliAnalysisTaskKinkResonance::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
   
   fKinkResonance->Analyse(fESD, fmcEventH);
   
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
  //TH1D* h=(TH1D *) fList->At(15);
  //TCanvas *c1 = new TCanvas("AliAnalysisTaskKinkResonance","Pt MC",10,10,510,510);
  //c1->cd(1);
  //if (h) h->DrawCopy("E");

}

