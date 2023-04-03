/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Archita Rani Dash                                              *
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

#include "TChain.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

//*****************************//
#include "AliAnalysisJetValidation.h"
//*****************************//
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODTrack.h"
#include <math.h>
#include <AliPIDResponse.h>
#include "AliMultSelection.h"
#include "AliEmcalJet.h"
#include "AliEmcalList.h"
#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliParticleContainer.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
//#include "AliAnalysisTaskParticleInJet.h"
#include "AliAODMCParticle.h"
#include "AliEmcalTrackSelection.h"
#include "AliVTrack.h"

#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include "AliAnalysisHelperJetTasks.h"
#include "AliGenPythiaEventHeader.h"
#include "TChain.h"
#include <map>

class AliAnalysisTaskEmcalJetValidation;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskEmcalJetValidation) // classimp: necessary for root

AliAnalysisTaskEmcalJetValidation::AliAnalysisTaskEmcalJetValidation() : AliAnalysisTaskEmcalJet(),
     fAOD(0), fOutputList(0), fHistJetPt(0),fHistJetPhi(0),fHistJetEta(0), jetconrec(0), jetcongen(0), fHistNEvents(0)
{

}

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetValidation::AliAnalysisTaskEmcalJetValidation(const char* name) : AliAnalysisTaskEmcalJet(name, kTRUE),
    fAOD(0), fOutputList(0), fHistJetPt(0),fHistJetPhi(0),fHistJetEta(0), fPIDResponse(0), jetconrec(0), jetcongen(0),fHistNEvents(0)
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskEmcalJetValidation::~AliAnalysisTaskEmcalJetValidation()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}

//_________________________________________________
void AliAnalysisTaskEmcalJetValidation::UserCreateOutputObjects()
{
	Printf("Check done %i",__LINE__);

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	if (man)
	{
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
	 }

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fHistPt = new TH1F("fHistPt", "fHistPt", 100, 0, 10);
  fHistNEvents = new TH1F("hNEvents", "Number of processed events", 15, -0.5, 14.5);
  //fHistjet = new TH1F("fHistjet", "fHistjet", 100, 0, 10);

  fHistJetPt = new TH1F("jetPt", "inclusive jetPt ; p_{T} (GeV/#it{c})", nBinsPt, 0, 100);
  fHistJetPhi = new TH1F("jetPhi", "inclusive jet #phi ; #phi ", nBinsPhi, 0, 6.4);
  fHistJetEta = new TH1F("jetEta", "inclusive jet #eta ; #eta ", nBinsEta, -0.9, 0.9);

  fOutputList->Add(fHistPt);
  fOutputList->Add(fHistNEvents);
  fOutputList->Add(fHistJetPt);
  fOutputList->Add(fHistJetPhi);
  fOutputList->Add(fHistJetEta);

  //fOutputList->Add(fHistNEvents);

  PostData(1, fOutputList);

}

void AliAnalysisTaskEmcalJetValidation::UserExec(Option_t *)
{
	Printf("Check done %i",__LINE__);
/*  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    //ESD??
  if(!fAOD) return;
*/
  AliVEvent* vevt = InputEvent();
  AliESDEvent* esd = (AliESDEvent*)(vevt);
  if (!esd) {
    printf("AliAnalysisTaskEmcalJetValidation::UserExec(): bad ESD\n");
    return;
}

  fHistNEvents->Fill(1);

  jetcongen = static_cast<AliJetContainer*>(fJetCollArray.At(1));   //jet container for mc generated jets = particle level
  jetcongen->ResetCurrentID();

  jetconrec = static_cast<AliJetContainer*>(fJetCollArray.At(0));   //jet container for reconstructed jets =detector level
  jetconrec->ResetCurrentID();

  //_____________________________________________
  //JET CONTAINERS
  AliEmcalJet * jetrec = 0;
  AliEmcalJet * jetmatched = 0;

  //_________________________________________________________
  //      LOOP OVER JETS  DETECTOR LEVEL JETS
  //_________________________________________________________
  while ((jetrec = jetconrec->GetNextAcceptJet())) 		//start of jet loop:
  {
   if(!jetrec) continue;
   if(jetrec->GetNumberOfTracks()==0) continue;

   fJetRecPt = jetrec->Pt();	//jet pt-distribution


   printf("Generated: fJetRecPt=%f\n",fJetRecPt);

   fHistJetPt->Fill(jetrec->Pt());
  }

  PostData(1, fOutputList);

}
//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetValidation::Terminate(Option_t *)
{
  PostData(1, fOutput);

  fOutput = dynamic_cast<AliEmcalList*> (GetOutputData(1)); // '1' refers to the output slot
  if(!fOutput) {
    printf("ERROR: Output list not available\n");
    return;
   }
}
