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

//-----------------------------------------------------------------------
// Example of task running on AliEn (CAF?)
// which provides standard way of calculating acceptance and efficiency
// between different steps of the procedure.
// The ouptut of the task is a AliCFContainer from which the efficiencies
// can be calculated
//-----------------------------------------------------------------------
// Author : R. Vernet, Consorzio Cometa - Catania (it)
//-----------------------------------------------------------------------


#ifndef ALICFV0TASK_CXX
#define ALICFV0TASK_CXX
#include <TROOT.h>
#include <TInterpreter.h>

#include "AliCFV0Task.h"
#include "TCanvas.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TH1I.h"
#include "TChain.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliCFManager.h"
#include "AliCFCutBase.h"
#include "AliCFContainer.h"
#include "TChain.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliESDv0.h"
#include "AliV0vertexer.h"
#include "AliCFPair.h"

//__________________________________________________________________________
AliCFV0Task::AliCFV0Task() :
  fRebuildV0s(0),
  fV0PDG(310),
  fChain(0x0),
  fESD(0x0),
  fCFManager(0x0),
  fHistEventsProcessed(0x0)
{
//Defual ctor
}
//___________________________________________________________________________
AliCFV0Task::AliCFV0Task(const Char_t* name) :
  AliAnalysisTask(name,"AliCFV0Task"),
  fRebuildV0s(0),
  fV0PDG(310),
  fChain(0x0),
  fESD(0x0),
  fCFManager(0x0),
  fHistEventsProcessed(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFV0Task","Calling Constructor");
  DefineInput (0,TChain::Class());
  DefineOutput(0,TH1I::Class());
  DefineOutput(1,AliCFContainer::Class());
//   DefineOutput(2,TList::Class());
}

//___________________________________________________________________________
AliCFV0Task& AliCFV0Task::operator=(const AliCFV0Task& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTask::operator=(c) ;
    fRebuildV0s = c.fRebuildV0s;
    fV0PDG      = c.fV0PDG;
    fChain      = c.fChain;
    fESD        = c.fESD;
    fCFManager  = c.fCFManager;
    fHistEventsProcessed = c.fHistEventsProcessed;
  }
  return *this;
}

//___________________________________________________________________________
AliCFV0Task::AliCFV0Task(const AliCFV0Task& c) :
  AliAnalysisTask(c),
  fRebuildV0s(c.fRebuildV0s),
  fV0PDG(c.fV0PDG),
  fChain(c.fChain),
  fESD(c.fESD),
  fCFManager(c.fCFManager),
  fHistEventsProcessed(c.fHistEventsProcessed)
{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliCFV0Task::~AliCFV0Task() {
  //
  //destructor
  //
  Info("~AliCFV0Task","Calling Destructor");
  if (fChain)               delete fChain ;
  if (fESD)                 delete fESD ;
  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
}

//___________________________________________________________________________

void AliCFV0Task::Init()
{

}
//_________________________________________________
void AliCFV0Task::Exec(Option_t *)
{
  //
  // Main loop function
  //
  Info("Exec","") ;
  // Get the mc truth
  AliMCEventHandler* mcTruth = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());

  if (!mcTruth) Error("Exec","NO MC INFO FOUND... EXITING");

  // transform possible old AliESD into AliESDEvent

  if (fESD->GetAliESDOld()) fESD->CopyFromOldESD(); //transition to new ESD format

  //pass the MC evt handler to the cuts that need it 

  fCFManager->SetEventInfo(mcTruth);

  // Get the MC event 
  AliMCEvent* mcEvent = mcTruth->MCEvent();

  // MC-event selection
  Double_t containerInput[2] ;
        
  //loop on the MC event
  Info("Exec","Looping on MC event");
  for (Int_t ipart=0; ipart<mcEvent->GetNumberOfTracks(); ipart++) { 
    AliMCParticle *mcPart  = mcEvent->GetTrack(ipart);

    //check the MC-level cuts
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) continue;
    containerInput[0] = mcPart->Pt();
    containerInput[1] = mcPart->Eta() ;
    //fill the container for Gen-level selection
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepGenerated);
    
    //check the Acceptance-level cuts
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mcPart)) continue;
    //fill the container for Acceptance-level selection
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructible);
  }    


  //Now go to rec level
  Info("Exec","Looping on ESD event");

  if (fRebuildV0s) RebuildV0s() ;
  printf("There are %d V0s in event\n",fESD->GetNumberOfV0s());
  for (Int_t iV0 = 0; iV0<fESD->GetNumberOfV0s(); iV0++) {

    AliESDv0* esdV0 = fESD->GetV0(iV0);

    //check if mother reconstructed V0 can be associated to a MC V0
    Int_t labMCV0 = IsMcV0(esdV0,fESD,mcEvent->Stack()) ;
    if (labMCV0 <0) continue;

    esdV0->ChangeMassHypothesis(fV0PDG); //important to do that before entering the cut check

    AliCFPair pair(esdV0,fESD);
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,&pair)) continue;

    //check if associated MC v0 passes the cuts
    AliMCParticle* mcV0 = mcEvent->GetTrack(labMCV0);
    if (!mcV0) continue;

    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcV0)) continue; 
    
    //fill the container
    Double_t mom[3];
    esdV0->GetPxPyPz(mom[0],mom[1],mom[2]);
    Double_t pt2 = mom[0]*mom[0]+mom[1]*mom[1] ;
    Double_t pt  = TMath::Sqrt(pt2);
    Double_t energy  = TMath::Sqrt(pt2 + mom[2]*mom[2] + TMath::Power(TDatabasePDG::Instance()->GetParticle(fV0PDG)->Mass(),2));

    containerInput[0] = pt ;
    containerInput[1] = GetRapidity(energy,mom[2]);
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed) ;   
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,&pair)) continue ;
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepSelected);
  }
  
  fHistEventsProcessed->Fill(0);
  PostData(0,fHistEventsProcessed) ;
  PostData(1,fCFManager->GetParticleContainer()) ;
  
//   TList * list = new TList();
//   fCFManager->AddQAHistosToList(list);
//   PostData(2,list) ;
}


//___________________________________________________________________________
void AliCFV0Task::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  Info("Terminate","");
  AliAnalysisTask::Terminate();


  Double_t max1 = fCFManager->GetParticleContainer()->ShowProjection(0,0)->GetMaximum();
  Double_t max2 = fCFManager->GetParticleContainer()->ShowProjection(1,0)->GetMaximum();

  fCFManager->GetParticleContainer()->ShowProjection(0,0)->GetYaxis()->SetRangeUser(0,max1*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(0,1)->GetYaxis()->SetRangeUser(0,max1*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(0,2)->GetYaxis()->SetRangeUser(0,max1*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(0,3)->GetYaxis()->SetRangeUser(0,max1*1.2);

  fCFManager->GetParticleContainer()->ShowProjection(1,0)->GetYaxis()->SetRangeUser(0,max2*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(1,1)->GetYaxis()->SetRangeUser(0,max2*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(1,2)->GetYaxis()->SetRangeUser(0,max2*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(1,3)->GetYaxis()->SetRangeUser(0,max2*1.2);

  fCFManager->GetParticleContainer()->ShowProjection(0,0)->SetMarkerStyle(23) ;
  fCFManager->GetParticleContainer()->ShowProjection(0,1)->SetMarkerStyle(24) ;
  fCFManager->GetParticleContainer()->ShowProjection(0,2)->SetMarkerStyle(25) ;
  fCFManager->GetParticleContainer()->ShowProjection(0,3)->SetMarkerStyle(26) ;

  fCFManager->GetParticleContainer()->ShowProjection(1,0)->SetMarkerStyle(23) ;
  fCFManager->GetParticleContainer()->ShowProjection(1,1)->SetMarkerStyle(24) ;
  fCFManager->GetParticleContainer()->ShowProjection(1,2)->SetMarkerStyle(25) ;
  fCFManager->GetParticleContainer()->ShowProjection(1,3)->SetMarkerStyle(26) ;

  TCanvas * c =new TCanvas("c","",1400,800);
  c->Divide(4,2);

  c->cd(1);
  fCFManager->GetParticleContainer()->ShowProjection(0,0)->Draw("p");
  c->cd(2);
  fCFManager->GetParticleContainer()->ShowProjection(0,1)->Draw("p");
  c->cd(3);
  fCFManager->GetParticleContainer()->ShowProjection(0,2)->Draw("p");
  c->cd(4);
  fCFManager->GetParticleContainer()->ShowProjection(0,3)->Draw("p");
  c->cd(5);
  fCFManager->GetParticleContainer()->ShowProjection(1,0)->Draw("p");
  c->cd(6);
  fCFManager->GetParticleContainer()->ShowProjection(1,1)->Draw("p");
  c->cd(7);
  fCFManager->GetParticleContainer()->ShowProjection(1,2)->Draw("p");
  c->cd(8);
  fCFManager->GetParticleContainer()->ShowProjection(1,3)->Draw("p");

  c->SaveAs("plots.eps");

  delete fHistEventsProcessed ;
}

//___________________________________________________________________________
void AliCFV0Task::ConnectInputData(Option_t *) {
  //
  // Initialize branches.
  //
  Info("ConnectInputData","ConnectInputData of task %s\n",GetName());

  fChain = (TChain*)GetInputData(0);
  fChain->SetBranchStatus("*FMD*",0);
  fChain->SetBranchStatus("*CaloClusters*",0);
  fESD = new AliESDEvent();
  fESD->ReadFromTree(fChain);
}

//___________________________________________________________________________
void AliCFV0Task::CreateOutputObjects() {
  //HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  Info("CreateOutputObjects","CreateOutputObjects of task %s\n", GetName());

  //slot #0
  fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;
}

//___________________________________________________________________________
Int_t AliCFV0Task::GetV0Label(UInt_t lab1, UInt_t lab2, AliStack* stack) const {
  //
  // returns the label of the V0, given the labels of the 2 daughter tracks
  // returns -1 if the V0 is fake
  //
  
  TParticle* part1 = stack->Particle(lab1) ;
  TParticle* part2 = stack->Particle(lab2) ;
  
  Int_t part1MotherLab=part1->GetFirstMother();
  Int_t part2MotherLab=part2->GetFirstMother();
  
  if (part1MotherLab==-1 || part2MotherLab==-1) return -1 ;
  if (part1MotherLab != part2MotherLab )        return -1 ;
  if (stack->Particle(part1MotherLab)->GetPdgCode() != fV0PDG ) return -1 ;

  switch (fV0PDG) {
  case kK0Short : 
    if ( (part1->GetPdgCode()==  -211 && part2->GetPdgCode()==  211)  || 
	 (part1->GetPdgCode()==   211 && part2->GetPdgCode()== -211)   )   return part1MotherLab ;
    break ;
  case kLambda0 :
    if ( (part1->GetPdgCode()==  -211 && part2->GetPdgCode()== 2212)  || 
	 (part1->GetPdgCode()==  2212 && part2->GetPdgCode()== -211)   )   return part1MotherLab ;
    break ;
  case kLambda0Bar :
    if ( (part1->GetPdgCode()==   211 && part2->GetPdgCode()== -2212) || 
	 (part1->GetPdgCode()== -2212 && part2->GetPdgCode()==   211)  )   return part1MotherLab ;
    break ;
  default :
    return -1;
    break ;
  }
  
  return -1 ;
}

//___________________________________________________________________________
Double_t AliCFV0Task::GetRapidity(Double_t energy, Double_t pz) {
  if (energy == pz || energy == -pz) {
    printf("GetRapidity : ERROR : rapidity for 4-vector with E = Pz -- infinite result");
    return 999;
  }
  if (energy < pz) {
    printf("GetRapidity : ERROR : rapidity for 4-vector with E = Pz -- infinite result");
    return 999;
  }
  Double_t y = 0.5 * log((energy + pz) / (energy - pz));
  return y;
} 

//___________________________________________________________________________
Int_t AliCFV0Task::IsMcV0(AliESDv0* v0, AliESDEvent* esd, AliStack* stack) const {
  Int_t nindex = v0->GetNindex();
  Int_t pindex = v0->GetPindex();
  AliESDtrack *nTrack = esd->GetTrack(nindex) ;
  AliESDtrack *pTrack = esd->GetTrack(pindex) ;
  
  if (!nTrack || !pTrack) return -1 ;

  Int_t nlab  = nTrack->GetLabel() ;
  Int_t plab  = pTrack->GetLabel() ;
  
  if (nlab <0 || plab <0) return -1 ;

  Int_t v0Label = GetV0Label((UInt_t)nlab,(UInt_t)plab,stack) ;
  return v0Label ;
}


//___________________________________________________________________________
void AliCFV0Task::RebuildV0s() {
  
  fESD->ResetV0s();

  //These are pp cuts : to change if Pb+Pb !
  Double_t cuts[]={33,  // max. allowed chi2
		   0.1,// min. allowed negative daughter's impact parameter 
		   0.1,// min. allowed positive daughter's impact parameter 
		   0.1,// max. allowed DCA between the daughter tracks
		   0.999,// max. allowed cosine of V0's pointing angle
		   0.9,  // min. radius of the fiducial volume
		   100.   // max. radius of the fiducial volume
  };
  //   AliV0vertexer* v0Vertexer = new AliV0vertexer(cuts);
  //   v0Vertexer->SetVertex(primVertexPos);
  AliV0vertexer* v0Vertexer = new AliV0vertexer();
  v0Vertexer->SetCuts(cuts) ;
  v0Vertexer->Tracks2V0vertices(fESD);
  delete v0Vertexer ;
}

#endif
