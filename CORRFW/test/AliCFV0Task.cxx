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

#include "AliCFV0Task.h"
#include "TCanvas.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TH1I.h"
#include "AliMCEvent.h"
#include "AliCFManager.h"
#include "AliCFContainer.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliESDv0.h"
#include "AliV0vertexer.h"
#include "AliCFPair.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "TChain.h"
#include "AliCFParticleGenCuts.h"
#include "AliAODv0.h"

//__________________________________________________________________________
AliCFV0Task::AliCFV0Task() :
  AliAnalysisTaskSE(),
  fRebuildV0s(0),
  fV0PDG(310),
  fCFManager(0x0),
  fHistEventsProcessed(0x0)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliCFV0Task::AliCFV0Task(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fRebuildV0s(0),
  fV0PDG(310),
  fCFManager(0x0),
  fHistEventsProcessed(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFV0Task","Calling Constructor");
  /*
    DefineInput(0) and DefineOutput(0)
    are taken care of by AliAnalysisTaskSE constructor
  */
  DefineOutput(1,TH1I::Class());
  DefineOutput(2,AliCFContainer::Class());
//   DefineOutput(3,TList::Class());
}

//___________________________________________________________________________
AliCFV0Task& AliCFV0Task::operator=(const AliCFV0Task& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
    fRebuildV0s = c.fRebuildV0s;
    fV0PDG      = c.fV0PDG;
    fCFManager  = c.fCFManager;
    fHistEventsProcessed = c.fHistEventsProcessed;
  }
  return *this;
}

//___________________________________________________________________________
AliCFV0Task::AliCFV0Task(const AliCFV0Task& c) :
  AliAnalysisTaskSE(c),
  fRebuildV0s(c.fRebuildV0s),
  fV0PDG(c.fV0PDG),
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
  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
}

//_________________________________________________
void AliCFV0Task::UserExec(Option_t *)
{
  //
  // Main loop function
  //
  Info("UserExec","") ;

  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }

  if (!fMCEvent) Error("UserExec","NO MC INFO FOUND!");
  fCFManager->SetEventInfo(fMCEvent);

  Bool_t isESDEvent = strcmp(fInputEvent->ClassName(),"AliESDEvent") == 0 ? kTRUE : kFALSE ;
  Bool_t isAODEvent = strcmp(fInputEvent->ClassName(),"AliAODEvent") == 0 ? kTRUE : kFALSE ;

  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(fInputEvent);
  AliAODEvent* fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
  

  // MC-event selection
  Double_t containerInput[2] ;
        
  //loop on the MC event
  Info("UserExec","Looping on MC event");
  for (Int_t ipart=0; ipart<fMCEvent->GetNumberOfTracks(); ipart++) { 
    AliMCParticle *mcPart  = fMCEvent->GetTrack(ipart);

    //check the MC-level cuts
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) continue;
    containerInput[0] = mcPart->Pt();
    containerInput[1] = mcPart->Y() ;
    //fill the container for Gen-level selection
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepGenerated);
    
    //check the Acceptance-level cuts
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mcPart)) continue;
    //fill the container for Acceptance-level selection
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructible);
  }    


  //Now go to rec level
  Info("UserExec","Looping on %s",fInputEvent->ClassName());
 
  //SET THE ESD AS EVENT INFO IN RECONSTRUCTION CUTS
  TObjArray* fCutsReco = fCFManager->GetParticleCutsList(AliCFManager::kPartRecCuts);
  TObjArrayIter iter1(fCutsReco);
  AliCFCutBase *cut = 0;
  while ( (cut = (AliCFCutBase*)iter1.Next()) ) {
    cut->SetEvtInfo(fInputEvent);
  }

  if (isESDEvent && fRebuildV0s) RebuildV0s(fESD) ;

  Info("UserExec","There are %d V0s in event",fInputEvent->GetNumberOfV0s());
  
  AliESDv0*  esdV0 = 0x0;
  AliAODv0*  aodV0 = 0x0;
  AliCFPair* pair  = 0x0;
  Int_t      labMCV0 = 0;

  for (Int_t iV0 = 0; iV0<fInputEvent->GetNumberOfV0s(); iV0++) {

    if (isESDEvent) {
      esdV0 = fESD->GetV0(iV0);
      //check if mother reconstructed V0 can be associated to a MC V0
      labMCV0 = IsMcV0(esdV0);
      if (labMCV0 <0) continue;
      pair = new AliCFPair(esdV0,fESD);
    }
    else if (isAODEvent) {
      aodV0 = fAOD->GetV0(iV0);
      labMCV0 = IsMcV0(aodV0);
      if (labMCV0 <0) continue;
      pair = new AliCFPair(aodV0);
    }
    else {
      Error("UserExec","Error: input data file is not an ESD nor an AOD");
      return;
    }
    
    pair->SetV0PDG(fV0PDG);
    pair->SetLabel(labMCV0);
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,pair)) continue;
    
    //check if associated MC v0 passes the cuts
    AliMCParticle* mcV0 = fMCEvent->GetTrack(labMCV0);
    if (!mcV0) continue;

    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcV0)) continue; 
    
    //fill the container
    Double_t pt, rapidity;
    
    if (isESDEvent) {
      Double_t mom[3];
      esdV0->GetPxPyPz(mom[0],mom[1],mom[2]);
      Double_t pt2 = mom[0]*mom[0]+mom[1]*mom[1] ;
      pt  = TMath::Sqrt(pt2);
      Double_t energy  = TMath::Sqrt(pt2 + mom[2]*mom[2] + TMath::Power(TDatabasePDG::Instance()->GetParticle(fV0PDG)->Mass(),2));
      rapidity = GetRapidity(energy,mom[2]) ;
    }
    else {
      pt = aodV0->Pt();
      rapidity = aodV0->Y(fV0PDG);
    }
    
    containerInput[0] = pt ;
    containerInput[1] = rapidity ;
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed) ;   
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,pair)) continue ;
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepSelected);
    
    delete pair;
  }
  
  fHistEventsProcessed->Fill(0);
  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
}


//___________________________________________________________________________
void AliCFV0Task::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  Info("Terminate","");
  AliAnalysisTaskSE::Terminate();


  //draw some example plots....

  AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));

  TH1D* h00 =   cont->ShowProjection(0,0) ;
  TH1D* h01 =   cont->ShowProjection(0,1) ;
  TH1D* h02 =   cont->ShowProjection(0,2) ;
  TH1D* h03 =   cont->ShowProjection(0,3) ;

  TH1D* h10 =   cont->ShowProjection(1,0) ;
  TH1D* h11 =   cont->ShowProjection(1,1) ;
  TH1D* h12 =   cont->ShowProjection(1,2) ;
  TH1D* h13 =   cont->ShowProjection(1,3) ;

  Double_t max1 = h00->GetMaximum();
  Double_t max2 = h10->GetMaximum();

  h00->GetYaxis()->SetRangeUser(0,max1*1.2);
  h01->GetYaxis()->SetRangeUser(0,max1*1.2);
  h02->GetYaxis()->SetRangeUser(0,max1*1.2);
  h03->GetYaxis()->SetRangeUser(0,max1*1.2);

  h10->GetYaxis()->SetRangeUser(0,max2*1.2);
  h11->GetYaxis()->SetRangeUser(0,max2*1.2);
  h12->GetYaxis()->SetRangeUser(0,max2*1.2);
  h13->GetYaxis()->SetRangeUser(0,max2*1.2);

  h00->SetMarkerStyle(23) ;
  h01->SetMarkerStyle(24) ;
  h02->SetMarkerStyle(25) ;
  h03->SetMarkerStyle(26) ;

  h10->SetMarkerStyle(23) ;
  h11->SetMarkerStyle(24) ;
  h12->SetMarkerStyle(25) ;
  h13->SetMarkerStyle(26) ;

  TCanvas * c =new TCanvas("c","",1400,800);
  c->Divide(4,2);

  c->cd(1);
  h00->Draw("p");
  c->cd(2);
  h01->Draw("p");
  c->cd(3);
  h02->Draw("p");
  c->cd(4);
  h03->Draw("p");
  c->cd(5);
  h10->Draw("p");
  c->cd(6);
  h11->Draw("p");
  c->cd(7);
  h12->Draw("p");
  c->cd(8);
  h13->Draw("p");

  c->SaveAs("plots.eps");
}

//___________________________________________________________________________
void AliCFV0Task::UserCreateOutputObjects() {
  //HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());

  //slot #1
  OpenFile(1);
  fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;
}

//___________________________________________________________________________
Int_t AliCFV0Task::GetV0Label(UInt_t lab1, UInt_t lab2) const {
  //
  // returns the label of the V0, given the labels of the 2 daughter tracks
  // returns -1 if the V0 is fake
  //
  
  AliStack* stack = fMCEvent->Stack();
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
  //
  // calculates rapidity, checking energy is larger than pz, otherwise returns 999.
  //
  
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
Int_t AliCFV0Task::IsMcV0(AliESDv0* v0) const {
  //
  // check if the passed V0 is associated to a MC one, 
  //   and returns the corresponding geant label.
  // returns -1 if the V0 is fake (i.e. label<0).
  //

  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(fInputEvent);

  Int_t nindex = v0->GetNindex();
  Int_t pindex = v0->GetPindex();
  AliESDtrack *nTrack = fESD->GetTrack(nindex) ;
  AliESDtrack *pTrack = fESD->GetTrack(pindex) ;
  
  if (!nTrack || !pTrack) return -1 ;

  Int_t nlab  = nTrack->GetLabel() ;
  Int_t plab  = pTrack->GetLabel() ;
  
  if (nlab <0 || plab <0) return -1 ;

  return GetV0Label((UInt_t)nlab,(UInt_t)plab) ;
}

//___________________________________________________________________________
Int_t AliCFV0Task::IsMcV0(AliAODv0* v0) const {
  //
  // check if the passed V0 is associated to a MC one, 
  //   and returns the corresponding geant label.
  // returns -1 if the V0 is fake (i.e. label<0).
  //

  AliAODVertex * vtx = v0->GetSecondaryVtx();

  AliAODTrack *nTrack = (AliAODTrack*)vtx->GetDaughter(1); //neg is filled after pos in AliAnalysisTaskStrange
  AliAODTrack *pTrack = (AliAODTrack*)vtx->GetDaughter(0);
  
  if (!nTrack || !pTrack) return -1 ;

  Int_t nlab  = nTrack->GetLabel() ;
  Int_t plab  = pTrack->GetLabel() ;

  if (nlab <0 || plab <0) return -1 ;

  return GetV0Label((UInt_t)nlab,(UInt_t)plab) ;
}


//___________________________________________________________________________
void AliCFV0Task::RebuildV0s(AliESDEvent* fESD) {
  
  fESD->ResetV0s();

  //These are pp cuts : to change if Pb+Pb !
  Double_t cuts[]={33,  // max. allowed chi2
		   0.01,// min. allowed negative daughter's impact parameter 
		   0.01,// min. allowed positive daughter's impact parameter 
		   0.5,// max. allowed DCA between the daughter tracks
		   0.98,// max. allowed cosine of V0's pointing angle
		   0.2,  // min. radius of the fiducial volume
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
