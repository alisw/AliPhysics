//#ifndef ALIANALYSISTASKMUONAODCREATION_CXX
//#define ALIANALYSISTASKMUONAODCREATION_CXX

/* $Id$ */

#include <TChain.h>
#include <TTree.h>
#include <TList.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TH1.h>

#include "AliAnalysisTaskMuonAODCreation.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliVEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskMuonAODCreation)

//__________________________________________________________________________
AliAnalysisTaskMuonAODCreation::AliAnalysisTaskMuonAODCreation() :
  fOutput(0x0),
  fTree(0x0),
  fOutputAOD(0x0)
{
}
//___________________________________________________________________________
AliAnalysisTaskMuonAODCreation::AliAnalysisTaskMuonAODCreation(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fOutput(0x0),
  fTree(0x0),
  fOutputAOD(0x0)
{
  // Constructor. Initialization of Inputs and Outputs
  //
    
  DefineOutput(1,TList::Class());

}

//___________________________________________________________________________
AliAnalysisTaskMuonAODCreation& AliAnalysisTaskMuonAODCreation::operator=(const AliAnalysisTaskMuonAODCreation& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskMuonAODCreation::AliAnalysisTaskMuonAODCreation(const AliAnalysisTaskMuonAODCreation& c) :
  AliAnalysisTaskSE(c),
  fOutput(c.fOutput),
  fTree(c.fTree),
  fOutputAOD(c.fOutputAOD)
 {
  //
  // Copy Constructor										
  //
}

//___________________________________________________________________________
AliAnalysisTaskMuonAODCreation::~AliAnalysisTaskMuonAODCreation() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskMuonAODCreation","Calling Destructor");
}

//___________________________________________________________________________
void AliAnalysisTaskMuonAODCreation::UserCreateOutputObjects(){
	
  AliAODHandler* handler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
 
  fOutputAOD   = handler->GetAOD();
  fTree = handler->GetTree();

  fOutput = new TList();
  fOutput->SetOwner(); 
 
  TH1D *pt_alltracks    = new TH1D("pt_alltracks","pt_alltracks",10,0,20);	
  TH1D *pt_muontracks    = new TH1D("pt_muontracks","pt_muontracks",10,0,20);	
	
  fOutput->Add(pt_alltracks); 	
  fOutput->Add(pt_muontracks); 	
  fOutput->ls(); 
} 



//_________________________________________________
void AliAnalysisTaskMuonAODCreation::UserExec(Option_t *)
{
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! aod ) {
    AliError("Cannot get AOD event");
    return;
  }  
  
  Int_t nMuons=0;
 
  for (Int_t j = 0; j<aod->GetNumberOfTracks(); j++) { 
    AliAODTrack *track = aod->GetTrack(j);
    ((TH1D*)(fOutput->FindObject("pt_alltracks")))->Fill(track->Pt()); 
    if(track->IsMuonTrack()) {
      nMuons++;
      ((TH1D*)(fOutput->FindObject("pt_muontracks")))->Fill(track->Pt()); 
    }  
  }
         
  AliAODHandler* outputHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());	     
   
  outputHandler->SetFillAOD(kFALSE);
   
  if(nMuons>0) { 
    outputHandler->SetFillAOD(kTRUE);
   
    PostData(0,fTree);
    PostData(1,fOutput);
  } else {
    return;
  }
}


//________________________________________________________________________
void AliAnalysisTaskMuonAODCreation::Terminate(Option_t *) 
{
  TCanvas *c = new TCanvas("c","plots",20,20,600,600);
  c->Divide(2,2);
  
  TH1D *h_pt_all = dynamic_cast<TH1D*> (fOutput->FindObject("pt_alltracks"));  
  TH1D *h_pt_muons = dynamic_cast<TH1D*> (fOutput->FindObject("pt_muontracks"));  
  c->cd(1);
  h_pt_all->Draw();
  c->cd(2);
  h_pt_muons->Draw();
  c->cd(3);
  h_pt_all->Draw();
  h_pt_muons->Draw("same");
}

