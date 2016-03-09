/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                       *
* Author: Baldo Sahlmueller, Friederike Bock                     *
* Version 1.0                                 *
*                                       *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims    *
* about the suitability of this software for any purpose. It is      *
* provided "as is" without express or implied warranty.               *
**************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on electromagnetic cocktail output
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliAnalysisTaskGammaCocktailMC.h"
#include "AliVParticle.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskGammaCocktailMC)

//________________________________________________________________________
AliAnalysisTaskGammaCocktailMC::AliAnalysisTaskGammaCocktailMC(): AliAnalysisTaskSE(),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fHistNEvents(NULL),
  fHistPtYPi0(NULL),
  fHistPtYPiPl(NULL),
  fHistPtYPiMi(NULL),
  fHistPtYEta(NULL),
  fHistPtYEtaPrim(NULL),
  fHistPtYOmega(NULL),
  fIsMC(1)
{
  
}

//________________________________________________________________________
AliAnalysisTaskGammaCocktailMC::AliAnalysisTaskGammaCocktailMC(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fHistNEvents(NULL),
  fHistPtYPi0(NULL),
  fHistPtYPiPl(NULL),
  fHistPtYPiMi(NULL),
  fHistPtYEta(NULL),
  fHistPtYEtaPrim(NULL),
  fHistPtYOmega(NULL),
  fIsMC(1)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaCocktailMC::~AliAnalysisTaskGammaCocktailMC()
{

}

//________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::UserCreateOutputObjects(){
  
  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer          = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }
  
  fHistNEvents                = new TH1F("NEvents", "NEvents", 1, 0, 1);
  fHistNEvents->Sumw2();
  fOutputContainer->Add(fHistNEvents);
  
  fHistPtYPi0                 = new TH2F("Pt_Y_Pi0","Pt_Y_Pi0", 500,0, 50, 400, -2.0, 2.0);
  fHistPtYPi0->Sumw2();
  fOutputContainer->Add(fHistPtYPi0);

  fHistPtYPiPl                = new TH2F("Pt_Y_PiPl","Pt_Y_PiPl", 500,0, 50, 400, -2.0, 2.0);
  fHistPtYPiPl->Sumw2();
  fOutputContainer->Add(fHistPtYPiPl);

  fHistPtYPiMi                = new TH2F("Pt_Y_PiMi","Pt_Y_PiMi", 500,0, 50, 400, -2.0, 2.0);
  fHistPtYPiMi->Sumw2();
  fOutputContainer->Add(fHistPtYPiMi);
  
  fHistPtYEta                 = new TH2F("Pt_Y_Eta","Pt_Y_Eta", 500,0, 50, 400, -2.0, 2.0);
  fHistPtYEta->Sumw2();
  fOutputContainer->Add(fHistPtYEta);

  fHistPtYEtaPrim             = new TH2F("Pt_Y_EtaPrim","Pt_Y_EtaPrim", 500,0, 50, 400, -2.0, 2.0);
  fHistPtYEtaPrim->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrim);

  fHistPtYOmega               = new TH2F("Pt_Y_Omega","Pt_Y_Omega", 500,0, 50, 400, -2.0, 2.0);
  fHistPtYOmega->Sumw2();
  fOutputContainer->Add(fHistPtYOmega);
   
  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::UserExec(Option_t *)
{

  fInputEvent = InputEvent();
//   cout << "I found an Event" << endl;
  
  fMCEvent = MCEvent();
  if(fMCEvent == NULL) fIsMC = 0;
  
  if (fIsMC==0) return;
//   cout << "I found an MC header" << endl;
    
  fMCStack = fMCEvent->Stack();
  if(fMCStack == NULL) fIsMC = 0;
  if (fIsMC==0) return;
  
  fHistNEvents->Fill(0.5);
//   cout << "the stack is intact" << endl;
  ProcessMCParticles();

  PostData(1, fOutputContainer);
}


//________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::ProcessMCParticles(){

  // Loop over all primary MC particle  
  for(UInt_t i = 0; i < fMCStack->GetNtrack(); i++) {
    // fill primary histograms
    TParticle* particle         = NULL;
    particle                    = (TParticle *)fMCStack->Particle(i);
    if (!particle) continue;
//     Bool_t hasMother            = kFALSE;
// //     cout << i << "\t"<< particle->GetMother(0) << endl;
//     if (particle->GetMother(0)>-1) 
//       hasMother                 = kTRUE;
//     TParticle* motherParticle   = NULL;
//     if( hasMother ) 
//       motherParticle            = (TParticle *)fMCStack->Particle(particle->GetMother(0));
//     if (motherParticle) 
//       hasMother                 = kTRUE;
//     else 
//       hasMother                 = kFALSE;

//     if (!(abs(particle->GetPdgCode()) == 111 || abs(particle->GetPdgCode()) == 221 || abs(particle->GetPdgCode()) == 331 ||
//       abs(particle->GetPdgCode()) == 223 || abs(particle->GetPdgCode()) == 211 )  )
//       continue;
    
    if (!(fabs(particle->Energy()-particle->Pz())>0.)) continue;
    Double_t yPre = (particle->Energy()+particle->Pz())/(particle->Energy()-particle->Pz());
//     cout << i << "\t"<< particle->GetPdgCode() << "\t"<< particle->Pz() << "\t" << particle->Energy()<< "\t" << particle->Energy()-particle->Pz() << "\t"<< yPre << endl;
    if (yPre == 0.) continue;
    
    Double_t y = 0.5*TMath::Log(yPre); 
//     if (y > 0.800) continue;
    switch(particle->GetPdgCode()){
      case 111:
        fHistPtYPi0->Fill(particle->Pt(), particle->Y());
        break;
      case 221:
        fHistPtYEta->Fill(particle->Pt(), particle->Y());
        break;
      case 331:
        fHistPtYEtaPrim->Fill(particle->Pt(), particle->Y());
        break;
      case 223:
        fHistPtYOmega->Fill(particle->Pt(), particle->Y());
        break;
      case 211:
        fHistPtYPiPl->Fill(particle->Pt(), particle->Y());
        break;
      case -211:
        fHistPtYPiMi->Fill(particle->Pt(), particle->Y());
        break;
    }    
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::Terminate(const Option_t *)
{
  
  //fOutputContainer->Print(); // Will crash on GRID
}


//_________________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::SetLogBinningXTH1(TH1* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::SetLogBinningXTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}
