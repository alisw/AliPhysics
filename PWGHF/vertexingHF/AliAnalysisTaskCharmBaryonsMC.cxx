
//----------------------------------------------------------------
//      Implementation of Class AliAnalysisTaskCharmBaryonsMC
//
// Task used to analize simulations at generation level (i.e. only
// needs galice.root and Kinematics.root).
// 
// This task make pt and rapidity histograms of Lc/D0/Xic0
//
// Author: Yosuke Watanabe
//
//----------------------------------------------------------------


#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"

#include "AliAnalysisTaskCharmBaryonsMC.h"
#include "TGraphErrors.h"
#include "AliLog.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"

#include <iostream>
#include "AliPDG.h"
#include "AliGenDPMjetEventHeader.h"
#include "TFile.h"

using namespace std;

ClassImp(AliAnalysisTaskCharmBaryonsMC)

  AliAnalysisTaskCharmBaryonsMC::AliAnalysisTaskCharmBaryonsMC() 
: AliAnalysisTaskSE(), fMyOut(0), 
  fHistEvt(0),
  fHistMult(0),
  fHistPtPromptD0(0),
  fHistPtPromptLc(0),
  fHistPtPromptXic0(0),
  fHistPtFeeddownD0(0),
  fHistPtFeeddownLc(0),
  fHistPtFeeddownXic0(0),
  fHistPtInclusiveD0(0),
  fHistPtInclusiveLc(0),
  fHistPtInclusiveXic0(0),
  fHistPtvsRapidityPromptD0(0),
  fHistPtvsRapidityPromptLc(0),
  fHistPtvsRapidityPromptXic0(0),
  fHistPtvsRapidityFeeddownD0(0),
  fHistPtvsRapidityFeeddownLc(0),
  fHistPtvsRapidityFeeddownXic0(0),
  fHistPtvsRapidityInclusiveD0(0),
  fHistPtvsRapidityInclusiveLc(0),
  fHistPtvsRapidityInclusiveXic0(0)
{

  //
  // default constructor
  //
}

//________________________________________________________________________
  AliAnalysisTaskCharmBaryonsMC::AliAnalysisTaskCharmBaryonsMC(const char *name) 
: AliAnalysisTaskSE(name), fMyOut(0),
  fHistEvt(0),
  fHistMult(0),
  fHistPtPromptD0(0),
  fHistPtPromptLc(0),
  fHistPtPromptXic0(0),
  fHistPtFeeddownD0(0),
  fHistPtFeeddownLc(0),
  fHistPtFeeddownXic0(0),
  fHistPtInclusiveD0(0),
  fHistPtInclusiveLc(0),
  fHistPtInclusiveXic0(0),
  fHistPtvsRapidityPromptD0(0),
  fHistPtvsRapidityPromptLc(0),
  fHistPtvsRapidityPromptXic0(0),
  fHistPtvsRapidityFeeddownD0(0),
  fHistPtvsRapidityFeeddownLc(0),
  fHistPtvsRapidityFeeddownXic0(0),
  fHistPtvsRapidityInclusiveD0(0),
  fHistPtvsRapidityInclusiveLc(0),
  fHistPtvsRapidityInclusiveXic0(0)
{
  //
  // constructor
  //


  AliPDG::AddParticlesToPdgDataBase();

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskCharmBaryonsMC::~AliAnalysisTaskCharmBaryonsMC() {

  // destructor

  if(fMyOut) {
    // fMyOut owns the histos
    delete fMyOut;
    fMyOut = 0;
  }

}


//________________________________________________________________________
void AliAnalysisTaskCharmBaryonsMC::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fMyOut = new TList();
  BookHistograms();

  fMyOut->SetOwner();

  // Suppress annoying printout
  AliLog::SetGlobalLogLevel(AliLog::kError);

}

void AliAnalysisTaskCharmBaryonsMC::BookHistograms() {
  //
  //
  //
  fHistEvt = new TH1F("fHistEvt","Number of events", 1, -0.5, 0.5); 
  fMyOut->Add(fHistEvt);
  fHistMult = new TH1F("fHistMult","Multiplicity", 100, 0.0, 1000.); 
  fMyOut->Add(fHistMult);

  fHistPtPromptD0 = new TH1D("fHistPtPromptD0","Prompt D0", 100, 0., 20.); 
  fMyOut->Add(fHistPtPromptD0);
  fHistPtPromptLc = new TH1D("fHistPtPromptLc","Prompt Lc", 100, 0., 20.); 
  fMyOut->Add(fHistPtPromptLc);
  fHistPtPromptXic0 = new TH1D("fHistPtPromptXic0","Prompt Xic0", 100, 0., 20.); 
  fMyOut->Add(fHistPtPromptXic0);
  fHistPtFeeddownD0 = new TH1D("fHistPtFeeddownD0","Feeddown D0", 100, 0., 20.); 
  fMyOut->Add(fHistPtFeeddownD0);
  fHistPtFeeddownLc = new TH1D("fHistPtFeeddownLc","Feeddown Lc", 100, 0., 20.); 
  fMyOut->Add(fHistPtFeeddownLc);
  fHistPtFeeddownXic0 = new TH1D("fHistPtFeeddownXic0","Feeddown Xic0", 100, 0., 20.); 
  fMyOut->Add(fHistPtFeeddownXic0);
  fHistPtInclusiveD0 = new TH1D("fHistPtInclusiveD0","Inclusive D0", 100, 0., 20.); 
  fMyOut->Add(fHistPtInclusiveD0);
  fHistPtInclusiveLc = new TH1D("fHistPtInclusiveLc","Inclusive Lc", 100, 0., 20.); 
  fMyOut->Add(fHistPtInclusiveLc);
  fHistPtInclusiveXic0 = new TH1D("fHistPtInclusiveXic0","Inclusive Xic0", 100, 0., 20.); 
  fMyOut->Add(fHistPtInclusiveXic0);

  fHistPtvsRapidityPromptD0 = new TH2D("fHistPtvsRapidityPromptD0","Prompt D0", 100, 0., 20.,20,-5.,5.); 
  fMyOut->Add(fHistPtvsRapidityPromptD0);
  fHistPtvsRapidityPromptLc = new TH2D("fHistPtvsRapidityPromptLc","Prompt Lc", 100, 0., 20.,20,-5.,5.); 
  fMyOut->Add(fHistPtvsRapidityPromptLc);
  fHistPtvsRapidityPromptXic0 = new TH2D("fHistPtvsRapidityPromptXic0","Prompt Xic0", 100, 0., 20.,20,-5.,5.); 
  fMyOut->Add(fHistPtvsRapidityPromptXic0);
  fHistPtvsRapidityFeeddownD0 = new TH2D("fHistPtvsRapidityFeeddownD0","Feeddown D0", 100, 0., 20.,20,-5.,5.); 
  fMyOut->Add(fHistPtvsRapidityFeeddownD0);
  fHistPtvsRapidityFeeddownLc = new TH2D("fHistPtvsRapidityFeeddownLc","Feeddown Lc", 100, 0., 20.,20,-5.,5.); 
  fMyOut->Add(fHistPtvsRapidityFeeddownLc);
  fHistPtvsRapidityFeeddownXic0 = new TH2D("fHistPtvsRapidityFeeddownXic0","Feeddown Xic0", 100, 0., 20.,20,-5.,5.); 
  fMyOut->Add(fHistPtvsRapidityFeeddownXic0);
  fHistPtvsRapidityInclusiveD0 = new TH2D("fHistPtvsRapidityInclusiveD0","Inclusive D0", 100, 0., 20.,20,-5.,5.); 
  fMyOut->Add(fHistPtvsRapidityInclusiveD0);
  fHistPtvsRapidityInclusiveLc = new TH2D("fHistPtvsRapidityInclusiveLc","Inclusive Lc", 100, 0., 20.,20,-5.,5.); 
  fMyOut->Add(fHistPtvsRapidityInclusiveLc);
  fHistPtvsRapidityInclusiveXic0 = new TH2D("fHistPtvsRapidityInclusiveXic0","Inclusive Xic0", 100, 0., 20.,20,-5.,5.); 
  fMyOut->Add(fHistPtvsRapidityInclusiveXic0);

}

//________________________________________________________________________
void AliAnalysisTaskCharmBaryonsMC::UserExec(Option_t *) 
{
  //
  // Main loop
  // Called for each event
  //

  // also a AliEvent...
  //  AliVEvent* mcEvent = MCEvent();
  AliMCEvent* mcEvent = MCEvent();
  //  AliGenEventHeader * htmp = mcEvent->GenEventHeader();

  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }

  fHistEvt->Fill(0);
  fHistMult->Fill(mcEvent->GetNumberOfTracks());

  // Track loop
  for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
    AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    //Bool_t isPrimary = mcEvent->Stack()->IsPhysicalPrimary(iTrack);
    Int_t pdgcode = abs(track->PdgCode());

    if(pdgcode==421 || pdgcode==4122 || pdgcode==4132){
      Int_t pdgarray_chad[100];
      Int_t labelarray_chad[100];
      Int_t ngen_chad=0;
      GetHistory(track,mcEvent,pdgarray_chad,labelarray_chad,ngen_chad);
      //cout<<isPrimary<<" "<<pdgcode<<" "<<pdgarray_chad[0]<<" "<<pdgarray_chad[1]<<endl;

      Bool_t isbfd = FromBottom(pdgarray_chad);

      if(pdgcode==421) fHistPtvsRapidityInclusiveD0->Fill(track->Pt(),track->Y());
      if(pdgcode==4122) fHistPtvsRapidityInclusiveLc->Fill(track->Pt(),track->Y());
      if(pdgcode==4132) fHistPtvsRapidityInclusiveXic0->Fill(track->Pt(),track->Y());
      if(isbfd){
        if(pdgcode==421) fHistPtvsRapidityFeeddownD0->Fill(track->Pt(),track->Y());
        if(pdgcode==4122) fHistPtvsRapidityFeeddownLc->Fill(track->Pt(),track->Y());
        if(pdgcode==4132) fHistPtvsRapidityFeeddownXic0->Fill(track->Pt(),track->Y());
      }else{
        if(pdgcode==421) fHistPtvsRapidityPromptD0->Fill(track->Pt(),track->Y());
        if(pdgcode==4122) fHistPtvsRapidityPromptLc->Fill(track->Pt(),track->Y());
        if(pdgcode==4132) fHistPtvsRapidityPromptXic0->Fill(track->Pt(),track->Y());
      }

      if(TMath::Abs(track->Y())<0.5){
        if(pdgcode==421) fHistPtInclusiveD0->Fill(track->Pt());
        if(pdgcode==4122) fHistPtInclusiveLc->Fill(track->Pt());
        if(pdgcode==4132) fHistPtInclusiveXic0->Fill(track->Pt());
        if(isbfd){
          if(pdgcode==421) fHistPtFeeddownD0->Fill(track->Pt());
          if(pdgcode==4122) fHistPtFeeddownLc->Fill(track->Pt());
          if(pdgcode==4132) fHistPtFeeddownXic0->Fill(track->Pt());
        }else{
          if(pdgcode==421) fHistPtPromptD0->Fill(track->Pt());
          if(pdgcode==4122) fHistPtPromptLc->Fill(track->Pt());
          if(pdgcode==4132) fHistPtPromptXic0->Fill(track->Pt());
        }
      }
    }
  } //track loop 


  // Post output data.
  PostData(1, fMyOut);
}      

//________________________________________________________________________
void AliAnalysisTaskCharmBaryonsMC::Terminate(Option_t *) 
{
  //
  // Draw result to the screen
  // Called once at the end of the query   
  //

  fMyOut  = dynamic_cast<TList*> (GetOutputData(1));

  Finalize();
}



void AliAnalysisTaskCharmBaryonsMC::Finalize() 
{

  //
  // Finalize
  //

}

void AliAnalysisTaskCharmBaryonsMC::GetHistory(AliMCParticle *part, AliMCEvent *mcevt, Int_t *pdgarray, Int_t *labelarray, Int_t &ngen)
{
  //
  // Get history of a particle
  //

	for(Int_t i=0;i<100;i++){
		pdgarray[i] = -9999;
		labelarray[i] = -9999;
	}
	ngen = 0;

	AliMCParticle *mcprim = part;
	while(mcprim->GetMother()>=0) {
		Int_t lab_prim=mcprim->GetMother();

		AliMCParticle *tmcprim = (AliMCParticle*)mcevt->GetTrack(lab_prim);
		if(!tmcprim) {
			break;
		}
		if((TMath::Abs(tmcprim->PdgCode())<10) || (TMath::Abs(tmcprim->PdgCode())==21)) break;

		mcprim = tmcprim;

		pdgarray[ngen] = mcprim->PdgCode();
		labelarray[ngen] = lab_prim;

		ngen ++;
		if(ngen == 100) break;
	}
}

Bool_t AliAnalysisTaskCharmBaryonsMC::Match(Float_t a, Float_t b){
  //
  // Match float
  //
	if(TMath::Abs(a-b)<0.00001) return kTRUE;
	if(TMath::Abs(a+b)<0.00001) return kTRUE;
	else return kFALSE;
}
Bool_t AliAnalysisTaskCharmBaryonsMC::Match(Int_t a, Int_t b){
  //
  // Match Int_t
  //
	if((a-b) == 0) return kTRUE;
	if((a+b) == 0) return kTRUE;
	else return kFALSE;
}

Bool_t AliAnalysisTaskCharmBaryonsMC::FromBottom(Int_t *hitory){
  //
  // See if the particle is from bottom hadrons
  //
	for(int i=0;i<10;i++){
		if(Match(hitory[i],511)) return kTRUE;
		if(Match(hitory[i],521)) return kTRUE;
		if(Match(hitory[i],531)) return kTRUE;
		if(Match(hitory[i],5122)) return kTRUE;
		if(Match(hitory[i],5232)) return kTRUE;
		if(Match(hitory[i],5132)) return kTRUE;
		if(Match(hitory[i],5332)) return kTRUE;
	}
	return kFALSE;
}
