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

//-----------------------------------------------------------------------
// Example of task (running locally, on AliEn and CAF),
// which provides standard way of calculating acceptance and efficiency
// between different steps of the procedure.
// The ouptut of the task is a AliCFContainer from which the efficiencies
// can be calculated
//-----------------------------------------------------------------------
// Author : R. Vernet, Consorzio Cometa - Catania (it)
//-----------------------------------------------------------------------
// Modification done by X. Lopez - LPC Clermont (fr)
//-----------------------------------------------------------------------

#ifndef ALICFMUONSINGLETASK1_CXX
#define ALICFMUONSINGLETASK1_CXX

#include "AliCFMuonSingleTask1.h"
#include "AliHeader.h"
#include "AliESDHeader.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TH1I.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliCFManager.h"
#include "AliCFCutBase.h"
#include "AliCFContainer.h"
#include "TChain.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliESDMuonTrack.h"
#include "AliESDtrack.h"
#include "AliESDInputHandler.h"
#include "TCanvas.h"

ClassImp(AliCFMuonSingleTask1)

//__________________________________________________________________________
AliCFMuonSingleTask1::AliCFMuonSingleTask1() :
  fReadAODData(0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fHistEventsProcessed(0x0),
  fNevt(0),
  fIsMC(kFALSE)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliCFMuonSingleTask1::AliCFMuonSingleTask1(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fReadAODData(0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fHistEventsProcessed(0x0),
  fNevt(0),
  fIsMC(kFALSE)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFMuonSingleTask1","Calling Constructor");

  fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;

  DefineOutput(1,TH1I::Class());
  DefineOutput(2,AliCFContainer::Class());

}

//___________________________________________________________________________
AliCFMuonSingleTask1& AliCFMuonSingleTask1::operator=(const AliCFMuonSingleTask1& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
    fReadAODData = c.fReadAODData ;
    fCFManager  = c.fCFManager;
    fQAHistList = c.fQAHistList ;
    fHistEventsProcessed = c.fHistEventsProcessed;
    fNevt = c.fNevt ;
  }
  return *this;
}

//___________________________________________________________________________
AliCFMuonSingleTask1::AliCFMuonSingleTask1(const AliCFMuonSingleTask1& c) :
  AliAnalysisTaskSE(c),
  fReadAODData(c.fReadAODData),
  fCFManager(c.fCFManager),
  fQAHistList(c.fQAHistList),
  fHistEventsProcessed(c.fHistEventsProcessed),
  fNevt(c.fNevt),
  fIsMC(kFALSE)
{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliCFMuonSingleTask1::~AliCFMuonSingleTask1() {
  //
  //destructor
  //
  Info("~AliCFMuonSingleTask1","Calling Destructor");
  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
  if (fQAHistList) {fQAHistList->Clear(); delete fQAHistList;}
}

//_________________________________________________
void AliCFMuonSingleTask1::UserExec(Option_t *)
{
  //
  // Main loop function
  //  

  Info("UserExec"," ") ;
  if (!fMCEvent && fIsMC) {  
    Error("UserExec","NO MC EVENT FOUND!");
    return;
  }

  fNevt++;
  fCFManager->SetEventInfo(fMCEvent);
  Double_t containerInput[18] ;

////////
//// MC
////////

// loop on the MC part

  if(fIsMC){

      for (Int_t ipart=0; ipart<fMCEvent->GetNumberOfTracks(); ipart++) { 
	  
	  AliMCParticle *mcPart  = (AliMCParticle*) fMCEvent->GetTrack(ipart);
	  TParticle *part = mcPart->Particle();

	  // rapidity and Pt cuts (default -4<y<-2.5 et 0<pt<20)
	  if (!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mcPart)) continue;

	  // Selection of mu
	  Int_t pdg=part->GetPdgCode();
	  if(TMath::Abs(pdg)==13){

	      Float_t emc = part->Energy();
	      Float_t pzmc = part->Pz();           
	      Float_t rapmc = Rap(emc,pzmc);
	      Float_t phimc = part->Phi(); 
	      phimc = Phideg(phimc);
	      Float_t ptmc = part->Pt();
	      Float_t pmc = part->P();    
	      Float_t zmc = part->Vz();
	      Float_t etamc = part->Eta();
	      Float_t chargemc=0;
	      if(pdg==13) chargemc=-1;
	      if(pdg==-13) chargemc=1;
	      Float_t thetamc;
	      thetamc = TMath::Abs(-TMath::Pi()+TMath::ASin(ptmc/pmc)); // convention for EVGEN
	      Float_t ymc = part->Vy();
	      Float_t xmc = part->Vx();
	      Float_t rmc = TMath::Sqrt(xmc*xmc+ymc*ymc);   
      
	      containerInput[0] = etamc ;   
	      containerInput[1] = rapmc ;   
	      containerInput[2] = phimc ;   
	      containerInput[3] = ptmc ;   
	      containerInput[4] = pmc ;
	      containerInput[14] = zmc ;
	      containerInput[16] = chargemc ;
	      containerInput[17] = thetamc ;

	      // 10 var calculated only for ESD .i.e set at 1 in MC step
	      for(Int_t i=5; i<=13; i++){
		  containerInput[i] = 1;
	      }

	      containerInput[7] = 3;
	      containerInput[10] = rmc ; // radius at production      
	      containerInput[15] = 175. ;  // rabsmc
	      
	      // fill the container at the first step
	      fCFManager->GetParticleContainer()->Fill(containerInput,0);
	  }
      }
  }
  	      

////////
//// ESD
////////

  Int_t trig1MuL = 1 ,trig1MuH = 2, trigUSL = 4, trigUSH = 8, trigLSL = 16, trigLSH = 32;
  Int_t trig = 0;
 
  AliESDEvent *fESD; 
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>
      (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  fESD = esdH->GetEvent();

// trigger 

  if (fESD->GetTriggerMask() & trig1MuL) trig = 1;
  if (fESD->GetTriggerMask() & trig1MuH) trig = 2;
  if (fESD->GetTriggerMask() & trigUSL)  trig = 3;
  if (fESD->GetTriggerMask() & trigUSH)  trig = 4;
  if (fESD->GetTriggerMask() & trigLSL)  trig = 5;
  if (fESD->GetTriggerMask() & trigLSH)  trig = 6;
 
// vertex 

  Float_t vx = -200 , vy = -200 , vz = -200, vt=-200;
  
  Double_t Ncont = 0.; 
  AliESDVertex* vertex = (AliESDVertex*) fESD->GetVertex();
  if (vertex->GetNContributors()) {
      vz = vertex->GetZv();
      vy = vertex->GetYv();
      vx = vertex->GetXv();
      vt=TMath::Sqrt(vx*vx+vy*vy);
      Ncont=vertex->GetNContributors();
  }

// tracks

  Int_t mult1 = fESD->GetNumberOfMuonTracks() ;

    for (Int_t j = 0; j<mult1; j++) {
	AliESDMuonTrack* mu1 = new AliESDMuonTrack(*(fESD->GetMuonTrack(j)));

	if (!mu1->ContainTrackerData()) { mu1=0; continue; }

	Float_t charger = mu1->Charge();
	Float_t er = mu1->E();
	Float_t pzr = mu1->Pz();
	Float_t rapr=Rap(er,pzr);
	Float_t phir = mu1->Phi(); 
	phir = Phideg(phir);
	Float_t ptr = mu1->Pt();
	Float_t pr = mu1->P();
	Float_t hit = mu1->GetNHit();	
	Float_t chi2 = mu1->GetChi2() / (2.0 * hit - 5);
	Int_t matchtrig = mu1->GetMatchTrigger();
	Float_t chi2match = mu1->GetChi2MatchTrigger();
	Float_t dcar = mu1->GetDCA();
	Float_t zr = mu1->GetZ();
	Float_t etar = mu1->Eta();
	Float_t rabs = mu1->GetRAtAbsorberEnd();
	Float_t thetar ;
	thetar = TMath::Abs(-TMath::Pi()+TMath::ASin(ptr/pr)); 

// rapidity and Pt cuts (default -4<y<-2.5 et 0<pt<20)
	if (!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mu1)) continue;

        Float_t thetaAbs = 180.*(1.-TMath::ATan(rabs/505.)/TMath::Pi());
        if (thetaAbs<171 || thetaAbs>178) continue;

	containerInput[0] = etar ;   
	containerInput[1] = rapr ;   
	containerInput[2] = phir ;   
	containerInput[3] = ptr ;
	containerInput[4] = pr ;
	containerInput[5] = hit ;   
	containerInput[6] = chi2 ;   
	containerInput[7] = matchtrig ;   
	containerInput[8] = chi2match ;
	containerInput[9] = Ncont ;   
	containerInput[10] = vt ;   
	containerInput[11] = vz ;   
	containerInput[12] = trig ;
	containerInput[13] = dcar ;
	containerInput[14] = zr ;
	containerInput[15] = thetaAbs;
	containerInput[16] = charger ;
	containerInput[17] = thetar ;

// fill the container at the second step (for simu, first for data)
	if(fIsMC){
	    fCFManager->GetParticleContainer()->Fill(containerInput,1);
	}
	else{
	    fCFManager->GetParticleContainer()->Fill(containerInput,0);
	}
	

    }      // mu Loop

//  ----------
  fHistEventsProcessed->Fill(0);
  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;

}
//________________________________________________________________________
Float_t AliCFMuonSingleTask1::Rap(Float_t e, Float_t pz) 
{
// calculate rapidity
    Float_t rap;
    if(e!=pz){
	rap = 0.5*TMath::Log((e+pz)/(e-pz));
	return rap;
    }
    else{
	rap = -200;
	return rap;
    }
}
//________________________________________________________________________
Float_t AliCFMuonSingleTask1::Phideg(Float_t phi) 
{
// calculate Phi in range [-180,180] 
    Float_t phideg;
    
	phideg = phi-TMath::Pi();
	phideg = phideg*57.296;
	return phideg;
}
//________________________________________________________________________
void AliCFMuonSingleTask1::Terminate(Option_t *) 
{
  // project pt (var[3]) from the two steps MC(0) and ESD(1)

    AliCFContainer *cont = dynamic_cast<AliCFContainer*> (GetOutputData(2));

    TH1D *kpt = cont->ShowProjection(16,0);
    TH1D *rpt = cont->ShowProjection(16,1);

    TCanvas *c1 = new TCanvas("AliCFMuonSingleTask1"," MC & ESD",10,10,510,510);
    c1->Divide(1,2);
    c1->cd(1);
    kpt->Draw("HIST");
    c1->cd(2);
    rpt->Draw("HIST");

}
//________________________________________________________________________

#endif
