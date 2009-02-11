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


#ifndef ALICFMUONRESTASK1_CXX
#define ALICFMUONRESTASK1_CXX

#include "AliCFMuonResTask1.h"
#include "AliHeader.h"
#include "AliESDHeader.h"
#include "AliStack.h"
#include "TParticle.h"
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

ClassImp(AliCFMuonResTask1)

//__________________________________________________________________________
AliCFMuonResTask1::AliCFMuonResTask1() :
  fReadAODData(0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fHistEventsProcessed(0x0),
  fNevt(0)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliCFMuonResTask1::AliCFMuonResTask1(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fReadAODData(0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fHistEventsProcessed(0x0),
  fNevt(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFMuonResTask1","Calling Constructor");

  fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;

  DefineOutput(1,TH1I::Class());
  DefineOutput(2,AliCFContainer::Class());

}

//___________________________________________________________________________
AliCFMuonResTask1& AliCFMuonResTask1::operator=(const AliCFMuonResTask1& c) 
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
AliCFMuonResTask1::AliCFMuonResTask1(const AliCFMuonResTask1& c) :
  AliAnalysisTaskSE(c),
  fReadAODData(c.fReadAODData),
  fCFManager(c.fCFManager),
  fQAHistList(c.fQAHistList),
  fHistEventsProcessed(c.fHistEventsProcessed),
  fNevt(c.fNevt) 
{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliCFMuonResTask1::~AliCFMuonResTask1() {
  //
  //destructor
  //
  Info("~AliCFMuonResTask1","Calling Destructor");
  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
  if (fQAHistList) {fQAHistList->Clear(); delete fQAHistList;}
}

//_________________________________________________
void AliCFMuonResTask1::UserExec(Option_t *)
{
  //
  // Main loop function
  // 
  
  Info("UserExec","") ;
  if (!fMCEvent) {
    Error("UserExec","NO MC EVENT FOUND!");
    return;
  }

  fNevt++; 
  Double_t containerInput[9] ;
 
////////
//// MC
//////// 

  fCFManager->SetEventInfo(fMCEvent);
  AliStack* stack = fMCEvent->Stack();

  // loop on the MC event
  for (Int_t ipart=0; ipart<fMCEvent->GetNumberOfTracks(); ipart++) { 
    AliMCParticle *mcPart  = fMCEvent->GetTrack(ipart);
 
    TParticle *part = mcPart->Particle(); 
    TParticle *part0 = mcPart->Particle();
    TParticle *part1 = mcPart->Particle();

    // Selection of the resonance
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) continue;

    // Mother kinematics
    Float_t e = part->Energy();
    Float_t pz = part->Pz();           
    Float_t py = part->Py();
    Float_t px = part->Px();
    Float_t rapmc = Rap(e,pz);
    Float_t mass = part->GetCalcMass();

    // Decays kinematics
    Int_t p0 = part->GetDaughter(0);
    part0 = stack->Particle(p0);
    Int_t pdg0 = part0->GetPdgCode();
    Float_t e0 = part0->Energy();
    Float_t pz0 = part0->Pz();
    Float_t py0 = part0->Py();
    Float_t px0 = part0->Px();
    Float_t phi0 = part0->Phi(); // Warning in TParticle Phi = pi + ATan2(Py,Px) = [0,2pi] 
    phi0 = Phideg(phi0);    
    Float_t rapmc0=Rap(e0,pz0);
    AliMCParticle *mcpart0 = new AliMCParticle(part0);

    // selection of the rapidity for first muon
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,mcpart0)) continue;                

    Int_t p1 = part->GetDaughter(1);
    part1 = stack->Particle(p1);
    Int_t pdg1 = part1->GetPdgCode();
    Float_t e1 = part1->Energy();
    Float_t pz1 = part1->Pz();
    Float_t py1 = part1->Py();
    Float_t px1 = part1->Px();
    Float_t phi1 = part1->Phi();
    phi1 = Phideg(phi1);
    Float_t rapmc1=Rap(e1,pz1);
    AliMCParticle *mcpart1 = new AliMCParticle(part1);

    // selection of the rapidity for second muon
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,mcpart1)) continue;

    if(pdg0==13 || pdg1==13) { 

	Float_t pmc = TMath::Sqrt((px0+px1)*(px0+px1)+(py0+py1)*(py0+py1)+
				   (pz0+pz1)*(pz0+pz1));
	Float_t ptmc = TMath::Sqrt((px0+px1)*(px0+px1)+(py0+py1)*(py0+py1));
	Float_t imassmc = Imass(e0,px0,py0,pz0,e1,px1,py1,pz1);
	Float_t rapmc=Rap((e0+e1),(pz0+pz1));

	containerInput[0] = fNevt ;   
	containerInput[1] = rapmc0 ;   
	containerInput[2] = phi0 ;   
	containerInput[3] = rapmc1 ;   
	containerInput[4] = phi1 ;   
	containerInput[5] = imassmc ;   
	containerInput[6] = rapmc ;   
	containerInput[7] = ptmc;
	containerInput[8] = pmc;	

	// fill the container at the first step
	fCFManager->GetParticleContainer()->Fill(containerInput,0);

    } // one muon is positive
  }    

////////
//// ESD
////////
  
  AliESDEvent *fESD; 
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>
      (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  fESD = esdH->GetEvent();

  Int_t mult1 = fESD->GetNumberOfMuonTracks() ;

    for (Int_t j = 0; j<mult1; j++) { 
	AliESDMuonTrack* mu1 = new AliESDMuonTrack(*(fESD->GetMuonTrack(j)));
	Float_t zr1 = mu1->Charge();
// Select mu-
	if(zr1<0){
	    Float_t pxr1 = mu1->Px();
	    Float_t pyr1 = mu1->Py();
	    Float_t pzr1 = mu1->Pz();
	    Float_t phir1 = mu1->Phi(); // Warning in AliESDMuonTrack Phi = ATan2(Py,Px) = [-pi,pi]
	    phir1 = phir1 * 57.296;
 	    Float_t er1 = mu1->E();
	    Float_t rap1=Rap(er1,pzr1);
	    // selection of the rapidity for first muon
	    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,mu1)) continue;

	    for (Int_t jj = 0; jj<mult1; jj++) {
		AliESDMuonTrack* mu2 = new AliESDMuonTrack(*(fESD->GetMuonTrack(jj)));
		Float_t zr2 = mu2->Charge();
// Select mu+
		if(zr2>0 && jj !=j){
		    Float_t pxr2 = mu2->Px();
		    Float_t pyr2 = mu2->Py();
		    Float_t pzr2 = mu2->Pz();
		    Float_t phir2 = mu2->Phi();
		    phir2 = phir2 * 57.296;
		    Float_t er2 = mu2->E();
		    Float_t rap2=Rap(er2,pzr2);
		    // selection of the rapidity for second muon
		    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,mu2)) continue;

		    Float_t prec = TMath::Sqrt((pxr1+pxr2)*(pxr1+pxr2)+(pyr1+pyr2)*(pyr1+pyr2)+
						(pzr1+pzr2)*(pzr1+pzr2));
		    Float_t ptrec = TMath::Sqrt((pxr1+pxr2)*(pxr1+pxr2)+(pyr1+pyr2)*(pyr1+pyr2));
		    Float_t raprec=Rap((er1+er2),(pzr1+pzr2));
		    Float_t imassrec = Imass(er1,pxr1,pyr1,pzr1,er2,pxr2,pyr2,pzr2);

		    containerInput[0] = fNevt ;   
		    containerInput[1] = rap1 ;   
		    containerInput[2] = phir1 ;   
		    containerInput[3] = rap2 ;
		    containerInput[4] = phir2 ;   
		    containerInput[5] = imassrec ;   
		    containerInput[6] = raprec ;   
		    containerInput[7] = ptrec;
		    containerInput[8] = prec;
		    
		    // fill the container at the second step
		    fCFManager->GetParticleContainer()->Fill(containerInput,1);

		}  // mu+ Selection
	    }      // second mu Loop
	}          // mu- Selection
    }        

//  ----------
  fHistEventsProcessed->Fill(0);
  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
}
//________________________________________________________________________
const Float_t AliCFMuonResTask1::Imass(Float_t e1, Float_t px1, Float_t py1, Float_t pz1,
				   Float_t e2, Float_t px2, Float_t py2, Float_t pz2) 
{
// invariant mass calculation
    Float_t imassrec = TMath::Sqrt((e1+e2)*(e1+e2)-((px1+px2)*(px1+px2)+
                                    (py1+py2)*(py1+py2)+(pz1+pz2)*(pz1+pz2)));
    return imassrec;
}
//________________________________________________________________________
const Float_t AliCFMuonResTask1::Rap(Float_t e, Float_t pz) 
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
const Float_t AliCFMuonResTask1::Phideg(Float_t phi) 
{
// calculate Phi from TParticle in range [-180,180] 
    Float_t phideg;
    
	phideg = phi-TMath::Pi();
	phideg = phideg*57.296;
	return phideg;
}
//________________________________________________________________________
void AliCFMuonResTask1::Terminate(Option_t *) 
{
  // draw result of the Invariant mass MC and ESD

    AliCFContainer *cont = dynamic_cast<AliCFContainer*> (GetOutputData(2));   

    TH1D *kmass = cont->ShowProjection(5,0);
    TH1D *rmass = cont->ShowProjection(5,1);

    TCanvas *c1 = new TCanvas("AliCFMuonResTask1","JPSI MC & ESD",10,10,510,510);
    c1->Divide(1,2);
    c1->cd(1);
    kmass->Draw("HIST");
    c1->cd(2);
    rmass->Draw("HIST");
}
//________________________________________________________________________

#endif
