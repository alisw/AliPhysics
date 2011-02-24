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


#ifndef ALICFMUONRESTASK1_CXX
#define ALICFMUONRESTASK1_CXX

#include "AliCFMuonResTask1.h"
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
  
  Info("UserExec"," ") ;
  if (!fMCEvent) {
    Error("UserExec","NO MC EVENT FOUND!");
    return;
  }

  fNevt++; 
  Double_t containerInput[13] ;
  Double_t beamEnergy=7000;
 
////////
//// MC
//////// 

  fCFManager->SetMCEventInfo(fMCEvent);  
  AliStack* stack = fMCEvent->Stack();

  // loop on the MC event
  for (Int_t ipart=0; ipart<fMCEvent->GetNumberOfTracks(); ipart++) { 
    AliMCParticle *mcPart  = (AliMCParticle*) fMCEvent->GetTrack(ipart);
 
    TParticle *part = mcPart->Particle(); 
 
    // Selection of the resonance
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) continue;

    // Mother kinematics
    Float_t e = part->Energy();
    Float_t pz = part->Pz();           
    //Float_t py = part->Py();
    //Float_t px = part->Px();
    Float_t rapmc = Rap(e,pz);
    //Float_t mass = part->GetCalcMass();

    // Decays kinematics

    Int_t p0 = part->GetDaughter(0);
    TParticle *part0 = stack->Particle(p0); 
   // selection of the rapidity for first muon
    AliMCParticle *mcpart0 = new AliMCParticle(part0);
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mcpart0)) continue;
    Int_t pdg0 = part0->GetPdgCode();

    Int_t p1 = part->GetDaughter(1);
    TParticle *part1 = stack->Particle(p1);
   // selection of the rapidity for second muon
    AliMCParticle *mcpart1 = new AliMCParticle(part1);
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mcpart1)) continue;
    Int_t pdg1 = part1->GetPdgCode();
 
   // 0 mu- 1 mu+
    Float_t e0=-999, pz0=-999, py0=-999, px0=-999, phi0=-999, rapmc0=-999;
    Float_t e1=-999, pz1=-999, py1=-999, px1=-999, phi1=-999, rapmc1=-999;
    Double_t charge0=-999, charge1=-999;

    // ordering sign: first = mu-
    if(pdg0==13){
	e0 = part0->Energy();
	pz0 = part0->Pz();
	py0 = part0->Py();
	px0 = part0->Px();
	phi0 = part0->Phi(); 
	phi0 = Phideg(phi0);    
	rapmc0=Rap(e0,pz0);
	charge0 = part0->GetPDG()->Charge()/3;

	e1 = part1->Energy();
	pz1 = part1->Pz();
	py1 = part1->Py();
	px1 = part1->Px();
	phi1 = part1->Phi();
	phi1 = Phideg(phi1);
	rapmc1=Rap(e1,pz1);
	charge1 = part1->GetPDG()->Charge()/3;
    }
    else{
	if(pdg0==-13){
	    e1 = part0->Energy();
	    pz1 = part0->Pz();
	    py1 = part0->Py();
	    px1 = part0->Px();
	    phi1 = part0->Phi();
	    phi1 = Phideg(phi1);    
	    rapmc1=Rap(e1,pz1);
	    charge1 = part0->GetPDG()->Charge()/3;
	    
	    e0 = part1->Energy();
	    pz0 = part1->Pz();
	    py0 = part1->Py();
	    px0 = part1->Px();
	    phi0 = part1->Phi();
	    phi0 = Phideg(phi0);
	    rapmc0=Rap(e0,pz0); 
	    charge0 = part1->GetPDG()->Charge()/3;
	}
    }
    
    if(pdg0==13 || pdg1==13) { 

	Float_t pmc = TMath::Sqrt((px0+px1)*(px0+px1)+(py0+py1)*(py0+py1)+
				   (pz0+pz1)*(pz0+pz1));
	Float_t ptmc = TMath::Sqrt((px0+px1)*(px0+px1)+(py0+py1)*(py0+py1));
	Float_t imassmc = Imass(e0,px0,py0,pz0,e1,px1,py1,pz1);
	//Float_t rapmc_check=Rap((e0+e1),(pz0+pz1));
	
	Double_t costCS = CostCS((Double_t) px0,(Double_t) py0,(Double_t) pz0,(Double_t) e0, (Double_t)charge0,(Double_t) px1,(Double_t) py1,(Double_t) pz1,(Double_t) e1, beamEnergy);
	Double_t costHE = CostHE((Double_t) px0,(Double_t) py0,(Double_t) pz0,(Double_t) e0, (Double_t)charge0,(Double_t) px1,(Double_t) py1,(Double_t) pz1,(Double_t) e1);
	Double_t phiCS = PhiCS((Double_t) px0,(Double_t) py0,(Double_t) pz0,(Double_t) e0, (Double_t)charge0,(Double_t) px1,(Double_t) py1,(Double_t) pz1,(Double_t) e1, beamEnergy);
	Double_t phiHE = PhiHE((Double_t) px0,(Double_t) py0,(Double_t) pz0,(Double_t) e0, (Double_t)charge0,(Double_t) px1,(Double_t) py1,(Double_t) pz1,(Double_t) e1, beamEnergy);

	containerInput[0] = fNevt+0.5 ;   
	containerInput[1] = rapmc0 ;   
	containerInput[2] = phi0 ;   
	containerInput[3] = rapmc1 ;   
	containerInput[4] = phi1 ;   
	containerInput[5] = imassmc ;   
	containerInput[6] = rapmc ;   
	containerInput[7] = ptmc;
	containerInput[8] = pmc;	
	containerInput[9] = costCS;	
	containerInput[10] = phiCS;	
	containerInput[11] = costHE;	
	containerInput[12] = phiHE;	

	// fill the container at the first step
	fCFManager->GetParticleContainer()->Fill(containerInput,0);
    }
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
	    Float_t phir1 = mu1->Phi(); 
	    phir1=Phideg(phir1);
 	    Float_t er1 = mu1->E();
	    Float_t rap1=Rap(er1,pzr1);
	    // selection of the rapidity for first muon
	    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mu1)) continue;

	    for (Int_t jj = 0; jj<mult1; jj++) {
		AliESDMuonTrack* mu2 = new AliESDMuonTrack(*(fESD->GetMuonTrack(jj)));
		Float_t zr2 = mu2->Charge();
// Select mu+
		if(zr2>0 && jj !=j){
		    Float_t pxr2 = mu2->Px();
		    Float_t pyr2 = mu2->Py();
		    Float_t pzr2 = mu2->Pz();
		    Float_t phir2 = mu2->Phi();
		    phir2 = Phideg(phir2);
		    Float_t er2 = mu2->E();
		    Float_t rap2=Rap(er2,pzr2);
		    // selection of the rapidity for second muon
		    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mu2)) continue;

		    Float_t prec = TMath::Sqrt((pxr1+pxr2)*(pxr1+pxr2)+(pyr1+pyr2)*(pyr1+pyr2)+
						(pzr1+pzr2)*(pzr1+pzr2));
		    Float_t ptrec = TMath::Sqrt((pxr1+pxr2)*(pxr1+pxr2)+(pyr1+pyr2)*(pyr1+pyr2));
		    Float_t raprec=Rap((er1+er2),(pzr1+pzr2));
		    Float_t imassrec = Imass(er1,pxr1,pyr1,pzr1,er2,pxr2,pyr2,pzr2);
		    
		    Double_t costCSrec = CostCS((Double_t) pxr1,(Double_t) pyr1,(Double_t)pzr1,(Double_t) er1, (Double_t)zr1,(Double_t) pxr2,(Double_t) pyr2,(Double_t)pzr2,(Double_t) er2, beamEnergy);
		    Double_t costHErec = CostHE((Double_t) pxr1,(Double_t) pyr1,(Double_t)pzr1,(Double_t) er1, (Double_t)zr1,(Double_t) pxr2,(Double_t) pyr2,(Double_t)pzr2,(Double_t) er2);
		    Double_t phiCSrec = PhiCS((Double_t) pxr1,(Double_t) pyr1,(Double_t)pzr1,(Double_t) er1, (Double_t)zr1,(Double_t) pxr2,(Double_t) pyr2,(Double_t)pzr2,(Double_t) er2, beamEnergy);
		    Double_t phiHErec = PhiHE((Double_t) pxr1,(Double_t) pyr1,(Double_t)pzr1,(Double_t) er1, (Double_t)zr1,(Double_t) pxr2,(Double_t) pyr2,(Double_t)pzr2,(Double_t) er2, beamEnergy);

		    containerInput[0] = fNevt+0.5 ;   
		    containerInput[1] = rap1 ;   
		    containerInput[2] = phir1 ;   
		    containerInput[3] = rap2 ;
		    containerInput[4] = phir2 ;   
		    containerInput[5] = imassrec ;   
		    containerInput[6] = raprec ;   
		    containerInput[7] = ptrec;
		    containerInput[8] = prec;
		    containerInput[9] = costCSrec;	
		    containerInput[10] = phiCSrec;	
		    containerInput[11] = costHErec;	
		    containerInput[12] = phiHErec;	
		    
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
Float_t AliCFMuonResTask1::Imass(Float_t e1, Float_t px1, Float_t py1, Float_t pz1,
				   Float_t e2, Float_t px2, Float_t py2, Float_t pz2) const
{
// invariant mass calculation
    Float_t imassrec = TMath::Sqrt((e1+e2)*(e1+e2)-((px1+px2)*(px1+px2)+
                                    (py1+py2)*(py1+py2)+(pz1+pz2)*(pz1+pz2)));
    return imassrec;
}
//________________________________________________________________________
Float_t AliCFMuonResTask1::Rap(Float_t e, Float_t pz) const
{
// calculate rapidity
    Float_t rap;
    if(TMath::Abs(e-pz)>1e-7){
	rap = 0.5*TMath::Log((e+pz)/(e-pz));
	return rap;
    }
    else{
	rap = -200;
	return rap;
    }
}
//________________________________________________________________________
Float_t AliCFMuonResTask1::Phideg(Float_t phi) const
{
// calculate Phi in range [-180,180] 
    Float_t phideg;
    
	phideg = phi-TMath::Pi();
	phideg = phideg*57.296;
	return phideg;
}
//________________________________________________________________________
Double_t AliCFMuonResTask1::CostCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2,
Double_t Energy)
{
// calculate costCS
  TLorentzVector pMu1CM, pMu2CM, pProjCM, pTargCM, pDimuCM; // In the CM. frame
  TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
  TVector3 beta,zaxisCS;
  Double_t mp=0.93827231;
  //
  // --- Fill the Lorentz vector for projectile and target in the CM frame
  //
  pProjCM.SetPxPyPzE(0.,0.,-Energy,TMath::Sqrt(Energy*Energy+mp*mp)); 
  pTargCM.SetPxPyPzE(0.,0.,Energy,TMath::Sqrt(Energy*Energy+mp*mp)); 
  //
  // --- Get the muons parameters in the CM frame 
  //
  pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
  pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  pDimuCM=pMu1CM+pMu2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  beta=(-1./pDimuCM.E())*pDimuCM.Vect();
  pMu1Dimu=pMu1CM;
  pMu2Dimu=pMu2CM;
  pProjDimu=pProjCM;
  pTargDimu=pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the CS angle 
  //
  zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //				     
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  Double_t cost;
  
  if(charge1>0) cost = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
  else cost = zaxisCS.Dot((pMu2Dimu.Vect()).Unit());
  return cost;
}
//________________________________________________________________________

//________________________________________________________________________
Double_t AliCFMuonResTask1::CostHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
// calculate costHE
  TLorentzVector pMu1CM, pMu2CM, pDimuCM; // In the CM frame 
  TLorentzVector pMu1Dimu, pMu2Dimu; // In the dimuon rest frame
  TVector3 beta,zaxisCS;
  //
  // --- Get the muons parameters in the CM frame
  //
  pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
  pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  pDimuCM=pMu1CM+pMu2CM;
  //
  // --- Translate the muon parameters in the dimuon rest frame
  //
  beta=(-1./pDimuCM.E())*pDimuCM.Vect();
  pMu1Dimu=pMu1CM;
  pMu2Dimu=pMu2CM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis;
  zaxis=(pDimuCM.Vect()).Unit();
  
  // --- Calculation of the polarization angle (angle between mu+ and the z axis defined above)
  Double_t cost;
  if(charge1>0) cost = zaxis.Dot((pMu1Dimu.Vect()).Unit()); 
  else cost = zaxis.Dot((pMu2Dimu.Vect()).Unit()); 
  
  return cost;
}
//________________________________________________________________________

//________________________________________________________________________
Double_t AliCFMuonResTask1::PhiCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2,
Double_t Energy)
{
// calculate phiCS
   TLorentzVector pMu1CM, pMu2CM, pProjCM, pTargCM, pDimuCM; // In the CM frame
   TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
   TVector3 beta,yaxisCS, xaxisCS, zaxisCS;
   Double_t mp=0.93827231;
   //
   // --- Fill the Lorentz vector for projectile and target in the CM frame
   //
   pProjCM.SetPxPyPzE(0.,0.,-Energy,TMath::Sqrt(Energy*Energy+mp*mp)); 
   pTargCM.SetPxPyPzE(0.,0.,Energy,TMath::Sqrt(Energy*Energy+mp*mp)); 
   //
   // --- Get the muons parameters in the CM frame 
   //
   pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
   pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
   //
   // --- Obtain the dimuon parameters in the CM frame
   //
   pDimuCM=pMu1CM+pMu2CM;
   //
   // --- Translate the dimuon parameters in the dimuon rest frame
   //
   beta=(-1./pDimuCM.E())*pDimuCM.Vect();
   pMu1Dimu=pMu1CM;
   pMu2Dimu=pMu2CM;
   pProjDimu=pProjCM;
   pTargDimu=pTargCM;
   pMu1Dimu.Boost(beta);
   pMu2Dimu.Boost(beta);
   pProjDimu.Boost(beta);
   pTargDimu.Boost(beta);
   //
   // --- Determine the z axis for the CS angle 
   //
   zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
   yaxisCS=(((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
   xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit();
 
   Double_t phi;
   if(charge1>0) phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
   else phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxisCS),((pMu2Dimu.Vect()).Dot(xaxisCS)));
     
   return phi;
}
//________________________________________________________________________

//________________________________________________________________________
Double_t AliCFMuonResTask1::PhiHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2, Double_t Energy)
{
// calculate phiHE
   TLorentzVector pMu1CM, pMu2CM, pProjCM, pTargCM, pDimuCM; // In the CM frame 
   TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
   TVector3 beta,xaxis,yaxis,zaxis;
   Double_t mp=0.93827231;
 
   //
   // --- Get the muons parameters in the CM frame
   //
   pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
   pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
   //
   // --- Obtain the dimuon parameters in the CM frame
   //
   pDimuCM=pMu1CM+pMu2CM;
   //
   // --- Translate the muon parameters in the dimuon rest frame
   // 
   zaxis=(pDimuCM.Vect()).Unit();
 
   beta=(-1./pDimuCM.E())*pDimuCM.Vect();
 
   pProjCM.SetPxPyPzE(0.,0.,-Energy,TMath::Sqrt(Energy*Energy+mp*mp));
   pTargCM.SetPxPyPzE(0.,0.,Energy,TMath::Sqrt(Energy*Energy+mp*mp)); 
 
   pProjDimu=pProjCM;
   pTargDimu=pTargCM;
 
   pProjDimu.Boost(beta);
   pTargDimu.Boost(beta);
   
   yaxis=((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
   xaxis=(yaxis.Cross(zaxis)).Unit();
   
   pMu1Dimu=pMu1CM;
   pMu2Dimu=pMu2CM;
   pMu1Dimu.Boost(beta);
   pMu2Dimu.Boost(beta);
 
   Double_t phi;
   if(charge1>0) phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
   else phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxis),(pMu2Dimu.Vect()).Dot(xaxis));
   
   return phi;
}

//________________________________________________________________________
void AliCFMuonResTask1::Terminate(Option_t *) 
{
  // draw result of the Invariant mass MC and ESD

    AliCFContainer *cont = dynamic_cast<AliCFContainer*> (GetOutputData(2));   

    TH1D *kmass = cont->ShowProjection(5,0);
    TH1D *rmass = cont->ShowProjection(5,1);
    TH1D *kcostCS = cont->ShowProjection(9,0);
    TH1D *rcostCS = cont->ShowProjection(9,1);

    TCanvas *c1 = new TCanvas("AliCFMuonResTask1","JPSI MC & ESD",10,10,510,510);
    c1->Divide(2,2);
    c1->cd(1);
    kmass->Draw("HIST");
    c1->cd(2);
    rmass->Draw("HIST");
    c1->cd(3);
    kcostCS->Draw("HIST");
    c1->cd(4);
    rcostCS->Draw("HIST");
}
//________________________________________________________________________

#endif
