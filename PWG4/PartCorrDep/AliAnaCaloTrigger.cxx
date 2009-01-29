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

//_________________________________________________________________________
// An analysis task to check the trigger data in ESD
// Creates an ntuple for 2x2 and NxN triggers
// Each ntuple connects the maximum trigger amplitudes 
// and its positions with reconstructed clusters
// and if MC stack available, with pt of parent.
//
//*-- Yves Schutz (CERN) & Gustavo Conesa Balbastre (INFN-LNF)
//////////////////////////////////////////////////////////////////////////////
//Root
#include <TNtuple.h>
#include <TVector3.h> 

//Aliroot
#include "AliAnaCaloTrigger.h" 
#include "AliStack.h"
#include "AliLog.h"
#include "AliESDCaloCluster.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"

//______________________________________________________________________________
AliAnaCaloTrigger::AliAnaCaloTrigger() :  
  fOutputContainer(0),
  fCalorimeter("PHOS"),
  fNtTrigger22(0), 
  fNtTriggerNN(0)

{
  // Default Constructor.

}

//______________________________________________________________________________
AliAnaCaloTrigger::AliAnaCaloTrigger(const char *name) : 
  AliAnalysisTaskSE(name),  
  fOutputContainer(0),
  fCalorimeter("PHOS"),
  fNtTrigger22(0), 
  fNtTriggerNN(0)

{
  // Constructor.
  // Output slot 
  DefineOutput(1,  TList::Class()) ; 
}

//____________________________________________________________________________
AliAnaCaloTrigger::AliAnaCaloTrigger(const AliAnaCaloTrigger & ct) : 
  AliAnalysisTaskSE(ct.GetName()),
  fOutputContainer(ct.fOutputContainer), fCalorimeter(ct. fCalorimeter),
  fNtTrigger22(ct.fNtTrigger22), fNtTriggerNN(ct.fNtTriggerNN)
{

  // cpy ctor
   
}

//_________________________________________________________________________
AliAnaCaloTrigger & AliAnaCaloTrigger::operator = (const AliAnaCaloTrigger & source)
{
  // assignment operator
  
  //if(&source == this) return *this;
  this->~AliAnaCaloTrigger();
  new(this) AliAnaCaloTrigger(source);

  fOutputContainer = source.fOutputContainer ;
  fCalorimeter = source. fCalorimeter ;
  fNtTrigger22 = source.fNtTrigger22 ;
  fNtTriggerNN = source.fNtTriggerNN ;

  return *this;
  
}

//______________________________________________________________________________
AliAnaCaloTrigger::~AliAnaCaloTrigger()
{
  // dtor
	if(fOutputContainer){
		fOutputContainer->Clear() ; 
		delete fOutputContainer ;
	}
}


//________________________________________________________________________

void AliAnaCaloTrigger::UserCreateOutputObjects()
{  

  // Create the outputs containers
  OpenFile(1) ;

  // create histograms 
  fNtTrigger22 = new TNtuple(fCalorimeter+"trigger22", "Trigger data 2x2 patch", "a22:a220:ptGen:enMax:phEnMax:eta22:phi22:etaMax:phiMax:phEtaMax:phPhiMax");
  fNtTriggerNN = new TNtuple(fCalorimeter+"triggerNN", "Trigger data NxN patch", "aNN:aNN0:ptGen:enMax:phEnMax:etaNN:phiNN:etaMax:phiMax:phEtaMax:phPhiMax");
  
  // create output container
  
  fOutputContainer = new TList() ; 
  fOutputContainer->SetName(GetName()) ; 
  
  fOutputContainer->AddAt(fNtTrigger22,             0) ; 
  fOutputContainer->AddAt(fNtTriggerNN,             1) ; 

}

//______________________________________________________________________________
void AliAnaCaloTrigger::UserExec(Option_t *) 
{
	// Processing of one event
	
	if ( !((Entry()-1)%100) ) 
		printf(" Processing event # %lld\n",  Entry()) ; 
	AliESDEvent* esd = (AliESDEvent*)InputEvent();
	
	//Get MC data, if available
	AliStack* stack = 0x0; 
	if(MCEvent())
		stack = MCEvent()->Stack();
	
	// Get trigger information of fCalorimeter 
	TArrayF * triggerAmplitudes = 0x0 ;
	TArrayF * triggerPosition   = 0x0 ;
	Int_t numberOfCaloClusters  =  esd->GetNumberOfCaloClusters() ;
	
	if(fCalorimeter == "PHOS"){
		triggerAmplitudes      = esd->GetPHOSTriggerAmplitudes();
		triggerPosition        = esd->GetPHOSTriggerPosition();
	}
	else if(fCalorimeter == "EMCAL"){
		triggerAmplitudes    = esd->GetEMCALTriggerAmplitudes();
		triggerPosition      = esd->GetEMCALTriggerPosition();
	}
	
	if( triggerAmplitudes && triggerPosition ){
		// trigger amplitudes
		const Float_t ka22    = static_cast<Float_t>(triggerAmplitudes->At(0)) ; 
		const Float_t ka22O   = static_cast<Float_t>(triggerAmplitudes->At(1)) ; 
		const Float_t kaNN    = static_cast<Float_t>(triggerAmplitudes->At(2)) ; 
		const Float_t kaNNO   = static_cast<Float_t>(triggerAmplitudes->At(3)) ; 
		
		// trigger position
		const Float_t kx22  =  static_cast<Float_t>(triggerPosition->At(0)) ; 
		const Float_t ky22  =  static_cast<Float_t>(triggerPosition->At(1)) ;
		const Float_t kz22  =  static_cast<Float_t>(triggerPosition->At(2)) ;
		const Float_t kxNN  =  static_cast<Float_t>(triggerPosition->At(3)) ; 
		const Float_t kyNN  =  static_cast<Float_t>(triggerPosition->At(4)) ;
		const Float_t kzNN  =  static_cast<Float_t>(triggerPosition->At(5)) ; 
		
		//printf("ka22 %f, ka220 %f, kaNN %f, kaNN0 %f\n",ka22,ka22O,kaNN,kaNNO);
		//printf("kx22 %f, ky22 %f, kz22 %f, kxNN %f, kyNN %f, kzNN %f \n",kx22,ky22,kz22,kxNN,kyNN,kzNN);
		
		Float_t enMax       = 0. ;
		Float_t phEnMax     = 0. ;
		Float_t etaMax      = 0.5 ;
		Float_t phiMax      = 0. ; 
		Float_t phEtaMax    = 0.5 ;
		Float_t phPhiMax    = 0. ; 
		
		TVector3 vpos22(kx22, ky22, kz22) ;
		TVector3 vposNN(kxNN, kyNN, kzNN) ;
		Float_t eta22 = vpos22.Eta() ; 
		Float_t phi22 = vpos22.Phi() * TMath::RadToDeg() + 360. ; 
		Float_t etaNN = vposNN.Eta() ; 
		Float_t phiNN = vposNN.Phi() * TMath::RadToDeg() + 360. ; 
		
		
		Int_t      icaloCluster = 0 ; 
		Int_t      labelmax     = -1 ;
		// loop over the Calorimeters Clusters
		
		for(icaloCluster = 0 ; icaloCluster < numberOfCaloClusters ; icaloCluster++) {
			
			AliESDCaloCluster * cluster = esd->GetCaloCluster(icaloCluster) ;
			
			if (cluster && ( (fCalorimeter == "PHOS" && cluster->IsPHOS())  ||  
							(fCalorimeter == "EMCAL" && cluster->IsEMCAL()))) {
				
				Float_t cluEnergy = cluster->E() ; 
				Float_t pos[3] ;
				TVector3 vpos ;
				
				cluster->GetPosition( pos ) ;
				
				if ( cluEnergy > enMax) { 
					enMax    = cluEnergy ; 
					vpos.SetXYZ(pos[0], pos[1], pos[2]) ; 
					etaMax   = vpos.Eta() ; 
					phiMax   = vpos.Phi() ; 
					labelmax = cluster->GetLabel();
				}
				
				Double_t * pid = cluster->GetPid() ;
				
				if(pid[AliPID::kPhoton] > 0.9) {
					if ( cluEnergy > phEnMax) { 
						phEnMax   = cluEnergy ; 
						vpos.SetXYZ(pos[0], pos[1], pos[2]) ; 
						phEtaMax = vpos.Eta() ; 
						phPhiMax = vpos.Phi() ; 
					}
				}
			}//if cluster
			
		    Float_t ptGen = -1;
			if(stack && labelmax < stack->GetNtrack() && labelmax >= 0 ){
				TParticle * particle = stack->Particle(labelmax); 
				ptGen = particle->Energy();
			}
			
			fNtTrigger22->Fill(ka22, ka22O, ptGen, enMax, phEnMax, eta22, phi22, etaMax, phiMax * TMath::RadToDeg() + 360., phEtaMax, phPhiMax * TMath::RadToDeg() + 360.);
			fNtTriggerNN->Fill(kaNN, kaNNO, ptGen, enMax, phEnMax, etaNN, phiNN, etaMax, phiMax * TMath::RadToDeg() + 360., phEtaMax, phPhiMax * TMath::RadToDeg() + 360.);
		
		}//CaloCluster loop
		
	}//If trigger arrays filled
	
	PostData(1, fOutputContainer);
	
}

//______________________________________________________________________________
//void AliAnaCaloTrigger::Terminate(Option_t *) const
//{
//  // Processing when the event loop is ended
//
//}
