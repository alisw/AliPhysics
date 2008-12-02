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
//
//*-- Yves Schutz (CERN) & Gustavo Conesa Balbastre (INFN-LNF)
//////////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TFile.h> 
#include <TNtuple.h>
#include <TVector3.h> 

#include "AliAnaCaloTrigger.h" 
#include "AliAnalysisManager.h"
#include "AliESDEvent.h" 
#include "AliLog.h"
#include "AliESDCaloCluster.h"

//______________________________________________________________________________
AliAnaCaloTrigger::AliAnaCaloTrigger() :  
  fChain(0),
  fESD(0), 
  fOutputContainer(0),
  fCalorimeter("PHOS"),
  fNtTrigger22(0), 
  fNtTriggerNN(0)

{
  // Default Constructor.

}

//______________________________________________________________________________
AliAnaCaloTrigger::AliAnaCaloTrigger(const char *name) : 
  AliAnalysisTask(name,"AnaCaloTrigger"),  
  fChain(0),
  fESD(0), 
  fOutputContainer(0),
  fCalorimeter("PHOS"),
  fNtTrigger22(0), 
  fNtTriggerNN(0)

{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}
//____________________________________________________________________________
AliAnaCaloTrigger::AliAnaCaloTrigger(const AliAnaCaloTrigger & ct) : 
  AliAnalysisTask(ct), fChain(ct.fChain), fESD(ct.fESD),
  fOutputContainer(ct.fOutputContainer), fCalorimeter(ct. fCalorimeter),
  fNtTrigger22(ct.fNtTrigger22), fNtTriggerNN(ct.fNtTriggerNN)
{

  // cpy ctor
  SetName (ct.GetName()) ; 
  SetTitle(ct.GetTitle()) ;
 
}

//_________________________________________________________________________
AliAnaCaloTrigger & AliAnaCaloTrigger::operator = (const AliAnaCaloTrigger & source)
{
  // assignment operator
  
  if(&source == this) return *this;

  fChain = source.fChain ; 
  fESD = source.fESD ;
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
  //fOutputContainer->Clear() ; 
  //delete fOutputContainer ;

}


//______________________________________________________________________________
void AliAnaCaloTrigger::ConnectInputData(const Option_t*)
{
  // Initialisation of branch container and histograms 
    
  AliInfo(Form("*** Initialization of %s", GetName())) ; 
  
  // Get input data
  fChain = dynamic_cast<TChain *>(GetInputData(0)) ;
  if (!fChain) {
    AliError(Form("Input 0 for %s not found\n", GetName()));
    return ;
  }
  
  fESD = new AliESDEvent();
  fESD->ReadFromTree(fChain);

}

//________________________________________________________________________

void AliAnaCaloTrigger::CreateOutputObjects()
{  

  // Create the outputs containers
  OpenFile(0) ;

  // create histograms 
  fNtTrigger22 = new TNtuple(fCalorimeter+"trigger22", "Trigger data 2x2 patch", "a22:a220:enMax:phEnMax:eta22:phi22:etaMax:phiMax:phEtaMax:phPhiMax");
  fNtTriggerNN = new TNtuple(fCalorimeter+"triggerNN", "Trigger data NxN patch", "aNN:aNN0:enMax:phEnMax:etaNN:phiNN:etaMax:phiMax:phEtaMax:phPhiMax");
  
  // create output container
  
  fOutputContainer = new TObjArray(2) ; 
  fOutputContainer->SetName(GetName()) ; 
  
  fOutputContainer->AddAt(fNtTrigger22,             0) ; 
  fOutputContainer->AddAt(fNtTriggerNN,             1) ; 

}

//______________________________________________________________________________
void AliAnaCaloTrigger::Exec(Option_t *) 
{
  // Processing of one event
  
  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 
  
  // Get trigger information of fCalorimeter 
  TArrayF * triggerAmplitudes = 0x0 ;
  TArrayF * triggerPosition   = 0x0 ;
  Int_t numberOfCaloClusters  =  fESD->GetNumberOfCaloClusters() ;

  if(fCalorimeter == "PHOS"){
    triggerAmplitudes      = fESD->GetPHOSTriggerAmplitudes();
    triggerPosition        = fESD->GetPHOSTriggerPosition();
  }
  else if(fCalorimeter == "EMCAL"){
    triggerAmplitudes    = fESD->GetEMCALTriggerAmplitudes();
    triggerPosition      = fESD->GetEMCALTriggerPosition();
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

  Int_t      icaloCluster ; 
  
  // loop over the Calorimeters Clusters
  
  for(icaloCluster = 0 ; icaloCluster < numberOfCaloClusters ; icaloCluster++) {
    
    AliESDCaloCluster * cluster = fESD->GetCaloCluster(icaloCluster) ;
    
    if (cluster && ( (fCalorimeter == "PHOS" && cluster->IsPHOS())  ||  
		     (fCalorimeter == "EMCAL" && cluster->IsEMCAL()))) {
	  
      Float_t cluEnergy = cluster->E() ; 
      Float_t pos[3] ;
      TVector3 vpos ;
      
      cluster->GetPosition( pos ) ;
      
      if ( cluEnergy > enMax) { 
	enMax = cluEnergy ; 
	vpos.SetXYZ(pos[0], pos[1], pos[2]) ; 
	etaMax = vpos.Eta() ; 
	phiMax = vpos.Phi() ; 
      }

      Double_t * pid = cluster->GetPid() ;
      
      if(pid[AliPID::kPhoton] > 0.9) {
	if ( cluEnergy > phEnMax) { 
	  phEnMax = cluEnergy ; 
	  vpos.SetXYZ(pos[0], pos[1], pos[2]) ; 
	  phEtaMax = vpos.Eta() ; 
	  phPhiMax = vpos.Phi() ; 
	}
      }
    }//if cluster
    
    fNtTrigger22->Fill(ka22, ka22O, enMax, phEnMax, eta22, phi22, etaMax, phiMax * TMath::RadToDeg() + 360., phEtaMax, phPhiMax * TMath::RadToDeg() + 360.);
    fNtTriggerNN->Fill(kaNN, kaNNO, enMax, phEnMax, etaNN, phiNN, etaMax, phiMax * TMath::RadToDeg() + 360., phEtaMax, phPhiMax * TMath::RadToDeg() + 360.);
  }//CaloCluster loop
  
  }//If trigger arrays filled
  
  PostData(0, fOutputContainer);
  
}

//______________________________________________________________________________
//void AliAnaCaloTrigger::Terminate(Option_t *) const
//{
//  // Processing when the event loop is ended
//
//}
