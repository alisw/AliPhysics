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
// An analysis task to check the PHOS/EMCAL simulated trigger
//
//*-- Yves Schutz & Gustavo Conesa Balbastre
//////////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TFile.h> 
#include <TNtuple.h>
#include <TVector3.h> 

#include "AliAnaCaloTrigger.h" 
#include "AliESD.h" 
#include "AliLog.h"

//______________________________________________________________________________
AliAnaCaloTrigger::AliAnaCaloTrigger(const char *name) : 
  AliAnalysisTask(name,""),  
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

//______________________________________________________________________________
AliAnaCaloTrigger::~AliAnaCaloTrigger()
{
  // dtor
  fOutputContainer->Clear() ; 
  delete fOutputContainer ;
  delete fNtTrigger22 ; 
  delete fNtTriggerNN ; 
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
  
  // One should first check if the branch address was taken by some other task
  char ** address = (char **)GetBranchAddress(0, "ESD");
  if (address) {
    fESD = (AliESD*)(*address);
  } else {
    fESD = new AliESD();
    SetBranchAddress(0, "ESD", &fESD);
  }
}

//________________________________________________________________________

void AliAnaCaloTrigger::CreateOutputObjects()
{  

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
  Int_t firstCaloCluster      = 0 ;
  Int_t numberOfCaloClusters  = 0 ;

  if(fCalorimeter == "PHOS"){
    triggerAmplitudes      = fESD->GetPHOSTriggerAmplitudes();
    triggerPosition        = fESD->GetPHOSTriggerPosition();
    firstCaloCluster       = fESD->GetFirstPHOSCluster() ;
    numberOfCaloClusters   = fESD->GetNumberOfPHOSClusters() ;
  }
  else if(fCalorimeter == "EMCAL"){
    triggerAmplitudes    = fESD->GetEMCALTriggerAmplitudes();
    triggerPosition      = fESD->GetEMCALTriggerPosition();
    firstCaloCluster     = fESD->GetFirstEMCALCluster() ;
    numberOfCaloClusters = fESD->GetNumberOfEMCALClusters() ;
  }
  
  // trigger amplitudes
  const Float_t a22    = static_cast<Float_t>(triggerAmplitudes->At(0)) ; 
  const Float_t a22O   = static_cast<Float_t>(triggerAmplitudes->At(1)) ; 
  const Float_t aNN    = static_cast<Float_t>(triggerAmplitudes->At(2)) ; 
  const Float_t aNNO   = static_cast<Float_t>(triggerAmplitudes->At(3)) ; 

  // trigger position
  const Float_t x22  =  static_cast<Float_t>(triggerPosition->At(0)) ; 
  const Float_t y22  =  static_cast<Float_t>(triggerPosition->At(1)) ;
  const Float_t z22  =  static_cast<Float_t>(triggerPosition->At(2)) ;
  const Float_t xNN  =  static_cast<Float_t>(triggerPosition->At(3)) ; 
  const Float_t yNN  =  static_cast<Float_t>(triggerPosition->At(4)) ;
  const Float_t zNN  =  static_cast<Float_t>(triggerPosition->At(5)) ; 
  
 
   
  Float_t enMax       = 0. ;
  Float_t phEnMax     = 0. ;
  Float_t etaMax      = 0.5 ;
  Float_t phiMax      = 0. ; 
  Float_t phEtaMax    = 0.5 ;
  Float_t phPhiMax    = 0. ; 
  
  TVector3 vpos22(x22, y22, z22) ;
  TVector3 vposNN(xNN, yNN, zNN) ;
  Float_t eta22 = vpos22.Eta() ; 
  Float_t phi22 = vpos22.Phi() * TMath::RadToDeg() + 360. ; 
  Float_t etaNN = vposNN.Eta() ; 
  Float_t phiNN = vposNN.Phi() * TMath::RadToDeg() + 360. ; 

  Int_t      icaloCluster ; 
  
  // loop over the Calorimeters Clusters
  
  for(icaloCluster = firstCaloCluster ; icaloCluster < firstCaloCluster + numberOfCaloClusters ; icaloCluster++) {
    AliESDCaloCluster * cluster = fESD->GetCaloCluster(icaloCluster) ;
    if (cluster) {

      Float_t cluEnergy = cluster->GetClusterEnergy() ; 
      Float_t pos[3] ;
      TVector3 vpos ;
      
      cluster->GetGlobalPosition( pos ) ;
      
      if ( cluEnergy > enMax) { 
	enMax = cluEnergy ; 
	vpos.SetXYZ(pos[0], pos[1], pos[2]) ; 
	etaMax = vpos.Eta() ; 
	phiMax = vpos.Phi() ; 
      }

      Float_t * pid = cluster->GetPid() ;
      
      if(pid[AliPID::kPhoton] > 0.9) {
	if ( cluEnergy > phEnMax) { 
	  phEnMax = cluEnergy ; 
	  vpos.SetXYZ(pos[0], pos[1], pos[2]) ; 
	  phEtaMax = vpos.Eta() ; 
	  phPhiMax = vpos.Phi() ; 
	}
      }
    }
    
    fNtTrigger22->Fill(a22, a22O, enMax, phEnMax, eta22, phi22, etaMax, phiMax * TMath::RadToDeg() + 360., phEtaMax, phPhiMax * TMath::RadToDeg() + 360.);
    fNtTriggerNN->Fill(aNN, aNNO, enMax, phEnMax, etaNN, phiNN, etaMax, phiMax * TMath::RadToDeg() + 360., phEtaMax, phPhiMax * TMath::RadToDeg() + 360.);
 }
  
  
  PostData(0, fOutputContainer);
  
}

//______________________________________________________________________________
void AliAnaCaloTrigger::Terminate(Option_t *)
{
  // Processing when the event loop is ended

}
