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
// An analysis task to check the PHOS photon data in simulated data
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TFile.h> 
#include <TNtuple.h>
#include <TVector3.h> 

#include "AliTriggerPHOS.h" 
#include "AliESD.h" 
#include "AliLog.h"

//______________________________________________________________________________
AliTriggerPHOS::AliTriggerPHOS(const char *name) : 
  AliAnalysisTask(name,""),  
  fChain(0),
  fESD(0), 
  fOutputContainer(0),
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
AliTriggerPHOS::~AliTriggerPHOS()
{
  // dtor
  fOutputContainer->Clear() ; 
  delete fOutputContainer ;
  delete fNtTrigger22 ; 
  delete fNtTriggerNN ; 
}


//______________________________________________________________________________
void AliTriggerPHOS::ConnectInputData(const Option_t*)
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

void AliTriggerPHOS::CreateOutputObjects()
{  

  // create histograms 
  fNtTrigger22 = new TNtuple("PHOStrigger22", "Trigger data 2x2 patch", "a22:a220:enMax:phEnMax:p22:eta22:phi22:etaMax:phiMax:phEtaMax:phPhiMax");
  fNtTriggerNN = new TNtuple("PHOStriggerNN", "Trigger data NxN patch", "aNN:aNN0:enMax:phEnMax:pNN:etaNN:phiNN:etaMax:phiMax:phEtaMax:phPhiMax");

  // create output container
  
  fOutputContainer = new TObjArray(2) ; 
  fOutputContainer->SetName(GetName()) ; 

  fOutputContainer->AddAt(fNtTrigger22,             0) ; 
  fOutputContainer->AddAt(fNtTriggerNN,             1) ; 

}

//______________________________________________________________________________
void AliTriggerPHOS::Exec(Option_t *) 
{
  // Processing of one event
  
  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 
  
  // ************************  PHOS *************************************
  // Get trigger information 
  
  // trigger amplitudes
  const TArrayF * triggerAmplitudes      = fESD->GetPHOSTriggerAmplitudes();
  const Float_t a22    = static_cast<Float_t>(triggerAmplitudes->At(0)) ; 
  const Float_t a22O   = static_cast<Float_t>(triggerAmplitudes->At(1)) ; 
  const Float_t aNN    = static_cast<Float_t>(triggerAmplitudes->At(2)) ; 
  const Float_t aNNO   = static_cast<Float_t>(triggerAmplitudes->At(3)) ; 

  // trigger position
  const TArrayI * triggerCells      = fESD->GetPHOSTriggerCells();
  const Float_t p22    =  static_cast<Float_t>(triggerCells->At(0)) ; 
  const Float_t phi22  =  static_cast<Float_t>(triggerCells->At(1)) ;
  const Float_t eta22  =  static_cast<Float_t>(triggerCells->At(2)) ;
  const Float_t pNN    =  static_cast<Float_t>(triggerCells->At(3)) ; 
  const Float_t phiNN  =  static_cast<Float_t>(triggerCells->At(4)) ;
  const Float_t etaNN  =  static_cast<Float_t>(triggerCells->At(5)) ; 
  
  Int_t       firstPhosCluster       = fESD->GetFirstPHOSCluster() ;
  const Int_t numberOfPhosClusters   = fESD->GetNumberOfPHOSClusters() ;
   
  Float_t enMax       = 0. ;
  Float_t phEnMax     = 0. ;
  Float_t etaMax      = 999. ;
  Float_t phiMax      = 999. ; 
  Float_t phEtaMax    = 999. ;
  Float_t phPhiMax    = 999. ; 

  Int_t      phosCluster ; 
  
  // loop over the PHOS Cluster
  
  for(phosCluster = firstPhosCluster ; phosCluster < firstPhosCluster + numberOfPhosClusters ; phosCluster++) {
    AliESDCaloCluster * caloCluster = fESD->GetCaloCluster(phosCluster) ;
    if (caloCluster) {

      Float_t cluEnergy = caloCluster->GetClusterEnergy() ; 
      Float_t pos[3] ;
      TVector3 vpos ;
      
      caloCluster->GetGlobalPosition( pos ) ;
      
      if ( cluEnergy > enMax) { 
	enMax = cluEnergy ; 
	vpos.SetXYZ(pos[0], pos[1], pos[2]) ; 
	etaMax = vpos.Eta() ; 
	phiMax = vpos.Phi() ; 
      }

      Float_t * pid = caloCluster->GetPid() ;
      
      if(pid[AliPID::kPhoton] > 0.9) {
	if ( cluEnergy > phEnMax) { 
	  phEnMax = cluEnergy ; 
	  vpos.SetXYZ(pos[0], pos[1], pos[2]) ; 
	  phEtaMax = vpos.Eta() ; 
	  phPhiMax = vpos.Phi() ; 
	}
      }
    }
    
    fNtTrigger22->Fill(a22, a22O, enMax, phEnMax, p22, eta22, phi22, etaMax, phiMax, phEtaMax, phPhiMax );
    fNtTriggerNN->Fill(aNN, aNNO, enMax, phEnMax, pNN, etaNN, phiNN, etaMax, phiMax, phEtaMax, phPhiMax );
 }
  
  
  PostData(0, fOutputContainer);
  
}

//______________________________________________________________________________
void AliTriggerPHOS::Terminate(Option_t *)
{
  // Processing when the event loop is ended

}
