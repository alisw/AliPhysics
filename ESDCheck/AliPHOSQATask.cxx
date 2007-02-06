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

//_________________________________________________________________________
// An analysis task to check the PHOS photon data in simulated data
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h> 
#include <TH1.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TLegend.h> 
#include <TNtuple.h>
#include <TROOT.h> 
#include <TVector3.h> 

#include "AliPHOSQATask.h" 
#include "AliESD.h" 
#include "AliLog.h"

//______________________________________________________________________________
AliPHOSQATask::AliPHOSQATask(const char *name) : 
  AliAnalysisTask(name,""),  
  fChain(0),
  fESD(0), 
  fhPHOS(0),
  fhPHOSEnergy(0),
  fhPHOSDigits(0),
  fhPHOSRecParticles(0),
  fhPHOSPhotons(0),
  fhPHOSInvariantMass(0),
  fhPHOSDigitsEvent(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}

//______________________________________________________________________________
AliPHOSQATask::~AliPHOSQATask()
{
  // dtor
  fOutputContainer->Clear() ; 
  delete fOutputContainer ;

  delete fhPHOSPos ;
  delete fhPHOS ;
  delete fhPHOSEnergy ;
  delete fhPHOSDigits ;
  delete fhPHOSRecParticles ;
  delete fhPHOSPhotons ;
  delete fhPHOSInvariantMass ;
  delete fhPHOSDigitsEvent ;
}


//______________________________________________________________________________
void AliPHOSQATask::ConnectInputData(const Option_t*)
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
void AliPHOSQATask::CreateOutputObjects()
{  
  // create histograms 
  
  fhPHOSPos            = new TNtuple("PHOSPos"         , "Position in PHOS"  , "x:y:z");
  fhPHOS               = new TNtuple("PHOS"            , "PHOS"  , "event:digits:clusters:photons");
  fhPHOSEnergy         = new TH1D("PHOSEnergy"         , "PHOSEnergy"        , 1000, 0., 10. ) ;
  fhPHOSDigits         = new TH1I("PHOSDigitsCluster"  , "PHOSDigits"        , 20 , 0 , 20  ) ;
  fhPHOSRecParticles   = new TH1D("PHOSRecParticles"   , "PHOSRecParticles" , 20 , 0., 20. ) ;
  fhPHOSPhotons        = new TH1I("PHOSPhotons"        , "PHOSPhotons"       , 20 , 0 , 20  ) ;
  fhPHOSInvariantMass  = new TH1D("PHOSInvariantMass"  , "PHOSInvariantMass" , 400, 0., 400.) ;
  fhPHOSDigitsEvent    = new TH1I("PHOSDigitsEvent"    , "PHOSDigitsEvent"   , 30 , 0 , 30  ) ;
  
  // create output container
  
  fOutputContainer = new TObjArray(8) ; 
  fOutputContainer->SetName(GetName()) ; 

  fOutputContainer->AddAt(fhPHOSPos,             0) ; 
  fOutputContainer->AddAt(fhPHOS,                1) ; 
  fOutputContainer->AddAt(fhPHOSEnergy,          2) ; 
  fOutputContainer->AddAt(fhPHOSDigits,          3) ; 
  fOutputContainer->AddAt(fhPHOSRecParticles,    4) ; 
  fOutputContainer->AddAt(fhPHOSPhotons,         5) ; 
  fOutputContainer->AddAt(fhPHOSInvariantMass,   6) ; 
  fOutputContainer->AddAt(fhPHOSDigitsEvent,     7) ; 
}

//______________________________________________________________________________
void AliPHOSQATask::Exec(Option_t *) 
{
  // Processing of one event
    
  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 
  
  //************************  PHOS *************************************
      
  Int_t       firstPhosCluster       = fESD->GetFirstPHOSCluster() ;
  const Int_t numberOfPhosClusters   = fESD->GetNumberOfPHOSClusters() ;
  
  TVector3 ** phosVector       = new TVector3*[numberOfPhosClusters] ;
  Float_t  * phosPhotonsEnergy = new Float_t[numberOfPhosClusters] ;
  Int_t      phosCluster ; 
  Int_t      numberOfPhotonsInPhos  = 0 ;
  Int_t      numberOfDigitsInPhos   = 0 ;
  
  // loop over the PHOS Cluster
  for(phosCluster = firstPhosCluster ; phosCluster < firstPhosCluster + numberOfPhosClusters ; phosCluster++) {
    AliESDCaloCluster * caloCluster = fESD->GetCaloCluster(phosCluster) ;
    if (caloCluster) {
      Float_t pos[3] ;
      caloCluster->GetGlobalPosition( pos ) ;
      fhPHOSEnergy->Fill( caloCluster->GetClusterEnergy() ) ;
      fhPHOSPos->Fill( pos[0], pos[1], pos[2] ) ;
      fhPHOSDigits->Fill(entry, caloCluster->GetNumberOfDigits() ) ;
      numberOfDigitsInPhos += caloCluster->GetNumberOfDigits() ;
      Float_t * pid = caloCluster->GetPid() ;
      if(pid[AliPID::kPhoton] > 0.9) {
	phosVector[numberOfPhotonsInPhos] = new TVector3(pos[0],pos[1],pos[2]) ;
	phosPhotonsEnergy[numberOfPhotonsInPhos]=caloCluster->GetClusterEnergy() ;
	numberOfPhotonsInPhos++;
      }
    }
  } //PHOS clusters
    
  fhPHOSRecParticles->Fill(numberOfPhosClusters);
  fhPHOSPhotons->Fill(numberOfPhotonsInPhos);
  fhPHOSDigitsEvent->Fill(numberOfDigitsInPhos);
  fhPHOS->Fill(entry, numberOfDigitsInPhos, numberOfPhosClusters, numberOfPhotonsInPhos) ; 

  // invariant Mass
  if (numberOfPhotonsInPhos > 1 ) {
    Int_t phosPhoton1, phosPhoton2 ; 
    for(phosPhoton1 = 0 ; phosPhoton1 < numberOfPhotonsInPhos ; phosPhoton1++) {
      for(phosPhoton2 = phosPhoton1 + 1 ; phosPhoton2 < numberOfPhotonsInPhos ; phosPhoton2++) {      
	Float_t tempMass = TMath::Sqrt( 2 * phosPhotonsEnergy[phosPhoton1] * phosPhotonsEnergy[phosPhoton2] *
					( 1 - TMath::Cos(phosVector[phosPhoton1]->Angle(*phosVector[phosPhoton2])) ) 
					);
	fhPHOSInvariantMass->Fill(tempMass*1000.);
      }
    }    
  }
  
  PostData(0, fOutputContainer);

  delete [] phosVector ; 
  delete [] phosPhotonsEnergy ; 
  
}

//______________________________________________________________________________
void AliPHOSQATask::Terminate(Option_t *)
{
  // Processing when the event loop is ended

  fOutputContainer = (TObjArray*)GetOutputData(0);  
  fhPHOSPos            = (TNtuple*)fOutputContainer->At(0);
  fhPHOS               = (TNtuple*)fOutputContainer->At(1);
  fhPHOSEnergy         = (TH1D*)fOutputContainer->At(2);
  fhPHOSDigits         = (TH1I*)fOutputContainer->At(3);
  fhPHOSRecParticles   = (TH1D*)fOutputContainer->At(4);
  fhPHOSPhotons        = (TH1I*)fOutputContainer->At(5);
  fhPHOSInvariantMass  = (TH1D*)fOutputContainer->At(6);
  fhPHOSDigitsEvent    = (TH1I*)fOutputContainer->At(7);

  printf("PHOSEnergy Mean         : %5.3f , RMS : %5.3f \n", fhPHOSEnergy->GetMean(),         fhPHOSEnergy->GetRMS()         ) ;
  printf("PHOSDigits Mean         : %5.3f , RMS : %5.3f \n", fhPHOSDigits->GetMean(),         fhPHOSDigits->GetRMS()         ) ;
  printf("PHOSRecParticles Mean   : %5.3f , RMS : %5.3f \n", fhPHOSRecParticles->GetMean(),   fhPHOSRecParticles->GetRMS()   ) ;
  printf("PHOSPhotons Mean        : %5.3f , RMS : %5.3f \n", fhPHOSPhotons->GetMean(),        fhPHOSPhotons->GetRMS()        ) ;
  printf("PHOSInvariantMass Mean  : %5.3f , RMS : %5.3f \n", fhPHOSInvariantMass->GetMean(),  fhPHOSInvariantMass->GetRMS()  ) ;
  printf("PHOSDigitsEvent Mean    : %5.3f , RMS : %5.3f \n", fhPHOSDigitsEvent->GetMean(),    fhPHOSDigitsEvent->GetRMS()    ) ;

  TCanvas  * cPHOS = new TCanvas("cPHOS", "PHOS ESD Test", 400, 10, 600, 700) ;
  cPHOS->Divide(3, 2);

  cPHOS->cd(1) ; 
  gPad->SetLogy();
  fhPHOSEnergy->SetAxisRange(0, 25.);
  fhPHOSEnergy->SetLineColor(2);
  fhPHOSEnergy->Draw();

  cPHOS->cd(2) ; 
  fhPHOSDigits->SetAxisRange(0,25.);
  fhPHOSDigits->SetLineColor(2);
  fhPHOSDigits->Draw();

  cPHOS->cd(3) ; 
  gPad->SetLogy();
  fhPHOSRecParticles->SetAxisRange(0, 25.);
  fhPHOSRecParticles->SetLineColor(2);
  fhPHOSRecParticles->Draw();

  cPHOS->cd(4) ; 
  gPad->SetLogy();
  fhPHOSPhotons->SetAxisRange(0,25.);
  fhPHOSPhotons->SetLineColor(2);
  fhPHOSPhotons->Draw();

  cPHOS->cd(5) ; 
  fhPHOSInvariantMass->SetLineColor(2);
  fhPHOSInvariantMass->Draw();
 
  cPHOS->cd(6) ; 
  gPad->SetLogy();
  fhPHOSDigitsEvent->SetAxisRange(0,40.);
  fhPHOSDigitsEvent->SetLineColor(2);
  fhPHOSDigitsEvent->Draw();
 
  cPHOS->Print("PHOS.eps");
 
  char line[1024] ; 
  sprintf(line, ".!tar -zcvf %s.tar.gz *.eps", GetName()) ; 
  gROOT->ProcessLine(line);
  sprintf(line, ".!rm -fR *.eps"); 
  gROOT->ProcessLine(line);
 
  AliInfo(Form("!!! All the eps files are in %s.tar.gz !!! \n", GetName())) ;
}
