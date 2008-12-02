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

/* $Id: */

//_________________________________________________________________________
// A basic analysis task to analyse photon detected by PHOS
// A basic analysis task to analyse photon detected by PHOS
// A basic analysis task to analyse photon detected by PHOS
// Adapted for AliAnalysisTaskSE and AOD production 
// by Gustavo Conesa
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////
//Root system
#include <TCanvas.h>
#include <TH1.h>
#include <TROOT.h>
#include <TLegend.h> 
#include <TNtuple.h>
#include <TVector3.h> 
#include <TTree.h>

//Analysis system
#include "AliAnalysisTaskPHOSExample.h" 
#include "AliAnalysisManager.h"
#include "AliESDEvent.h" 
#include "AliESDCaloCluster.h" 
#include "AliAODEvent.h"
#include "AliAODPhoton.h"
#include "AliLog.h"
#include "AliESDVertex.h"

//______________________________________________________________________________
AliAnalysisTaskPHOSExample::AliAnalysisTaskPHOSExample() : 
  fDebug(0),  
  fAODPhotons(0x0), 
  fPhotonsInPhos(0),
  fPhotonId(1.0),
  fOutputList(0x0), 
  fhPHOSPos(0),
  fhPHOS(0),
  fhPHOSEnergy(0),
  fhPHOSDigits(0),
  fhPHOSRecParticles(0),
  fhPHOSPhotons(0),
  fhPHOSInvariantMass(0),
  fhPHOSDigitsEvent(0)
{
  //Default constructor
}
//______________________________________________________________________________
AliAnalysisTaskPHOSExample::AliAnalysisTaskPHOSExample(const char *name) : 
  AliAnalysisTaskSE(name),  
  fDebug(0),
  fAODPhotons(0x0), 
  fPhotonsInPhos(0),
  fPhotonId(1.0),
  fOutputList(0x0), 
  fhPHOSPos(0),
  fhPHOS(0),
  fhPHOSEnergy(0),
  fhPHOSDigits(0),
  fhPHOSRecParticles(0),
  fhPHOSPhotons(0),
  fhPHOSInvariantMass(0),
  fhPHOSDigitsEvent(0)
{
  // Constructor.

  // Output slots 
  DefineOutput(1,  TList::Class()) ; 
}

//____________________________________________________________________________
AliAnalysisTaskPHOSExample::AliAnalysisTaskPHOSExample(const AliAnalysisTaskPHOSExample& ap) :
  AliAnalysisTaskSE(ap.GetName()),  
  fDebug(ap.fDebug),
  fAODPhotons(ap.fAODPhotons), 
  fPhotonsInPhos(ap.fPhotonsInPhos),
  fPhotonId(ap.fPhotonId),
  fOutputList(ap.fOutputList), 
  fhPHOSPos(ap.fhPHOSPos),
  fhPHOS(ap.fhPHOS),
  fhPHOSEnergy(ap.fhPHOSEnergy),
  fhPHOSDigits(ap.fhPHOSDigits),
  fhPHOSRecParticles(ap.fhPHOSRecParticles),
  fhPHOSPhotons(ap.fhPHOSPhotons),
  fhPHOSInvariantMass(ap.fhPHOSInvariantMass),
  fhPHOSDigitsEvent(ap.fhPHOSDigitsEvent)
{ 
  // cpy ctor
}

//_____________________________________________________________________________
AliAnalysisTaskPHOSExample& AliAnalysisTaskPHOSExample::operator = (const AliAnalysisTaskPHOSExample& ap)
{
// assignment operator

  this->~AliAnalysisTaskPHOSExample();
  new(this) AliAnalysisTaskPHOSExample(ap);
  
  fDebug = ap.fDebug;
  fAODPhotons = ap.fAODPhotons; 
  fPhotonsInPhos = ap.fPhotonsInPhos;
  fPhotonId = ap.fPhotonId;
  fOutputList = ap.fOutputList;
  fhPHOSPos = ap.fhPHOSPos;
  fhPHOS = ap.fhPHOS;
  fhPHOSEnergy = ap.fhPHOSEnergy;
  fhPHOSDigits = ap.fhPHOSDigits;
  fhPHOSRecParticles = ap.fhPHOSRecParticles;
  fhPHOSPhotons = ap.fhPHOSPhotons;
  fhPHOSInvariantMass = ap.fhPHOSInvariantMass;
  fhPHOSDigitsEvent = ap.fhPHOSDigitsEvent;
  
  return *this;
}

//______________________________________________________________________________
AliAnalysisTaskPHOSExample::~AliAnalysisTaskPHOSExample()
{
  // dtor
  if(fOutputList) {
    fOutputList->Clear() ; 
    delete fOutputList ;
  }
}


//________________________________________________________________________
void AliAnalysisTaskPHOSExample::UserCreateOutputObjects()
{  
  // Create the outputs containers
   //AODs
  fAODPhotons = new TClonesArray("AliAODPhoton", 0);
  fAODPhotons->SetName("Photons");
  AddAODBranch("TClonesArray", &fAODPhotons);

  OpenFile(1) ; 

  fhPHOSPos            = new TNtuple("PHOSPos"         , "Position in PHOS"  , "x:y:z");
  fhPHOS               = new TNtuple("PHOS"            , "PHOS"  , "event:digits:clusters:photons");
  fhPHOSEnergy         = new TH1D("PHOSEnergy"         , "PHOSEnergy"        , 100, 0., 100. ) ;
  fhPHOSDigits         = new TH1I("PHOSDigitsCluster"  , "PHOSDigits"        , 20 , 0 , 20  ) ;
  fhPHOSRecParticles   = new TH1D("PHOSRecParticles"   , "PHOSRecParticles" , 20 , 0., 20. ) ;
  fhPHOSPhotons        = new TH1I("PHOSPhotons"        , "PHOSPhotons"       , 20 , 0 , 20  ) ;
  fhPHOSInvariantMass  = new TH1D("PHOSInvariantMass"  , "PHOSInvariantMass" , 400, 0., 400.) ;
  fhPHOSDigitsEvent    = new TH1I("PHOSDigitsEvent"    , "PHOSDigitsEvent"   , 30 , 0 , 30  ) ;
  
  // create output container
  
  fOutputList = new TList() ; 
  fOutputList->SetName(GetName()) ; 

  fOutputList->AddAt(fhPHOSPos,             0) ; 
  fOutputList->AddAt(fhPHOS,                1) ; 
  fOutputList->AddAt(fhPHOSEnergy,          2) ; 
  fOutputList->AddAt(fhPHOSDigits,          3) ; 
  fOutputList->AddAt(fhPHOSRecParticles,    4) ; 
  fOutputList->AddAt(fhPHOSPhotons,         5) ; 
  fOutputList->AddAt(fhPHOSInvariantMass,   6) ; 
  fOutputList->AddAt(fhPHOSDigitsEvent,     7) ; 


}

//______________________________________________________________________________
void AliAnalysisTaskPHOSExample::UserExec(Option_t *) 
{
  // Processing of one event
  
  // if ( !((Entry()-1)%100) ) 
  AliInfo(Form(" Processing event # %lld",  Entry())) ; 
  AliESDEvent* esd = (AliESDEvent*)InputEvent();
  
  //************************  PHOS *************************************
      TRefArray * caloClustersArr  = new TRefArray();  
  esd->GetPHOSClusters(caloClustersArr);
  
  const Int_t kNumberOfPhosClusters   = caloClustersArr->GetEntries() ;  
  
  TVector3 ** phosVector       = new TVector3*[kNumberOfPhosClusters] ;
  Float_t  * phosPhotonsEnergy = new Float_t[kNumberOfPhosClusters] ;
  Int_t      phosCluster ; 
  Int_t      numberOfDigitsInPhos   = 0 ;
  Double_t v[3] ; //vertex ;
  esd->GetVertex()->GetXYZ(v) ;

  fPhotonsInPhos  = 0 ;
  // loop over the PHOS Cluster
  for(phosCluster = 0 ; phosCluster < kNumberOfPhosClusters ; phosCluster++) {
    AliESDCaloCluster * caloCluster = (AliESDCaloCluster *) caloClustersArr->At(phosCluster) ;
  
    //AliESDCaloCluster * caloCluster = fESD->GetCaloCluster(phosCluster) ;
    if (caloCluster) {
      Float_t pos[3] ;
      caloCluster->GetPosition( pos ) ;
      fhPHOSEnergy->Fill( caloCluster->E() ) ;
      fhPHOSPos->Fill( pos[0], pos[1], pos[2] ) ;
      fhPHOSDigits->Fill(Entry(), caloCluster->GetNCells() ) ;
      numberOfDigitsInPhos += caloCluster->GetNCells() ;
      Double_t * pid = caloCluster->GetPid() ;
      if(pid[AliPID::kPhoton] > GetPhotonId() ) {
	phosVector[fPhotonsInPhos] = new TVector3(pos[0],pos[1],pos[2]) ;
	phosPhotonsEnergy[fPhotonsInPhos]=caloCluster->E() ;
	TLorentzVector momentum ;
	caloCluster->GetMomentum(momentum, v);   
	new ((*fAODPhotons)[fPhotonsInPhos++]) AliAODPhoton (momentum);
      }
    }
  } //PHOS clusters
    
  fhPHOSRecParticles->Fill(kNumberOfPhosClusters);
  fhPHOSPhotons->Fill(fPhotonsInPhos);
  fhPHOSDigitsEvent->Fill(numberOfDigitsInPhos);
  fhPHOS->Fill(Entry(), numberOfDigitsInPhos, kNumberOfPhosClusters, fPhotonsInPhos) ; 

  // invariant Mass
  if (fPhotonsInPhos > 1 ) {
    Int_t phosPhoton1, phosPhoton2 ; 
    for(phosPhoton1 = 0 ; phosPhoton1 < fPhotonsInPhos ; phosPhoton1++) {
      for(phosPhoton2 = phosPhoton1 + 1 ; phosPhoton2 < fPhotonsInPhos ; phosPhoton2++) {      
	Float_t tempMass = TMath::Sqrt( 2 * phosPhotonsEnergy[phosPhoton1] * phosPhotonsEnergy[phosPhoton2] *
					( 1 - TMath::Cos(phosVector[phosPhoton1]->Angle(*phosVector[phosPhoton2])) ) 
					);
	fhPHOSInvariantMass->Fill(tempMass*1000.);
      }
    }    
  }
    
  PostData(1, fOutputList);

  delete [] phosVector ; 
  delete [] phosPhotonsEnergy ; 
  
}


//______________________________________________________________________________
void AliAnalysisTaskPHOSExample::Init()
{
  // Intialisation of parameters
  AliInfo("Doing initialisation") ; 
  fPhotonId = 0.9 ; 
}

//______________________________________________________________________________
void AliAnalysisTaskPHOSExample::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  
  Bool_t problem = kFALSE ; 
  AliInfo(Form(" *** %s Report:", GetName())) ; 
  printf("        PHOSEnergy Mean         : %5.3f , RMS : %5.3f \n", fhPHOSEnergy->GetMean(),         fhPHOSEnergy->GetRMS()         ) ;
  printf("        PHOSDigits Mean         : %5.3f , RMS : %5.3f \n", fhPHOSDigits->GetMean(),         fhPHOSDigits->GetRMS()         ) ;
  printf("        PHOSRecParticles Mean   : %5.3f , RMS : %5.3f \n", fhPHOSRecParticles->GetMean(),   fhPHOSRecParticles->GetRMS()   ) ;
  printf("        PHOSPhotons Mean        : %5.3f , RMS : %5.3f \n", fhPHOSPhotons->GetMean(),        fhPHOSPhotons->GetRMS()        ) ;
  printf("        PHOSInvariantMass Mean  : %5.3f , RMS : %5.3f \n", fhPHOSInvariantMass->GetMean(),  fhPHOSInvariantMass->GetRMS()  ) ;
  printf("        PHOSDigitsEvent Mean    : %5.3f , RMS : %5.3f \n", fhPHOSDigitsEvent->GetMean(),    fhPHOSDigitsEvent->GetRMS()    ) ;

  TCanvas  * cPHOS = new TCanvas("cPHOS", "PHOS ESD Test", 400, 10, 600, 700) ;
  cPHOS->Divide(3, 2);

  cPHOS->cd(1) ; 
  if ( fhPHOSEnergy->GetMaximum() > 0. ) 
    gPad->SetLogy();
  fhPHOSEnergy->SetAxisRange(0, 25.);
  fhPHOSEnergy->SetLineColor(2);
  fhPHOSEnergy->Draw();

  cPHOS->cd(2) ; 
  fhPHOSDigits->SetAxisRange(0,25.);
  fhPHOSDigits->SetLineColor(2);
  fhPHOSDigits->Draw();

  cPHOS->cd(3) ; 
  if ( fhPHOSRecParticles->GetMaximum() > 0. ) 
    gPad->SetLogy();
  fhPHOSRecParticles->SetAxisRange(0, 25.);
  fhPHOSRecParticles->SetLineColor(2);
  fhPHOSRecParticles->Draw();

  cPHOS->cd(4) ; 
  if ( fhPHOSPhotons->GetMaximum() > 0. ) 
    gPad->SetLogy();
  fhPHOSPhotons->SetAxisRange(0,25.);
  fhPHOSPhotons->SetLineColor(2);
  fhPHOSPhotons->Draw();

  cPHOS->cd(5) ; 
  fhPHOSInvariantMass->SetLineColor(2);
  fhPHOSInvariantMass->Draw();
 
  cPHOS->cd(6) ; 
  if ( fhPHOSDigitsEvent->GetMaximum() > 0. ) 
    gPad->SetLogy();
  fhPHOSDigitsEvent->SetAxisRange(0,40.);
  fhPHOSDigitsEvent->SetLineColor(2);
  fhPHOSDigitsEvent->Draw();
 
  cPHOS->Print("PHOS.eps");
 
  char line[1024] ; 
  sprintf(line, ".!tar -zcf %s.tar.gz *.eps", GetName()) ; 
  gROOT->ProcessLine(line);
  sprintf(line, ".!rm -fR *.eps"); 
  gROOT->ProcessLine(line);
 
  AliInfo(Form("!!! All the eps files are in %s.tar.gz !!!", GetName())) ;

  char * report ; 
  if(problem)
    report="Problems found, please check!!!";  
  else 
    report="OK";

  AliInfo(Form("*** %s Summary Report: %s \n",GetName(), report)) ; 
}
