/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

// Generator for particles according generic functions
// For high pt performance studies realistic distribution of the particles - local density in jets to be approximated 
//
//  
//  TF1 *   fF1Momentum;          // 1/momentum distribution function inGeV 
//  TF1 *   fFPhi;                // phi distribution function in rad     - if not set flat 0-2pi used
//  TF1 *   fFTheta;              // theta distribution function in rad   - if not set flat pi/4-3pi/4 used
//  TF3 *   fFPosition;           // position distribution function in cm - TO FIX if not set 0 used ()
//  TF1 *   fFPdg;                // pdg distribution function            - if not set flat jet pdg used 
//  We assume that the moment, postion and PDG code of particles are independent  
//  Only tracks/particle crossing the reference radius at given z range
//
// Origin: marian.ivanov@cern.ch


/*
  To test generator for particular setting run   
  AliGenPerformance::TestAliGenPerformance(Int_t nEvents, TF1 *f1pt, TF1 *fpdg){
  For distribution of 
    fF1Momentum=new TF1("f1pt","1-10*x",0,0.1);
    AliGenPerformance::TestAliGenPerformance(1000, fF1Momentum,0)
    pt distribution of charged primary particles described by powerlaw with slope -1.7
  
  gSystem->Load("libpythia6");
  gSystem->Load("libEGPythia6");
  gSystem->Load("liblhapdf");
  gSystem->Load("libAliPythia6");
  //
  AliGenPerformance::TestAliGenPerformance(5000,0);
  TFile * f = TFile::Open("testAliGenPerformance.root");
  testGener.Draw("pt>>his(100,10,100)","charge!=0&&abs(fKF)<2000","");
  f1pt = new TF1("f1pt","1/(0.01+x)",0,0.1);  
  his->Fit(f1);

*/

  


#include <TParticle.h>
#include <TF1.h>
#include <TF3.h>
#include <TDatabasePDG.h>

#include "AliRun.h"
#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliGenPerformance.h"
#include "AliGenEventHeader.h"
#include "TMCParticle.h"
#include <AliPythia.h>
#include <TTreeStream.h>

ClassImp(AliGenPerformance)

//-----------------------------------------------------------------------------
AliGenPerformance::AliGenPerformance():
  AliGenerator(),
  fNJets(1),                // mean number of jets per event
  fF1Momentum(0),           // momentum distribution function 
  fFPhi(0),                 // phi distribution function
  fFTheta(0),               // theta distribution function
  fFPosition(0),            // position distribution function 
  fFPdg(0),                 // pdg distribution function  
  fTestStream(0)            // test stream - used for tuning of parameters of generator
{
  //
  // Default constructor
  //
  SetNumberParticles(1);
}

AliGenPerformance::AliGenPerformance(const AliGenPerformance& func):
  AliGenerator(),
  fNJets(func.fNJets),                   //
  fF1Momentum(func.fF1Momentum),           // momentum distribution function 
  fFPhi(func.fFPhi),                     // phi distribution function
  fFTheta(func.fFTheta),                 // theta distribution function
  fFPosition(func.fFPosition),           // position distribution function 
  fFPdg(func.fFPdg)                     // pdg distribution function  
{
    // Copy constructor
    SetNumberParticles(1);
}

AliGenPerformance & AliGenPerformance::operator=(const AliGenPerformance& func)
{
    // Assigment operator
      if(&func == this) return *this;
      fNJets      = func.fNJets; 
      fF1Momentum = func.fF1Momentum; 
      fFPhi       = func.fFPhi;      
      fFTheta     = func.fFTheta;    
      fFPosition  = func.fFPosition; 
      fFPdg       = func.fFPdg;      
      return *this;
}


//-----------------------------------------------------------------------------
void AliGenPerformance::Generate()
{
  //
  // Generate one muon
  //
  Int_t naccepted =0;
  Int_t njets=gRandom->Poisson(fNJets);
  TDatabasePDG *databasePDG = TDatabasePDG::Instance();

  for (Int_t iparton=0; iparton<njets; iparton++){
    //
    //
    //
    Float_t mom[3];
    Float_t posf[3];
    Double_t pos[3];
    Int_t pdg;
    Double_t ptot, pt,  phi, theta; 
    //
    if (fF1Momentum){
      ptot     = 1./fF1Momentum->GetRandom();
    }else{
      ptot     = 0.001+fF1Momentum->GetRandom()*0.2;
      ptot/=ptot;
    }
    if (fFPhi){
      phi      = fFPhi->GetRandom();
    }else{
      phi      = gRandom->Rndm()*TMath::TwoPi();
    }
    if (fFTheta){
      theta    = fFTheta->GetRandom();
    }else{
      theta    =  (gRandom->Rndm()-0.5)*TMath::Pi()*0.5 +TMath::Pi()/2;
    }
    pt     = ptot*TMath::Sin(theta);
    mom[0] = pt*TMath::Cos(phi); 
    mom[1] = pt*TMath::Sin(phi); 
    mom[2] = ptot*TMath::Cos(theta);
    pos[0]=fOrigin[0]; pos[1]=fOrigin[1]; pos[2]=fOrigin[2];
    //
    if (fFPosition) fFPosition->GetRandom3(pos[0],pos[1],pos[2]); // ????
    //
    posf[0]=pos[0];
    posf[1]=pos[1];
    posf[2]=pos[2];
    if (fFPdg){
      pdg = TMath::Nint(fFPdg->GetRandom());
    }else{
      pdg = 1+TMath::Nint(gRandom->Rndm()*5.);
    }    
    Float_t polarization[3]= {0,0,0};
    Int_t nt;
    //
    AliPythia *py=AliPythia::Instance();
    py->Py1ent(-1, -pdg, ptot, theta, phi);
    py->Py1ent( 2,  pdg, ptot, theta, phi+TMath::Pi());
    py->Pyexec();
    TObjArray * array = py->GetPrimaries();
    Int_t nParticles=array->GetEntries();
    //array->Print();
    for (Int_t iparticle=2; iparticle<nParticles;  iparticle++){
      TMCParticle * mcParticle= (TMCParticle*)array->At(iparticle);
      Int_t flavour = mcParticle->GetKF();
      mom[0]=mcParticle->GetPx();
      mom[1]=mcParticle->GetPy();
      mom[2]=mcParticle->GetPz();
      //
      if (!fTestStream) PushTrack(fTrackIt,-1,flavour,mom, posf, polarization,0,kPPrimary,nt);
      if (fTestStream){
	TParticlePDG * pdgParticle=databasePDG->GetParticle(mcParticle->GetKF());
	if (pdgParticle){
	  Double_t charge=pdgParticle->Charge();
	  Double_t mass=pdgParticle->Mass();
	  Double_t  pt=TMath::Sqrt(mcParticle->GetPx()*mcParticle->GetPx()+mcParticle->GetPy()*mcParticle->GetPy());
	  (*fTestStream)<<"testGener"<<
	    "njets="<<njets<<
	    "ptot="<<ptot<<
	    "theta="<<theta<<
	    "phi="<<phi<<
	    "pdg="<<pdg<<
	    "charge="<<charge<<
	    "mass="<<mass<<
	    "pt="<<pt<<
	    "nParticles="<<nParticles<<
	    "ipart="<<iparticle<<
	    "mcParticle.="<<mcParticle<<
	    "\n";
	}
      }
      naccepted++;
    }
  }
  //  AliGenEventHeader* header = new AliGenEventHeader("THn");
  //gAlice->SetGenEventHeader(header);

  return;
}


//-----------------------------------------------------------------------------
void AliGenPerformance::Init()
{
  // 
  // Initialisation, check consistency of selected ranges
  //
  

  printf("************ AliGenPerformance ****************\n");
  printf("************************************************\n");
  if (!fF1Momentum){
    AliInfo("Momentum distribution function not specified");
  }
  if (!fFPhi){
    AliInfo("phi distribution function not specified");
  }
  if (!fFTheta){
    AliInfo("Theta distribution function not specified");
  }
  if (!fFPosition){
    AliInfo("Position distribution function not specified");
  }
  if (!fFPdg){
    AliInfo("PDG distribution function not specified");
  }

  return;
}

void AliGenPerformance::SetFunctions(TF1 * momentum, TF1 *fphi, TF1 *ftheta,TF3 * position, TF1* pdg){
  //
  // Set the function
  //
  fF1Momentum = momentum;
  fFPhi = fphi;
  fFTheta = ftheta;
  fFPosition = position;
  fFPdg = pdg;
}

void AliGenPerformance::TestAliGenPerformance(Int_t nEvents, TF1 *f1pt, TF1 *fpdg){
  //
  // test the genrator class - write particle  to the tree
  //
  AliGenPerformance *genPerformance= new AliGenPerformance;
  genPerformance->SetNJets(1);
  if (!f1pt) f1pt = new TF1("f1pt","1-10*x",0,0.1);
  if (!fpdg) fpdg = new TF1("f1pt","x",1,6);
  genPerformance->SetFunctions(f1pt,0,0,0,fpdg);
  TTreeSRedirector*pcstream = new TTreeSRedirector("testAliGenPerformance.root","recreate");
  genPerformance->SetStreamer(pcstream);
  for (Int_t i=0; i<nEvents;i++){
    genPerformance->Generate();
  }
  delete pcstream;
}

