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

//
// Wrapper for MEVSIM generator.
// It is using TMevSim to comunicate with fortarn code
// 
//	
// Sylwester Radomski <radomski@if.pw.edu.pl>
//

#include <Riostream.h>
#include <TClonesArray.h>
#include <TParticle.h>

#include "AliGenMevSim.h"
#include "AliMevSimConfig.h"
#include "AliMevSimParticle.h"
//#include "AliRun.h"
#include "TMevSim.h"

static TRandom * gAliRandom;

ClassImp(AliGenMevSim)

//____________________________________________________________________________
AliGenMevSim::AliGenMevSim() : AliGenerator(-1) 
{
  //
  // Default ctor
  //
  fConfig = new AliMevSimConfig();
  fMevSim = new TMevSim();
  gAliRandom = fRandom;
  
}
//____________________________________________________________________________
AliGenMevSim::AliGenMevSim(AliMevSimConfig *config): AliGenerator(-1) 
{
  //
  // Standard ctor
  //
  fConfig = config;
  fMevSim = new TMevSim(); 
  gAliRandom = fRandom;
}

//____________________________________________________________________________
AliGenMevSim::~AliGenMevSim() 
{
  //
  // Standard destructor
  //
  if (fMevSim) delete fMevSim;
}
//____________________________________________________________________________
void AliGenMevSim::SetConfig(AliMevSimConfig *config) 
{
  //
  // Sets the MevSim configuration
  //
  fConfig = config;
}

//____________________________________________________________________________
void AliGenMevSim::AddParticleType(AliMevSimParticle *type) 
{
  //
  // Add one particle type to MevSim
  //
  fMevSim->AddPartTypeParams((TMevSimPartTypeParams*)type);
}

//____________________________________________________________________________
void AliGenMevSim::Init() 
{
  //
  // Generator initialisation method
  //

  // fill data from AliMevSimConfig;

  TMevSim *mevsim = fMevSim;

  // geometry & momentum cut

  if (TestBit(kPtRange))  mevsim->SetPtCutRange(fPtMin, fPtMax);
  
  if (TestBit(kPhiRange)) // from radians to 0 - 360 deg
    mevsim->SetPhiCutRange( fPhiMin * 180 / TMath::Pi() , fPhiMax * 180 / TMath::Pi() );
  
  if (TestBit(kThetaRange)) // from theta to eta
  {
    mevsim->SetEtaCutRange( -TMath::Log( TMath::Tan(fThetaMax/2)) , -TMath::Log( TMath::Tan(fThetaMin/2)) );
  }  

  // mevsim specyfic parameters
  
  mevsim->SetModelType(fConfig->GetModelType());
  Int_t ctrl; Float_t psiRMean = 0; Float_t psiRStDev = 0;
  fConfig->GetRectPlane(ctrl,psiRMean,psiRStDev);
  mevsim->SetReacPlaneCntrl(ctrl);
  mevsim->SetPsiRParams(psiRMean, psiRStDev);
  Float_t mean; Float_t stDev;
  fConfig->GetMultFac(mean,stDev);
  mevsim->SetMultFacParams(mean,stDev);
  // mevsim initialisation

  mevsim->Initialize();
}

//____________________________________________________________________________
void AliGenMevSim::Generate() 
{
  //
  // Read the formatted output file and load the particles
  // Temporary solution
  //

  Int_t i;

  PDG_t pdg;
  Float_t orgin[3] = {0,0,0};
  Float_t polar[3] = {0,0,0};
  Float_t p[3] = {1,1,1};
  Float_t time = 0;
  
  const Int_t kParent = -1;
  Int_t id;

  // vertexing 

  VertexInternal();

  orgin[0] = fVertex[0];
  orgin[1] = fVertex[1];
  orgin[2] = fVertex[2];

  cout << "Vertex ";
  for (i =0; i<3; i++)
    cout << orgin[i] << "\t";
  cout << endl;

  Int_t nParticles = 0;

  TClonesArray *particles = new TClonesArray("TParticle");
  TParticle *particle;

  fMevSim->GenerateEvent();
  fNpart= fMevSim->ImportParticles(particles,"");

  cout << "Found " << fNpart << " particles ..." << endl;
  nParticles = fNpart;

  for (i=0; i<nParticles; i++) {
    
    particle = (TParticle*) (*particles)[i];

    pdg = (PDG_t) particle->GetPdgCode();
    p[0] = particle->Px();
    p[1] = particle->Py();
    p[2] = particle->Pz();
    
    PushTrack(fTrackIt, kParent, pdg, p, orgin, polar, time, kPPrimary, id);

  }  
 
  particles->Clear();
  if (particles) delete particles;
}

//____________________________________________________________________________
#ifndef WIN32
# define ran ran_
# define type_of_call
#else
# define ran RAN
# define type_of_call _stdcall
#endif

extern "C" Float_t type_of_call ran(Int_t &)
{
  //
  //  Replacement for package random number generator
  //
  return gAliRandom->Rndm(); 
}

