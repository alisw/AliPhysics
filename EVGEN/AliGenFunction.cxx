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
//  
//  TF1 *   fFMomentum;           // momentum distribution function inGeV
//  TF1 *   fFPhi;                // phi distribution function in rad
//  TF1 *   fFTheta;              // theta distribution function in rad
//  TF3 *   fFPosition;           // position distribution function in cm
//  TF1 *   fFPdg;                // pdg distribution function  
//  We assume that the moment, postion and PDG code of particles are independent  
//  Only tracks/particle crossing the reference radius at given z range
//
// Origin: marian.ivanov@cern.ch


/*
  Example generic function for cosmic generation:
 
//
  	AliGenFunction *generCosmic = new AliGenFunction;
	generCosmic->SetNumberParticles(100);
	TF1 *fmom   = new TF1("fmom","TMath::Exp(-x/8)",0,30);
	//
	TF1 *fphi   = new TF1("fphi","TMath::Gaus(x,-TMath::Pi()/2,0.3)",-3.14,3.14);
	TF1 *ftheta = new TF1("ftheta","TMath::Gaus(x,TMath::Pi()/2,0.3)",-3.14,3.14);
	TF3 *fpos = new TF3("fpos","1+(x+y+z)*0",-2000,2000,700,701,-2000,2000);
	TF1 * fpdg= new TF1("fpdg","1*(abs(x-13)<0.1)+1*(abs(x+13)<0.1)",-300,300);
	fpdg->SetNpx(10000);  // neccessary - number of bins for generation
	generCosmic->SetVertexSmear(kPerEvent);
	generCosmic->SetFunctions(fmom,fphi,ftheta,fpos,fpdg);
	generCosmic->SetCylinder(375,-375,375);
	generCosmic->SetBkG(0.2);
	gener->AddGenerator(generCosmic,"generCosmic",1);


  
*/




#include <TParticle.h>
#include <TF1.h>
#include <TF3.h>
#include <TDatabasePDG.h>

#include "AliRun.h"
#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliGenFunction.h"
#include "AliGenEventHeader.h"

ClassImp(AliGenFunction)

//-----------------------------------------------------------------------------
AliGenFunction::AliGenFunction():
  AliGenerator(),
  fBkG(0.),  
  fFMomentum(0),           // momentum distribution function 
  fFPhi(0),                // phi distribution function
  fFTheta(0),              // theta distribution function
  fFPosition(0),           // position distribution function 
  fFPdg(0),                // pdg distribution function  
  fRefRadius(0),           // reference radius to be crossed
  fZmin(0),                // minimal z at reference radius
  fZmax(0),                // z at reference radius
  fMaxTrial(10000)         // maximal number of attempts
{
  //
  // Default constructor
  //
  SetNumberParticles(1);
}

AliGenFunction::AliGenFunction(const AliGenFunction& func):
  AliGenerator(),
  fBkG(func.fBkG),  
  fFMomentum(func.fFMomentum),           // momentum distribution function 
  fFPhi(func.fFPhi),                     // phi distribution function
  fFTheta(func.fFTheta),                 // theta distribution function
  fFPosition(func.fFPosition),           // position distribution function 
  fFPdg(func.fFPdg),                     // pdg distribution function  
  fRefRadius(func.fRefRadius),           // reference radius to be crossed
  fZmin(func.fZmin),                     // minimal z at reference radius
  fZmax(func.fZmax),                     // z at reference radius
  fMaxTrial(10000)                       // maximal number of attempts
{
    // Copy constructor
    SetNumberParticles(1);
}

AliGenFunction & AliGenFunction::operator=(const AliGenFunction& func)
{
    // Assigment operator
      if(&func == this) return *this;
      fBkG       = func.fBkG;
      fFMomentum = func.fFMomentum; 
      fFPhi      = func.fFPhi;      
      fFTheta    = func.fFTheta;    
      fFPosition = func.fFPosition; 
      fFPdg      = func.fFPdg;      
      fRefRadius = func.fRefRadius;      
      fZmin      = func.fZmin;                
      fZmax      = func.fZmax;                
      fMaxTrial  = func.fMaxTrial;
      return *this;
}


//-----------------------------------------------------------------------------
void AliGenFunction::Generate()
{
  //
  // Generate one muon
  //
  Int_t naccepted =0;
  
  for (Int_t ipart=0; ipart<fMaxTrial && naccepted<fNpart; ipart++){
    //
    //
    //
    Float_t mom[3];
    Float_t posf[3];
    Double_t pos[3];
    Int_t pdg;
    Double_t ptot, pt,  phi, theta; 
    //
    ptot     = fFMomentum->GetRandom();
    phi      = fFPhi->GetRandom();
    theta    = fFTheta->GetRandom();
    pt     = ptot*TMath::Sin(theta);
    mom[0] = pt*TMath::Cos(phi); 
    mom[1] = pt*TMath::Sin(phi); 
    mom[2] = ptot*TMath::Cos(theta);

    //
    fFPosition->GetRandom3(pos[0],pos[1],pos[2]);
    posf[0]=pos[0];
    posf[1]=pos[1];
    posf[2]=pos[2];
    pdg = TMath::Nint(fFPdg->GetRandom());
    if (!IntersectCylinder(fRefRadius,fZmin, fZmax, pdg, posf, mom)) continue;
    //
    //
    Float_t polarization[3]= {0,0,0};
    Int_t nt;
    PushTrack(fTrackIt,-1,pdg,mom, posf, polarization,0,kPPrimary,nt);
    naccepted++;
  }

  AliGenEventHeader* header = new AliGenEventHeader("THn");
  gAlice->SetGenEventHeader(header);
  return;
}
//-----------------------------------------------------------------------------
void AliGenFunction::Init()
{
  // 
  // Initialisation, check consistency of selected ranges
  //
  printf("************ AliGenFunction ****************\n");
  printf("************************************************\n");
  if (!fFMomentum){
    AliFatal("Momentum distribution function not specified");
  }
  if (!fFPosition){
    AliFatal("Position distribution function not specified");
  }
  if (!fFPdg){
    AliFatal("PDG distribution function not specified");
  }
  if (fZmin==fZmax){
    AliFatal("Z range not specified");
  }

  return;
}

void AliGenFunction::SetFunctions(TF1 * momentum, TF1 *fphi, TF1 *ftheta,TF3 * position, TF1* pdg){
  //
  // Set the function
  //
  fFMomentum = momentum;
  fFPhi = fphi;
  fFTheta = ftheta;
  fFPosition = position;
  fFPdg = pdg;
}

void AliGenFunction::SetCylinder(Double_t refR, Double_t zmin, Double_t zmax){
  //
  // Set the cylinder geometry
  //
  fRefRadius = refR;          // reference radius to be crossed
  fZmin = zmin;               // minimal z at reference radius
  fZmax = zmax;               // maximal z at reference radius

}


//-----------------------------------------------------------------------------
Bool_t AliGenFunction::IntersectCylinder(Float_t r,Float_t zmin, Float_t zmax,Int_t pdg,
					 Float_t o[3],Float_t p[3]) const
{
  //
  // Intersection between muon and cylinder [-z,+z] with radius r
  //
  Double_t mass=0;
  if (TDatabasePDG::Instance()->GetParticle(pdg)){
    mass = TDatabasePDG::Instance()->GetParticle(13)->Mass();
  }

  Float_t en = TMath::Sqrt(mass*mass+p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  TParticle part(pdg,0,0,0,0,0,p[0],p[1],p[2],en,o[0],o[1],o[2],0);
  AliESDtrack track(&part);
  Double_t pos[3]={0.,0.,0.},sigma[3]={0.,0.,0.};
  AliESDVertex origin(pos,sigma);

  track.RelateToVertex(&origin,fBkG,10000.);

  Float_t d0z0[2],covd0z0[3];
  track.GetImpactParameters(d0z0,covd0z0);

  // check rphi 
  if(TMath::Abs(d0z0[0])>r) return kFALSE;
  // check z
  if(d0z0[1]>zmax) return kFALSE;
  if(d0z0[1]<zmin) return kFALSE;
  //
  return kTRUE;
}
