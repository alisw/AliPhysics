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

//___________________________________________________
////////////////////////////////////////////////////////////////////////////////
// 
// class AliHBTPositionRandomizer
//
// These class randomizes particle vertex positions
// Piotr.Skowronski@cern.ch
//
////////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TRandom.h>

#include "AliAOD.h"
#include "AliHBTPositionRandomizer.h"
#include "AliLog.h"
#include "AliVAODParticle.h"


ClassImp(AliHBTPositionRandomizer)

const Int_t AliHBTPositionRandomizer::fgkNumberOfPids = 10;

/*********************************************************************/

AliHBTPositionRandomizer::AliHBTPositionRandomizer():
 fReader(0x0),
 fDefaultRandomizer(0x0),
 fRandomizers(0x0),
 fNPid(0),
 fPids(0x0),
 fAddToExistingPos(kFALSE),
 fOnlyParticlesFromVertex(kFALSE),
 fRandomizeTracks(kFALSE),
 fVX(0.0),
 fVY(0.0),
 fVZ(0.0)
{
//constructor
}
/*********************************************************************/

AliHBTPositionRandomizer::AliHBTPositionRandomizer(AliReader* reader):
 fReader(reader),
 fRandomizers(new TObjArray(fgkNumberOfPids)),
 fNPid(1),
 fPids(new Int_t[fgkNumberOfPids]),
 fAddToExistingPos(kFALSE),
 fOnlyParticlesFromVertex(kFALSE),
 fRandomizeTracks(kFALSE),
 fVX(0.0),
 fVY(0.0),
 fVZ(0.0)
{
//constructor
 fRandomizers->AddAt(new AliHBTRndmGaussBall(8.0,0.0,0.0),0);
} 
/*********************************************************************/

AliHBTPositionRandomizer::AliHBTPositionRandomizer(const AliHBTPositionRandomizer& in):
 AliReader(in),
 fReader(),
 fRandomizers(0x0),
 fNPid(0),
 fPids(0x0),
 fAddToExistingPos(kFALSE),
 fOnlyParticlesFromVertex(kFALSE),
 fRandomizeTracks(kFALSE),
 fVX(0.0),
 fVY(0.0),
 fVZ(0.0)
{
  //cpy constructor
  in.Copy(*this);
}
/*********************************************************************/
AliHBTPositionRandomizer::~AliHBTPositionRandomizer()
{
  //dtor
  delete fReader;
  delete fRandomizers;
  delete [] fPids;
}
/*********************************************************************/
AliHBTPositionRandomizer& AliHBTPositionRandomizer::operator=(const AliHBTPositionRandomizer& in)
{
  //assigment operator
  in.Copy(*this);
  return *this;
}
/*********************************************************************/

AliAOD* AliHBTPositionRandomizer::GetEventSim() const
{
 // gets from fReader and randomizes current particle event
 if (fReader == 0x0) 
  {
    Error("GetEventSim","Reader is null");
    return 0x0;
  } 
 AliAOD *e =  fReader->GetEventSim();
 if (e) 
   if (e->IsRandomized() == kFALSE) 
     Randomize(e);
 return e;
}
/*********************************************************************/

AliAOD* AliHBTPositionRandomizer::GetEventRec() const
{
 // gets from fReader and randomizes current track event
 if (fReader == 0x0) 
  {
    Error("GetEventRec","Reader is null");
    return 0x0;
  }  
 AliAOD *e =  fReader->GetEventRec();
 if (fRandomizeTracks && e) if (e->IsRandomized() == kFALSE) Randomize(e);
 return e;
}
/*********************************************************************/

AliAOD* AliHBTPositionRandomizer::GetEventSim(Int_t n)
{
//returns event n
 if (fReader == 0x0) return 0x0;
 AliAOD *e =  fReader->GetEventSim(n);
 if (e->IsRandomized() == kFALSE) Randomize(e);
 return e;
}

/*********************************************************************/
void AliHBTPositionRandomizer::Randomize(AliAOD* event) const
{
// randomizes postions of all particles in the event
  static const Double_t kfmtocm = 1.e-13;
  AliDebug(5," ");
  if (event == 0x0) return;

  for (Int_t i = 0; i < event->GetNumberOfParticles(); i++)
   {
     AliVAODParticle* p = event->GetParticle(i);
     Double_t x,y,z,t=0.0;
     AliHBTRndm* r = GetRandomizer(p->GetPdgCode());
     r->Randomize(x,y,z,t,p);
     
     Double_t nx = x*kfmtocm;
     Double_t ny = y*kfmtocm;
     Double_t nz = z*kfmtocm;
     Double_t nt = t*kfmtocm;
     
     if (fAddToExistingPos)
      {
       nx += p->Vx();
       ny += p->Vy();
       nz += p->Vz();
       nt += p->T();
      }
     p->SetProductionVertex(nx,ny,nz,nt); 
   }
  event->SetRandomized();
}
/*********************************************************************/

AliHBTRndm* AliHBTPositionRandomizer::GetRandomizer(Int_t pdg) const
{
  //returns randomizer for a given pdg 
  Int_t idx = GetRandomizerIndex(pdg);//in most of cases 
  if (idx < 0) idx = 0;//if not found return a default one
  return (AliHBTRndm*)fRandomizers->At(idx);
}
/*********************************************************************/
Int_t AliHBTPositionRandomizer::GetRandomizerIndex(Int_t pdg) const
{
  //returns randomizer index for a given pdg 

  if (pdg == 0) return 0;
  
  for (Int_t i=1; i < fNPid; i++)
   {
     if (fPids[i] == pdg) 
      return i;
   }
   
  return -1;
}
/*********************************************************************/

void AliHBTPositionRandomizer::SetRandomizer(Int_t pid, AliHBTRndm* rndm)
{
 //sets the randomizer for a given particle type
  if (rndm == 0x0)
   {
     Error("SetRandomizer","Randomizer is null");
     return;
   }
   
  Int_t idx = GetRandomizerIndex(pid);
  if (idx >= 0) 
   {
     delete fRandomizers->At(idx);
     fRandomizers->AddAt(rndm,idx);
   }  
  
  if (fNPid == fgkNumberOfPids)
   {
     Error("SetRandomizer","There is no more space in the array");
     return;
   }

  fPids[fNPid] = pid;
  fRandomizers->AddAt(rndm,fNPid);
  fNPid++;
}
/*********************************************************************/

void AliHBTPositionRandomizer::SetGaussianBall(Int_t pid, Double_t r, Double_t meantime, Double_t sigmatime)
{
 //Sets Gaussian Ball Model
  SetGaussianBall(pid,r,r,r,meantime,sigmatime);
}
/*********************************************************************/

void AliHBTPositionRandomizer::SetGaussianBall(Int_t pid, Double_t rx, Double_t ry, Double_t rz, Double_t meantime, Double_t sigmatime)
{
 //Sets Gaussian Ball Model
  AliHBTRndm* rndm = new AliHBTRndmGaussBall(rx,ry,rz,meantime,sigmatime);
  SetRandomizer(pid,rndm);
}
/*********************************************************************/

void AliHBTPositionRandomizer::SetCyllinderSurface(Int_t pid, Double_t r, Double_t l)
{
 //Sets Cylinder Surface Model
  AliHBTRndm* rndm = new  AliHBTRndmCyllSurf(r,l);
  SetRandomizer(pid,rndm);
}
/*********************************************************************/

void AliHBTPositionRandomizer::SetEventVertex(Double_t x, Double_t y,Double_t z)
{
//sets event vertex position
  fVX = x;
  fVY = y;
  fVZ = z;
}


void AliHBTPositionRandomizer::SetEllipse(Int_t pid, Double_t rx, Double_t ryz)
{
//sets the ellipse randomization for the given pid
  AliHBTRndm* rndm = new AliHBTRndmEllipse(rx,ryz);
  SetRandomizer(pid,rndm);   
}

/*********************************************************************/
//_____________________________________________________________________
///////////////////////////////////////////////////////////////////////
//                                                                   //
//  class AliHBTRndmGaussBall                                        //
//                                                                   //
///////////////////////////////////////////////////////////////////////

AliHBTRndmGaussBall::AliHBTRndmGaussBall():
 fRx(0.0),
 fRy(0.0),
 fRz(0.0),
 fTmean(0.0),
 fTsigma(0.0)
{
  //constructor
}
/*********************************************************************/

AliHBTRndmGaussBall::AliHBTRndmGaussBall(Float_t r, Double_t meantime, Double_t sigmatime):
 fRx(r),
 fRy(r),
 fRz(r),
 fTmean(meantime),
 fTsigma(sigmatime)
{
  //constructor
}
/*********************************************************************/

AliHBTRndmGaussBall::AliHBTRndmGaussBall(Float_t rx, Float_t ry, Float_t rz, Double_t meantime, Double_t sigmatime):
 fRx(rx),
 fRy(ry),
 fRz(rz),
 fTmean(meantime),
 fTsigma(sigmatime)
{
  //constructor
}
/*********************************************************************/


AliHBTRndmEllipse::AliHBTRndmEllipse(Float_t rmin, Float_t rmax):
 fRmin(rmin),
 fRmax(rmax)
{
     //constructor
}

/*********************************************************************/

void AliHBTRndmGaussBall::Randomize(Double_t& x,Double_t& y,Double_t&z,Double_t&t, AliVAODParticle*/*particle*/) const
{
//randomizez gauss for each coordinate separately
  x = gRandom->Gaus(0.0,fRx);
  y = gRandom->Gaus(0.0,fRy);
  z = gRandom->Gaus(0.0,fRz);
  
  if (fTsigma == 0.0)
   {
     t = fTmean;
     return;
   }
  
  t = gRandom->Gaus(fTmean,fTsigma);
    
}
/*********************************************************************/
//_____________________________________________________________________
///////////////////////////////////////////////////////////////////////
//                                                                   //
//  class AliHBTRndmGaussBall                                        //
//                                                                   //
///////////////////////////////////////////////////////////////////////

void AliHBTRndmCyllSurf::Randomize(Double_t& x,Double_t& y,Double_t&z,Double_t&/*t*/, AliVAODParticle* particle) const
{
//Randomizes x,y,z
   Double_t r = fR + gRandom->Gaus(0.0, 1.0);
   Double_t sf = r/particle->Pt();//scaling factor for position transformation ->
                             //we move direction of string momentum but legth defined by r
   x = sf*particle->Px();
   y = sf*particle->Py();
   z = gRandom->Uniform(-fL,fL);
}

/*********************************************************************/
/*********************************************************************/

void AliHBTRndmEllipse::Randomize(Double_t& x, Double_t& y, Double_t& z,Double_t&/*t*/, AliVAODParticle*p) const
{
    // p=0; //workaround - fix this damn little thingy
   double R;
     double phi=p->Phi();
     
     R=fRmin+(fRmax-fRmin)*TMath::Sin(phi);
     x=R*TMath::Sin(phi);
     y=R*TMath::Cos(phi);
     z=z;
}
