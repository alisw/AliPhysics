#include "AliHBTPositionRandomizer.h"
#include <TRandom.h>
#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTParticle.h"


ClassImp(AliHBTPositionRandomizer)

/*********************************************************************/

AliHBTPositionRandomizer::AliHBTPositionRandomizer():
 fReader(0x0),
 fRandomizer(0x0),
 fModel(0),
 fAddToExistingPos(kFALSE),
 fOnlyParticlesFromVertex(kFALSE),
 fVX(0.0),
 fVY(0.0),
 fVZ(0.0)
{
//constructor
}
/*********************************************************************/

AliHBTPositionRandomizer::AliHBTPositionRandomizer(AliHBTReader* reader):
 fReader(reader),
 fRandomizer(new AliHBTRndmGaussBall(8.0)),
 fModel(0),
 fAddToExistingPos(kFALSE),
 fOnlyParticlesFromVertex(kFALSE),
 fVX(0.0),
 fVY(0.0),
 fVZ(0.0)
{
//constructor
} 
/*********************************************************************/

AliHBTEvent* AliHBTPositionRandomizer::GetParticleEvent() 
{
 // gets from fReader and randomizes current particle event
 if (fReader == 0x0) return 0x0;
 AliHBTEvent *e =  fReader->GetParticleEvent();
 if (e->IsRandomized() == kFALSE) Randomize(e);
 return e;
}
/*********************************************************************/

AliHBTEvent* AliHBTPositionRandomizer::GetTrackEvent() 
{
 // gets from fReader and randomizes current track event
 if (fReader == 0x0) return 0x0;
 AliHBTEvent *e =  fReader->GetTrackEvent();
 if (e->IsRandomized() == kFALSE) Randomize(e);
 return e;
}
/*********************************************************************/

Int_t AliHBTPositionRandomizer::Read(AliHBTRun* particles, AliHBTRun *tracks)
{
  //Reads all available events and randomizes them
  if (fReader == 0x0) return 1;
  Info("Randomize(AliHBTRun*)","");
  Int_t err = fReader->Read(particles,tracks);
  if (err) return err;
  Randomize(particles);
  return 0;
}
/*********************************************************************/

AliHBTEvent* AliHBTPositionRandomizer::GetParticleEvent(Int_t n)
{
//returns event n
 if (fReader == 0x0) return 0x0;
 AliHBTEvent *e =  fReader->GetParticleEvent(n);
 if (e->IsRandomized() == kFALSE) Randomize(e);
 return e;
}

/*********************************************************************/

void AliHBTPositionRandomizer::Randomize(AliHBTRun* run)
{
// randomizes postions of all particles in the run
  if (run == 0x0) return;
  Info("Randomize(AliHBTRun*)","");
  for (Int_t i = 0; i < run->GetNumberOfEvents(); i++)
   {
     Randomize(run->GetEvent(i));
   }
}
/*********************************************************************/
void AliHBTPositionRandomizer::Randomize(AliHBTEvent* event)
{
// randomizes postions of all particles in the event
  static const Double_t fmtocm = 1.e-13;
  Info("Randomize(AliHBTEvent*)","");
  if (event == 0x0) return;

  for (Int_t i = 0; i < event->GetNumberOfParticles(); i++)
   {
     AliHBTParticle* p = event->GetParticle(i);
     Double_t x,y,z,t=0.0;
     fRandomizer->Randomize(x,y,z,p);
     p->SetProductionVertex(x*fmtocm,y*fmtocm,z*fmtocm,t*fmtocm);
   }
  event->SetRandomized();
}
/*********************************************************************/

void AliHBTPositionRandomizer::SetGaussianBall(Double_t r)
{
  SetGaussianBall(r,r,r);
}
/*********************************************************************/

void AliHBTPositionRandomizer::SetGaussianBall(Double_t rx, Double_t ry, Double_t rz)
{
  delete fRandomizer;
  fRandomizer = new AliHBTRndmGaussBall(rx,ry,rz);
}
/*********************************************************************/

void AliHBTPositionRandomizer::SetCyllinderSurface(Double_t r, Double_t l)
{
  delete fRandomizer;
  fRandomizer = new  AliHBTRndmCyllSurf(r,l);
}
/*********************************************************************/

void AliHBTPositionRandomizer::SetEventVertex(Double_t x, Double_t y,Double_t z)
{
//sets event vertex position
  fVX = x;
  fVY = y;
  fVZ = z;
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
 fRz(0.0)
{
  //constructor
}
/*********************************************************************/

AliHBTRndmGaussBall::AliHBTRndmGaussBall(Float_t r):
 fRx(r),
 fRy(r),
 fRz(r)
{
  //constructor
}
/*********************************************************************/

AliHBTRndmGaussBall::AliHBTRndmGaussBall(Float_t rx, Float_t ry, Float_t rz):
 fRx(rx),
 fRy(ry),
 fRz(rz)
{
  //constructor
}
/*********************************************************************/

void AliHBTRndmGaussBall::Randomize(Double_t& x,Double_t& y,Double_t&z, AliHBTParticle*/*particle*/)
{
//randomizez gauss for each coordinate separately
  x = gRandom->Gaus(0.0,fRx);
  y = gRandom->Gaus(0.0,fRy);
  z = gRandom->Gaus(0.0,fRz);
}
/*********************************************************************/
//_____________________________________________________________________
///////////////////////////////////////////////////////////////////////
//                                                                   //
//  class AliHBTRndmGaussBall                                        //
//                                                                   //
///////////////////////////////////////////////////////////////////////

void AliHBTRndmCyllSurf::Randomize(Double_t& x,Double_t& y,Double_t&z, AliHBTParticle* particle)
{
   Double_t r = fR + gRandom->Gaus(0.0, 1.0);
   Double_t sf = r/particle->Pt();//scaling factor for position transformation ->
                             //we move direction of string momentum but legth defined by r
   x = sf*particle->Px();
   y = sf*particle->Py();
   z = gRandom->Uniform(-fL,fL);
  
}
