#include "AliHBTPositionRandomizer.h"
//___________________________________________________
////////////////////////////////////////////////////////////////////////////////
// 
// class AliHBTPositionRandomizer
//
// These class randomizes particle vertex positions
// Piotr.Skowronski@cern.ch
//
////////////////////////////////////////////////////////////////////////////////

#include <TRandom.h>
#include "AliAOD.h"
#include "AliVAODParticle.h"


ClassImp(AliHBTPositionRandomizer)

/*********************************************************************/

AliHBTPositionRandomizer::AliHBTPositionRandomizer():
 fReader(0x0),
 fRandomizer(0x0),
 fModel(0),
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
 fRandomizer(new AliHBTRndmGaussBall(8.0)),
 fModel(0),
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

AliHBTPositionRandomizer::AliHBTPositionRandomizer(const AliHBTPositionRandomizer& in):
 AliReader(in),
 fReader(),
 fRandomizer(0x0),
 fModel(0),
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
  delete fRandomizer;
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
 if (e->IsRandomized() == kFALSE) Randomize(e);
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
 if (fRandomizeTracks) if (e->IsRandomized() == kFALSE) Randomize(e);
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
  if (AliVAODParticle::GetDebug() > 5) Info("Randomize(AliAOD*)","");
  if (event == 0x0) return;

  for (Int_t i = 0; i < event->GetNumberOfParticles(); i++)
   {
     AliVAODParticle* p = event->GetParticle(i);
     Double_t x,y,z,t=0.0;
     fRandomizer->Randomize(x,y,z,p);
     
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

void AliHBTPositionRandomizer::SetGaussianBall(Double_t r)
{
 //Sets Gaussian Ball Model
  SetGaussianBall(r,r,r);
}
/*********************************************************************/

void AliHBTPositionRandomizer::SetGaussianBall(Double_t rx, Double_t ry, Double_t rz)
{
 //Sets Gaussian Ball Model
  delete fRandomizer;
  fRandomizer = new AliHBTRndmGaussBall(rx,ry,rz);
}
/*********************************************************************/

void AliHBTPositionRandomizer::SetCyllinderSurface(Double_t r, Double_t l)
{
 //Sets Cylinder Surface Model
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

void AliHBTRndmGaussBall::Randomize(Double_t& x,Double_t& y,Double_t&z, AliVAODParticle*/*particle*/) const
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

void AliHBTRndmCyllSurf::Randomize(Double_t& x,Double_t& y,Double_t&z, AliVAODParticle* particle) const
{
//Randomizes x,y,z
   Double_t r = fR + gRandom->Gaus(0.0, 1.0);
   Double_t sf = r/particle->Pt();//scaling factor for position transformation ->
                             //we move direction of string momentum but legth defined by r
   x = sf*particle->Px();
   y = sf*particle->Py();
   z = gRandom->Uniform(-fL,fL);
}
