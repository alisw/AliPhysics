#ifndef ALIHBTPOSITIONRANDOMIZER_H
#define ALIHBTPOSITIONRANDOMIZER_H

#include "AliHBTReader.h"

class AliHBTRndm;
class AliHBTEvent;
class AliHBTRun;
class AliHBTParticle;

class AliHBTPositionRandomizer: public AliHBTReader
{
 public:
   enum EModelTypes{kGausBall,kCylinder,kCylinderSurf};
   AliHBTPositionRandomizer();
   AliHBTPositionRandomizer(AliHBTReader* reader);
   
   Int_t  Next(){return (fReader)?fReader->Next():1;}
   void   Rewind(){if(fReader) fReader->Rewind();}
   
   Bool_t ReadsTracks() const {return (fReader)?fReader->ReadsTracks():kFALSE;}
   Bool_t ReadsParticles() const {return (fReader)?fReader->ReadsParticles():kFALSE;}
   
   Int_t  Read(AliHBTRun* particles, AliHBTRun *tracks);
   
   AliHBTEvent* GetParticleEvent() ;
   AliHBTEvent* GetTrackEvent() ;

   AliHBTEvent* GetParticleEvent(Int_t n);
   AliHBTEvent* GetTrackEvent(Int_t n){return (fReader)?fReader->GetTrackEvent(n):0x0;}
   Int_t GetNumberOfPartEvents(){return (fReader)?fReader->GetNumberOfPartEvents():0;}
   Int_t GetNumberOfTrackEvents(){return (fReader)?fReader->GetNumberOfTrackEvents():0;}
   
   void Randomize(AliHBTEvent* event);
   void Randomize(AliHBTRun* run);
   void SetEventVertex(Double_t x, Double_t y,Double_t z);
   
   void SetGaussianBall(Double_t r);
   void SetGaussianBall(Double_t rx, Double_t ry, Double_t rz);
   void SetCyllinderSurface(Double_t r, Double_t l);
   
 protected:
   void Randomize(Double_t& x,Double_t& y,Double_t&z, AliHBTParticle*p);
   Int_t ReadNext(){return (fReader)?fReader->Next():1;}
 private:
   AliHBTReader* fReader;
   AliHBTRndm*   fRandomizer;
   
   Int_t    fModel;
   
   Bool_t   fAddToExistingPos;
   Bool_t   fOnlyParticlesFromVertex;

   Double_t fVX; //vertex position
   Double_t fVY; //vertex position
   Double_t fVZ; //vertex position
   
   ClassDef(AliHBTPositionRandomizer,1)
};

class AliHBTRndm: public TObject
{
  public:
   AliHBTRndm(){}
   ~AliHBTRndm(){}
   virtual void Randomize(Double_t& x,Double_t& y,Double_t&z, AliHBTParticle*p) = 0;
   ClassDef(AliHBTRndm,1)
};

class AliHBTRndmGaussBall: public AliHBTRndm
{
  public:
   AliHBTRndmGaussBall();
   AliHBTRndmGaussBall(Float_t r);
   AliHBTRndmGaussBall(Float_t rx, Float_t ry, Float_t rz);
   ~AliHBTRndmGaussBall(){}
   void Randomize(Double_t& x,Double_t& y,Double_t&z, AliHBTParticle*/*particle*/);
  private:
   Float_t fRx;
   Float_t fRy;
   Float_t fRz;
   ClassDef(AliHBTRndmGaussBall,1)
};

class AliHBTRndmCyllSurf: public AliHBTRndm
{
  public:
   AliHBTRndmCyllSurf():fR(0.0){}
   AliHBTRndmCyllSurf(Float_t r, Float_t l):fR(r),fL(l){}
   ~AliHBTRndmCyllSurf(){}
   
   void Randomize(Double_t& x,Double_t& y,Double_t&z, AliHBTParticle* particle);
  private:
   Float_t fR;
   Float_t fL;
 
   ClassDef(AliHBTRndmCyllSurf,1)
};

#endif

