#ifndef ALIHBTPOSITIONRANDOMIZER_H
#define ALIHBTPOSITIONRANDOMIZER_H
//___________________________________________________
////////////////////////////////////////////////////////////////////////////////
// 
// class AliHBTPositionRandomizer
//
// These class randomizes particle vertex positions
// Piotr.Skowronski@cern.ch
// 
////////////////////////////////////////////////////////////////////////////////

#include "AliHBTReader.h"

class AliHBTRndm;
class AliHBTEvent;
class AliHBTRun;
class AliHBTParticle;
class TH1I;

class AliHBTPositionRandomizer: public AliHBTReader
{
 public:
   enum EModelTypes{kGausBall,kCylinder,kCylinderSurf};
   AliHBTPositionRandomizer();
   AliHBTPositionRandomizer(AliHBTReader* reader);
   AliHBTPositionRandomizer(const AliHBTPositionRandomizer& in);
   
   virtual ~AliHBTPositionRandomizer();

   const AliHBTPositionRandomizer& operator=(const AliHBTPositionRandomizer& in);
   
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
   virtual TH1I* GetTrackCounter() const {return (fReader)?fReader->GetTrackCounter():0x0;}
   virtual void  WriteTrackCounter() const {if(fReader) fReader->WriteTrackCounter();}

   void Randomize(AliHBTEvent* event) const;
   void Randomize(AliHBTRun* run) const;
   void SetEventVertex(Double_t x, Double_t y,Double_t z);
   
   void SetGaussianBall(Double_t r);
   void SetGaussianBall(Double_t rx, Double_t ry, Double_t rz);
   void SetCyllinderSurface(Double_t r, Double_t l);
   
 protected:
   void Randomize(Double_t& x,Double_t& y,Double_t&z, AliHBTParticle*p);
   Int_t ReadNext(){return (fReader)?fReader->Next():1;}
   
 private:
   AliHBTReader* fReader;      // Pointer to reader
   AliHBTRndm*   fRandomizer;  // Pointer to class that performs randomization according to some model
   
   Int_t    fModel;            //Defines what model is used
   
   Bool_t   fAddToExistingPos;  //Determines if randomized position should be added to previous one, or overwrite old one
   Bool_t   fOnlyParticlesFromVertex; //Determines if randomization should be performed for particles from vertex

   Double_t fVX; //vertex position
   Double_t fVY; //vertex position
   Double_t fVZ; //vertex position
   
   ClassDef(AliHBTPositionRandomizer,1)
};

class AliHBTRndm: public TObject
{
  public:
   AliHBTRndm(){}
   virtual ~AliHBTRndm(){}
   virtual void Randomize(Double_t& x,Double_t& y,Double_t&z, AliHBTParticle*p) const = 0;
   ClassDef(AliHBTRndm,1)
};

class AliHBTRndmGaussBall: public AliHBTRndm
{
  public:
   AliHBTRndmGaussBall();
   AliHBTRndmGaussBall(Float_t r);
   AliHBTRndmGaussBall(Float_t rx, Float_t ry, Float_t rz);
   virtual ~AliHBTRndmGaussBall(){}
   void Randomize(Double_t& x,Double_t& y,Double_t&z, AliHBTParticle*/*particle*/) const;
  private:
   Float_t fRx; //Dispertion in x 
   Float_t fRy; //Dispertion in y
   Float_t fRz; //Dispertion in z
   ClassDef(AliHBTRndmGaussBall,1)
};

class AliHBTRndmCyllSurf: public AliHBTRndm
{
  public:
   AliHBTRndmCyllSurf():fR(0.0){}
   AliHBTRndmCyllSurf(Float_t r, Float_t l):fR(r),fL(l){}
   virtual ~AliHBTRndmCyllSurf(){}
   
   void Randomize(Double_t& x,Double_t& y,Double_t&z, AliHBTParticle* particle) const;
  private:
   Float_t fR; //Redius of cylinder
   Float_t fL; //Length of cylinder
 
   ClassDef(AliHBTRndmCyllSurf,1)
};

#endif

