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

#include "AliReader.h"
class AliHBTRndm;
class TH1I;

class AliHBTPositionRandomizer: public AliReader
{
 public:
   enum EModelTypes{kGausBall,kCylinder,kCylinderSurf,kEllipse};
   AliHBTPositionRandomizer();
   AliHBTPositionRandomizer(AliReader* reader);
   AliHBTPositionRandomizer(const AliHBTPositionRandomizer& in);
   
   virtual ~AliHBTPositionRandomizer();

   AliHBTPositionRandomizer& operator=(const AliHBTPositionRandomizer& in);
   
   Int_t  Next(){return (fReader)?fReader->Next():1;}
   void   Rewind(){if(fReader) fReader->Rewind();}
   
   Bool_t ReadsRec() const {return (fReader)?fReader->ReadsRec():kFALSE;}
   Bool_t ReadsSim() const {return (fReader)?fReader->ReadsSim():kFALSE;}
   
   AliAOD* GetEventSim() const ;
   AliAOD* GetEventRec() const ;

   AliAOD* GetEventSim(Int_t n);
   AliAOD* GetEventRec(Int_t n){return (fReader)?fReader->GetEventRec(n):0x0;}
   
   Int_t GetNumberOfSimEvents(){return (fReader)?fReader->GetNumberOfSimEvents():0;}
   Int_t GetNumberOfRecEvents(){return (fReader)?fReader->GetNumberOfRecEvents():0;}
   virtual TH1I* GetTrackCounter() const {return (fReader)?fReader->GetTrackCounter():0x0;}
   virtual void  WriteTrackCounter() const {if(fReader) fReader->WriteTrackCounter();}

   void Randomize(AliAOD* event) const;
   void SetEventVertex(Double_t x, Double_t y,Double_t z);
   
   void SetRandomizer(Int_t pid,AliHBTRndm* rndm);
   
   void SetGaussianBall(Int_t pid, Double_t r, Double_t meantime, Double_t sigmatime);
   void SetGaussianBall(Int_t pid, Double_t rx, Double_t ry, Double_t rz, Double_t meantime, Double_t sigmatime);
   void SetCyllinderSurface(Int_t pid, Double_t r, Double_t l);
   void SetEllipse(Int_t pid, Double_t rmin, Double_t rmax);
   
   void AddToPosition(Bool_t flag){fAddToExistingPos = flag;}
   void RandomizeTracks(Bool_t flag){fRandomizeTracks = flag;}
   
   AliHBTRndm* GetRandomizer(Int_t pdg) const;
   Int_t GetRandomizerIndex(Int_t pdg) const;
   
 protected:
   void Randomize(Double_t& x,Double_t& y,Double_t&z,AliVAODParticle*p);
   Int_t ReadNext(){return (fReader)?fReader->Next():1;}
   
 private:
   AliReader*   fReader;      // Pointer to reader
   AliHBTRndm*  fDefaultRandomizer;  // Pointer to class that performs randomization according to some model - default one
   TObjArray*   fRandomizers;//array with randomizers - each particle type can have different randomization parameters/model
   Int_t        fNPid;//number of randomizers defined in fPid and fRandomizers
   Int_t*       fPids;//[fgkNumberOfPids]
   
   Bool_t       fAddToExistingPos;  //Determines if randomized position should be added to previous one, or overwrite old one
   Bool_t       fOnlyParticlesFromVertex; //Determines if randomization should be performed for particles from vertex

   Bool_t       fRandomizeTracks; //Determines if tracks should also be randimized 
   
   Double_t     fVX; //vertex position
   Double_t     fVY; //vertex position
   Double_t     fVZ; //vertex position
   
   static const Int_t fgkNumberOfPids;//size of fPid array
   
   ClassDef(AliHBTPositionRandomizer,1)
};

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

class AliHBTRndm: public TObject
{
  public:
   AliHBTRndm(){}
   virtual ~AliHBTRndm(){}
   virtual void Randomize(Double_t& x, Double_t& y, Double_t& z, Double_t& t, AliVAODParticle*p) const = 0;
     ClassDef(AliHBTRndm,1)
};

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

class AliHBTRndmGaussBall: public AliHBTRndm
{
  public:
   AliHBTRndmGaussBall();
   AliHBTRndmGaussBall(Float_t r, Double_t meantime, Double_t sigmatime);
   AliHBTRndmGaussBall(Float_t rx, Float_t ry, Float_t rz, Double_t meantime, Double_t sigmatime);
   virtual ~AliHBTRndmGaussBall(){}
   void Randomize(Double_t& x,Double_t& y,Double_t&z, Double_t& t, AliVAODParticle*/*particle*/) const;
  private:
   Float_t fRx; //Dispertion in x 
   Float_t fRy; //Dispertion in y
   Float_t fRz; //Dispertion in z
   Float_t fTmean; //Mean emision time in t
   Float_t fTsigma; //Dispertion in t
   ClassDef(AliHBTRndmGaussBall,1)
};

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

class AliHBTRndmCyllSurf: public AliHBTRndm
{
  public:
   AliHBTRndmCyllSurf():fR(0.0){}
   AliHBTRndmCyllSurf(Float_t r, Float_t l):fR(r),fL(l){}
   virtual ~AliHBTRndmCyllSurf(){}
   
   void Randomize(Double_t& x,Double_t& y,Double_t&z, Double_t& /*t*/t, AliVAODParticle* particle) const;
  private:
   Float_t fR; //Redius of cylinder
   Float_t fL; //Length of cylinder
 
   ClassDef(AliHBTRndmCyllSurf,1)
};

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

class AliHBTRndmEllipse: public AliHBTRndm
{
  public:
   AliHBTRndmEllipse():fRmin(0.),fRmax(0.){};
   AliHBTRndmEllipse(Float_t rmin, Float_t rmax);
   virtual ~AliHBTRndmEllipse(){}
   
   void Randomize(Double_t& x,Double_t& y, Double_t& z,  Double_t& , AliVAODParticle* particle) const;
  private:
   Float_t fRmin; //Radius in x direction
   Float_t fRmax; //Radius in y direction
 
   ClassDef(AliHBTRndmEllipse,2)
};



#endif

