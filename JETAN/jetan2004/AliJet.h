#ifndef ALIJET_H
#define ALIJET_H

/* $Id$ */

///////////////////////////////////////////////////////////////////
//
// class AliJet
//
// loizides@ikf.uni-frankfurt.de
//
///////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TString.h>
class TParticle;
class AliJetParticle;

class AliJet: public TObject
{
  public:
#if 0
  AliJetEvent(Int_t size=1000);
  AliJetEvent(const AliJetEvent& source);
  virtual ~AliJetEvent();
    
  void SetVertex(Float_t v[3]){fVertexX=v[0];fVertexY=v[1];fVertexZ=v[2];}
  void SetVertex(Float_t v1,Float_t v2, Float_t v3){fVertexX=v1;fVertexY=v2;fVertexZ=v3;}
  void SetHeader(TString& s){fHeader=s;}
  void Reset(Int_t size=-1); //deletes all entries
  
  //adds particle to the event
  void AddParticle(AliJetParticle* p);  
  void AddParticle(const AliJetParticle* p); 
  void AddParticle(const TParticle* part,Int_t idx=-1, Int_t l=0); 
  void AddParticle(Float_t px, Float_t py, Float_t pz, Float_t etot, Int_t idx=-1, Int_t l=0);
  void AddParticle(Float_t px, Float_t py, Float_t pz, Float_t etot, Int_t idx, Int_t l,
		   Float_t pt, Float_t phi, Float_t eta);

  const AliJetParticle* GetParticle(Int_t n) //gets particle without boundary check
    {return (const AliJetParticle*)fParticles->At(n);} 
  const AliJetParticle* GetParticleSafely(Int_t n); 
  Int_t GetNParticles()              const {return fNParticles;}
  const TClonesArray* GetParticles() const {return fParticles;}
  Float_t GetVertexX()               const {return fVertexX;}  
  Float_t GetVertexY()               const {return fVertexY;}  
  Float_t GetVertexZ()               const {return fVertexZ;}  

  void Print(Option_t *t="") const;

  protected:
  TString fHeader;          //   event description
  Int_t fNParticles;        //   number of particles read
  TClonesArray *fParticles; //-> particles in event

  Float_t fVertexX; //vertex x
  Float_t fVertexY; //vertex y
  Float_t fVertexZ; //vertex z
#endif
  ClassDef(AliJet,1) //class AliJet
};
#endif
