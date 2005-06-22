#ifndef ALIJETEVENT_H
#define ALIJETEVENT_H

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// class AliJetEvent                                             //
//                                                               //
// loizides@ikf.uni-frankfurt.de                                 //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <TString.h>
class TParticle;
class AliJetParticle;
class AliJet;

class AliJetEvent: public TObject
{
  public:
#if 0
  AliJetEvent(Int_t size=1000);
  AliJetEvent(const AliJetEvent& source);
  virtual ~AliJetEvent();
    
  void SetHeader(TString& s){fHeader=s;}
  //void Reset(Int_t size=-1); //deletes all entries
  
  //adds particle to the event
  void AddJet(AliJet* j);  
  void AddJet(const AliJet* j); 

  const AliJet* GetJet(Int_t n) //gets jet without boundary check
    {return (const AliJet*)fJets->At(n);} 
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
  ClassDef(AliJetEvent,1) //class AliJetEvent
};
#endif
