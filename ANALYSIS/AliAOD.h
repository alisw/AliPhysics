#ifndef ALIAOD_H
#define ALIAOD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// base class for AOD containers
//
/////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TClonesArray.h>
#include "AliVAODParticle.h"

class TParticle;

class AliAOD: public TObject {
public:
  AliAOD();
  virtual ~AliAOD();
  
  AliAOD(const AliAOD& in);
  virtual AliAOD& operator=(const AliAOD& in);
  virtual void             CopyData(AliAOD* aod);//Copys all data from aod, but leaves local type of particles
  
  virtual TClonesArray*    GetParticles() const {return fParticles;} 
  virtual void             SetParticleClassName(const char* classname);
  virtual void             SetParticleClass(TClass* pclass);
  
  virtual Int_t            GetNumberOfParticles() const  {return (fParticles)?fParticles->GetEntriesFast():0;}
  virtual AliVAODParticle* GetParticle(Int_t index) const {return  (fParticles)?(AliVAODParticle*)fParticles->At(index):0x0;}
  virtual void             AddParticle(AliVAODParticle* particle);
  virtual void             AddParticle(TParticle* part, Int_t idx); //adds particle to the event
  virtual void             AddParticle(Int_t pdg, Int_t idx, Double_t px, Double_t py, Double_t pz, Double_t etot,
                                       Double_t vx, Double_t vy, Double_t vz, Double_t time);
  
  virtual void             Reset();
  
  void                     SwapParticles(Int_t i, Int_t j);//swaps particles positions; used by AliReader::Blend
  Bool_t                   IsRandomized() const {return fIsRandomized;}
  void                     SetRandomized(Bool_t flag = kTRUE){fIsRandomized = flag;}
  
  void                     GetPrimaryVertex(Double_t&x, Double_t&y, Double_t&z);
  void                     SetPrimaryVertex(Double_t x, Double_t y, Double_t z);
  
  Int_t                    GetNumberOfCharged(Double_t etamin = -10.0, Double_t etamax = 10.0) const;
  void                     Move(Double_t x, Double_t y, Double_t z);//moves all spacial coordinates about this vector
  virtual void             SetOwner(Bool_t owner);
  virtual void             Print(Option_t* /*option*/ = 0);
  const TClass*            GetParticleClass() const {return fParticleClass;}
private:
  TClonesArray            *fParticles;   // array of AOD particles, AliAOD is owner of particles
  Bool_t                   fIsRandomized;//flag indicating if positions of particles were randomized - used by HBTAN
  Double_t                 fPrimaryVertexX;//X position of the primary vertex
  Double_t                 fPrimaryVertexY;//Y position of the primary vertex
  Double_t                 fPrimaryVertexZ;//Z position of the primary vertex
  TClass*                  fParticleClass;//object that defines type of the particle       
  
  ClassDef(AliAOD,1)  // base class for AOD containers
};

#endif
