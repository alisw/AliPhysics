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
#include <TObjArray.h>
#include "AliVAODParticle.h"

class TParticle;

class AliAOD: public TObject {
public:
  AliAOD();
  virtual ~AliAOD() { Reset(); }

  virtual void             SetOwner(Bool_t owner){fParticles.SetOwner(owner);}
  virtual TObjArray*       GetParticles() {return &fParticles;};
  virtual Int_t            GetNumberOfParticles() const  {return fParticles.GetEntriesFast();}
  virtual AliVAODParticle* GetParticle(Int_t index) const {return (AliVAODParticle*) fParticles[index];}
  virtual void             AddParticle(AliVAODParticle* particle)  {fParticles.Add(particle);};
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
private:
  TObjArray                fParticles;   // array of AOD particles, AliAOD is owner of particles
  Bool_t                   fIsRandomized;//flag indicating if positions of particles were randomized - used by HBTAN
  Double_t                 fPrimaryVertexX;//X position of the primary vertex
  Double_t                 fPrimaryVertexY;//Y position of the primary vertex
  Double_t                 fPrimaryVertexZ;//Z position of the primary vertex
  
  
  ClassDef(AliAOD,1)  // base class for AOD containers
};

#endif
