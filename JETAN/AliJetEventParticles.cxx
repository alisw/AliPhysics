// $Id$

//__________________________________________________________
///////////////////////////////////////////////////////////////////
//
// class AliJetEventParticles
//
// loizides@ikf.uni-frankfurt.de
//
///////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include "AliJetParticle.h"
#include "AliJetEventParticles.h"

ClassImp(AliJetEventParticles)

/**************************************************************************/ 
 
AliJetEventParticles::AliJetEventParticles(Int_t size) :
  fNParticles(0),
  fParticles(new TClonesArray("AliJetParticle",size)),
  fVertexX(0.),
  fVertexY(0.),
  fVertexZ(0.)
{
  //default constructor   
}

/**************************************************************************/ 
 
AliJetEventParticles::AliJetEventParticles(const AliJetEventParticles& source) :
  TObject(source), 
  fNParticles(source.fNParticles),
  fParticles(new TClonesArray("AliJetEventParticles",source.fNParticles)),
  fVertexX(0.),
  fVertexY(0.),
  fVertexZ(0.)
{
  //copy constructor
  for(Int_t i =0; i<fNParticles; i++)
    {
      const AliJetParticle *kjp=(const AliJetParticle *)source.fParticles->At(i);
      fParticles->AddAt(new AliJetParticle(*(kjp)),i);
    }
}

/**************************************************************************/ 

AliJetEventParticles::~AliJetEventParticles()
{
  //destructor   
  Reset();
  delete fParticles;
}


/**************************************************************************/ 

void  AliJetEventParticles::Reset(Int_t size)
{
  //deletes all particles from the event

  if(fParticles) fParticles->Clear();
  if(size>=0) fParticles->Expand(size);
  fNParticles = 0;
} 

/**************************************************************************/ 

void AliJetEventParticles::AddParticle(AliJetParticle* part)
{
  //Adds new particle to the event
  fParticles->AddAt(part,fNParticles++);
}

/**************************************************************************/ 

void AliJetEventParticles::AddParticle(const AliJetParticle* part)
{
  //Adds new particle to the event
  new((*fParticles)[fNParticles++]) AliJetParticle(*part);
}

/**************************************************************************/ 

void AliJetEventParticles::AddParticle(const TParticle* part,Int_t idx, Int_t l)
{
  //Adds new particle to the event
  new((*fParticles)[fNParticles++]) AliJetParticle(part,idx,l);
}

/**************************************************************************/ 

void AliJetEventParticles::AddParticle(Float_t px, Float_t py, Float_t pz, 
                              Float_t etot, Int_t idx, Int_t l)
{
  //Adds new particle to the event
  new((*fParticles)[fNParticles++]) AliJetParticle(px,py,pz,etot,idx,l); 
}

/**************************************************************************/ 

void AliJetEventParticles::AddParticle(Float_t px, Float_t py, Float_t pz, Float_t etot, Int_t idx, Int_t l,
		              Float_t pt, Float_t phi, Float_t eta)
{
  //Adds new particle to the event
  new((*fParticles)[fNParticles++]) AliJetParticle(px,py,pz,etot,idx,l,pt,phi,eta); 
}


/**************************************************************************/ 

const AliJetParticle* AliJetEventParticles::GetParticleSafely(Int_t n)
{
  //returns nth particle with range check
  if( (n<0) || (fNParticles<=n) ) return 0;
  return (const AliJetParticle*)fParticles->At(n);
}

/**************************************************************************/ 

void AliJetEventParticles::Print(Option_t* /*t*/) const
{
  cout << "--- AliJetEventParticles ---" << endl;
  if(fHeader.Length()) cout << fHeader.Data() << endl;
  cout << "no of particles: " << fNParticles << endl;
}
