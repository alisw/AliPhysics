#ifndef ALIJETEVENTPARTICLES_H
#define ALIJETEVENTPARTICLES_H

/* $Id$ */

///////////////////////////////////////////////////////////////////
//
// class AliJetEventParticles
//
// loizides@ikf.uni-frankfurt.de 
//
///////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>

class TParticle;
class TClonesArray;
class AliJetParticle;

class AliJetEventParticles: public TObject
{
  public:
  AliJetEventParticles(Int_t size=1000);
  AliJetEventParticles(const AliJetEventParticles& source);
  virtual ~AliJetEventParticles();
    
  void SetVertex(Float_t v[3]){fVertexX=v[0];fVertexY=v[1];fVertexZ=v[2];}
  void SetVertex(Float_t v1,Float_t v2, Float_t v3){fVertexX=v1;fVertexY=v2;fVertexZ=v3;}
  void SetHeader(TString& s){fHeader=s;}
  void Reset(Int_t size=-1); //deletes all entries
  void Clear(Option_t *option="") {TObject::Clear(option);Reset();}
  void Set(const AliJetEventParticles& source);
  
  //adds particle to the event
  void AddParticle(AliJetParticle* p);  
  void AddParticle(const AliJetParticle* p); 
  void AddParticle(const TParticle* part,Int_t idx=-1, Int_t l=0, Int_t ncl=0); 
  void AddParticle(Float_t px, Float_t py, Float_t pz, Float_t etot, Int_t idx=-1, Int_t l=0, Int_t ncl=0);
  void AddParticle(Float_t px, Float_t py, Float_t pz, Float_t etot, Int_t idx, Int_t l, Int_t ncl,
		   Float_t pt, Float_t phi, Float_t eta);

  void AddSignal(const AliJetEventParticles& source);

  const AliJetParticle* GetParticle(Int_t n) const //gets particle without boundary check
    {return (const AliJetParticle*)fParticles->At(n);} 
  const AliJetParticle* GetParticleSafely(Int_t n); 
  Int_t GetNParticles()              const {return fNParticles;}
  const TClonesArray* GetParticles() const {return fParticles;}
  Float_t GetVertexX()               const {return fVertexX;}  
  Float_t GetVertexY()               const {return fVertexY;}  
  Float_t GetVertexZ()               const {return fVertexZ;}  

  Int_t    Trials()         const {return fTrials;}
  Int_t    NTriggerJets()   const {return fNJets;}
  Int_t    NUQTriggerJets() const {return fNUQJets;}
  void     TriggerJet(Int_t i, Float_t p[4]) const;
  void     UQJet(Int_t i, Float_t p[4])      const;
  void     TriggerJet(Int_t i, Float_t &p1, Float_t &p2, Float_t &p3, Float_t &E)          const;
  void     UQJet(Int_t i, Float_t &p1, Float_t &p2, Float_t &p3, Float_t &E)               const;
  void     Hard(Int_t i, Float_t &p1, Float_t &p2, Float_t &p3, Float_t &E, Float_t &type) const;
  void     Hard(Int_t i, Float_t p[4], Float_t &type) const;

  void     AddJet(Float_t px, Float_t py, Float_t pz, Float_t e);
  void     AddUQJet(Float_t px, Float_t py, Float_t pz, Float_t e);
  void     AddJet(Float_t p[4]);
  void     AddUQJet(Float_t p[4]);
  void     AddHard(Int_t i,Float_t px, Float_t py, Float_t pz, Float_t e, Float_t type);

  Double_t GetXJet() const {return fXJet;}
  Double_t GetYJet() const {return fYJet;}    
  void     GetZQuench(Double_t z[4]) const;
  TString  GetHeader() const {return fHeader;}
  Float_t  GetImpact()  const {return fImpact;}
  Int_t    GetNhard()   const {return fNHardScatters;}
  Int_t    GetNpart()   const {return fNwNwColl;}
  Int_t    GetEventNr() const {return fEventNr;}

  void     SetXYJet(Double_t x, Double_t y); 
  void     SetZQuench(Double_t z[4]);
  void     SetTrials(Int_t trials) {fTrials = trials;}

  void     SetImpact(Float_t b){fImpact=b;}
  void     SetNhard(Int_t n){fNHardScatters=n;}
  void     SetNpart(Int_t n){fNwNwColl=n;}
  void SetEventNr(Int_t n){fEventNr=n;}
  void     Print(Option_t *t="") const;

  protected:
  TString fHeader;          //   event description
  Int_t fNParticles;        //   number of particles read
  TClonesArray *fParticles; //-> particles in event

  Float_t fVertexX;         // vertex x
  Float_t fVertexY;         // vertex y
  Float_t fVertexZ;         // vertex z

  Int_t    fTrials;         // Number of trials to fulfill trigger condition
  Int_t    fNJets;          // Number of triggered jets
  Int_t    fNUQJets;        // Number of unquenched
  Double_t fXJet;           // Jet production point (x)
  Double_t fYJet;           // Jet production point (y)
  Float_t  fJets[4][10];    // Trigger jets
  Float_t  fUQJets[4][10];  // Unquenched trigger jets
  Float_t  fHard[5][2];     // Hard partons
  Double_t fZquench[4];     // Quenching fraction

  Float_t fImpact;          // impact parameter
  Int_t   fNHardScatters;   // number of hard scatterings
  Int_t   fNwNwColl;        // number of Nwounded-Nwounded collisions

  Int_t fEventNr;           // simple counter to check event number

  ClassDef(AliJetEventParticles,6) //class AliJetEventParticles
};
#endif
