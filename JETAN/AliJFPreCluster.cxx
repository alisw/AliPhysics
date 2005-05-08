// $Id$

#include <Riostream.h>

#include <TMath.h>
#include <TParticle.h>

#include "AliJFPreCluster.h"

ClassImp(AliJFPreCluster)

AliJFPreCluster::AliJFPreCluster() : fPx(0),fPy(0),fPz(0),fE(0),fParticles("TParticle",0)
{
}

AliJFPreCluster::AliJFPreCluster(const AliJFPreCluster &copy) :
              fPx(0),fPy(0),fPz(0),fE(0),fParticles("TParticle",0)
{
  fPx=copy.GetPx();
  fPy=copy.GetPy();
  fPz=copy.GetPz();
  fE=copy.GetE();
  fParticles.Expand(copy.GetParticles()->GetEntries());
  fParticles.AddAll((TClonesArray*)copy.GetParticles());
}

AliJFPreCluster::AliJFPreCluster(const TParticle *p) :
              fPx(0),fPy(0),fPz(0),fE(0),fParticles("TParticle",0)
{
  fPx=p->Px();
  fPy=p->Py();
  fPz=p->Pz();
  fE=p->Energy();
  fParticles.Expand(1);
  fParticles[0]=new TParticle(*p);
}

AliJFPreCluster::AliJFPreCluster(Float_t px, Float_t py, Float_t pz, Float_t E, const TParticle *p) : 
              fPx(px),fPy(py),fPz(pz),fE(E),fParticles("TParticle",0)
{
  if(fE<0) fE=TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz);
  fParticles.Expand(1);
  fParticles[0]=new TParticle(*p);
}

AliJFPreCluster::AliJFPreCluster(Float_t px, Float_t py, Float_t pz, Float_t E, TClonesArray *parts) : 
              fPx(px),fPy(py),fPz(pz),fE(E),fParticles("TParticle",0)
{
  if(fE<0) fE=TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz);
  fParticles.Expand(parts->GetEntries());
  fParticles.AddAll(parts);
}

AliJFPreCluster::~AliJFPreCluster()
{
  fParticles.Delete();
}

AliJFPreCluster& AliJFPreCluster::operator=(const AliJFPreCluster &rhs)
{
  fPx=rhs.GetPx();
  fPy=rhs.GetPy();
  fPz=rhs.GetPz();
  fE=rhs.GetE();
  fParticles.Delete();
  fParticles.Expand(rhs.GetParticles()->GetEntries());
  fParticles.AddAll((TClonesArray*)rhs.GetParticles());

  return *this;
}

ostream& operator<<(ostream& o, const AliJFPreCluster &c)
{
  o << c.GetPx() << " " << c.GetPy() << " " << c.GetPz() << " " << c.GetE();

  return o;
}
    
/*
void AliJFPreCluster::SetValues(Float_t px, Float_t py, Float_t pz, Float_t E) 
{
  fPx=px;
  fPy=py;
  fPz=pz;
  fE=E;

  if(fE<0) fE=TMath::Sqrt(fPx*fPx+fPy+fPy+fPz*fPz);
}
*/
