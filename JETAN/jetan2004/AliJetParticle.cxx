// $Id$

//___________________________________________________________
/////////////////////////////////////////////////////////////
//
// class AliJetParticle
//
// loizides@ikf.uni-frankfurt.de
//
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TMath.h>
#include <TParticle.h>
#include "AliJetParticle.h"

ClassImp(AliJetParticle)

AliJetParticle::AliJetParticle() : 
  TObject()    
{
  Clear();
}

AliJetParticle::AliJetParticle(const AliJetParticle& in) :
  TObject(in)
{
  Clear();
  SetMomentum(in.fPx,in.fPy,in.fPz,in.fE);
  fIdxInEvent=in.fIdxInEvent;
  fLabel=in.fLabel;
  fType=in.fType;
  fNhits=in.fNhits;
}
 
AliJetParticle::AliJetParticle(const TParticle* p, Int_t idx, Int_t l, Int_t ncl) :
  TObject()
{
  Clear();
  SetMomentum(p->Px(),p->Py(),p->Pz(),p->Energy());
  fIdxInEvent=idx;
  fLabel=l;
  fType=(Int_t)p->GetWeight();
  fNhits=ncl;
}

AliJetParticle::AliJetParticle(Float_t px, Float_t py, Float_t pz, 
                               Float_t etot, Int_t idx, Int_t l, Int_t ncl) :
  TObject()
{
  Clear();
  SetMomentum(px,py,pz,etot);
  fIdxInEvent=idx;
  fLabel=l;
  fNhits=ncl;
}

AliJetParticle::AliJetParticle(Float_t px, Float_t py, Float_t pz, 
                               Float_t etot, Int_t idx, Int_t l, Int_t ncl,
			       Float_t pt, Float_t phi, Float_t eta) :
  TObject(),
  fPx(px),fPy(py),fPz(pz),
  fE(etot),fIdxInEvent(idx),
  fType(0),fLabel(l),fNhits(ncl),
  fPt(pt),fEta(eta),fPhi(phi)
{
}

TParticle* AliJetParticle::Particle() const
{
  TParticle *ret=new TParticle(0,0,0,0,0,0,fPx,fPy,fPz,fE,0,0,0,0);
  ret->SetWeight(fType);
  return ret;
}

void AliJetParticle::Calculate()
{
  //calculate values from px, py, pz
  //which are frequently needed

  fPt=TMath::Sqrt(fPx*fPx+fPy*fPy);
  const Float_t kp=P();
  fEta=0.5*TMath::Log((kp+fPz+1e-30)/(kp-fPz+1e-30)); 
  fPhi=TMath::Pi()+TMath::ATan2(-fPy,-fPx);
}

void AliJetParticle::Clear(Option_t* /*t*/)
{
  fPx=0.;
  fPy=0.;
  fPz=0.;
  fE=0.;
  fIdxInEvent=-1;
  fType=0;
  fLabel=-1;
  fNhits=0;
  fPt=0.;
  fEta=0.;
  fPhi=0.;
}

void AliJetParticle::Print(Option_t* /*t*/) const
{
  cout << fPt << " " << fEta << " " << fPhi << endl;
}

Int_t AliJetParticle::Compare(const TObject *obj) const
{
  Double_t val=((AliJetParticle*)obj)->Pt();

  if(fPt>val) return 1;
  else if (fPt<val) return -1;
  else return 0;
}
