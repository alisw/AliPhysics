//$Id$

#include <Riostream.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <TMath.h>

#include "AliTkTowerV2.h"
#include "AliTkConeJet.h"

#define Ntowers__ 100

AliTkConeJet::AliTkConeJet() : TObject(), 
                         fEta(-999),fPhi(-999),fEt(-999),fType(0),
                         fNTowers(0),fTowers(0),fPtCut(0.),
			 fCEta(-999),fCPhi(-999),fCEt(-999),
			 fPLength(0.),fXAxis(0.),fYAxis(0.),fZAxis(0.),
			 fPtLead(0.),fLeadPart(0),fNParts(0),fParts(0)
{
  fTowers = new TClonesArray("AliTkTowerV2",Ntowers__);
}

AliTkConeJet::AliTkConeJet(Float_t pt, Float_t eta, Float_t phi,Int_t type) 
                        : TObject(), 
			  fEta(eta),fPhi(phi),fEt(pt),fType(type), 
			  fNTowers(0),fTowers(0),fPtCut(0.),
			  fCEta(eta),fCPhi(phi),fCEt(pt),
			  fPLength(0.),fXAxis(0.),fYAxis(0.),fZAxis(0.),
			  fPtLead(0.),fLeadPart(0),fNParts(0),fParts(0)
{
  fTowers = new TClonesArray("AliTkTowerV2",0);
}

AliTkConeJet::AliTkConeJet(const AliTkConeJet &j) 
                        : TObject(), 
			  fEta(-999), fPhi(-999), fEt(-999),fType(0), 
			  fNTowers(0),fTowers(0),fPtCut(0.),
			  fCEta(-999),fCPhi(-999),fCEt(-999),
			  fPLength(0.),fXAxis(0.),fYAxis(0.),fZAxis(0.),
			  fPtLead(0.),fLeadPart(0),fNParts(0),fParts(0)
{
  fEta = j.getEta();
  fPhi = j.getPhi();
  fEt = j.getEt();
  fCEta = j.getCEta();
  fCPhi = j.getCPhi();
  fCEt = j.getCEt();
  fNTowers = j.getNTowers();
  fPtCut = j.getPtCut();
  fCEta = j.getCEta();
  fCPhi = j.getCPhi();
  fCEt = j.getCEt();
  fPLength = j.getPLength();
  fXAxis = j.getXAxis();   
  fYAxis = j.getYAxis();
  fZAxis = j.getZAxis();
  fPtLead = j.getPtLead();
  fLeadPart = new TParticle(*j.getLeadPart());
  fNParts = j.getNParts();

  // need to build a copy of the old TClonesArray...
  fTowers = new TClonesArray("AliTkTowerV2",Ntowers__);
  TClonesArray *otherTowers = j.getTowerList();
  AliTkTowerV2 *myTower;
  TIterator *iter = otherTowers->MakeIterator();
  Int_t i =0;
  while ((myTower = (AliTkTowerV2 *) iter->Next()) != NULL) {
    new ((*fTowers)[i]) AliTkTowerV2(*myTower);
    i++;
  }
  if(i!=fNTowers)
    cerr << "AliTkConeJet: should not happen " << i << " " << fNTowers << endl;
  delete iter;

  // need to build a copy of the old TClonesArray...
  fParts = new TClonesArray("TParticle",j.getNParts());
  TClonesArray *otherParts = j.getParts();
  TParticle *myParticle;
  iter = otherParts->MakeIterator();
  i =0;
  while ((myParticle = (TParticle *) iter->Next()) != NULL) {
    new ((*fParts)[i]) TParticle(*myParticle);
    i++;
  }
  if(i!=fNParts)
    cerr << "AliTkConeJet: should not happen " << i << " " << fNParts << endl;
  delete iter;
}
 
AliTkConeJet::~AliTkConeJet() 
{
  delete fTowers;
  if(fLeadPart) delete fLeadPart;
  if(fParts) delete fParts;
}

void AliTkConeJet::addTower(AliTkTowerV2 *tower) 
{
  new ((*fTowers)[fNTowers]) AliTkTowerV2(*tower);
  fNTowers++; 
}

void AliTkConeJet::Clear(Option_t *option) 
{
  TObject::Clear(option);
  fTowers->Delete();
  fNTowers=0;
  if(fParts) delete fParts;
  fParts=0;
  if(fLeadPart) delete fLeadPart;
  fLeadPart=0;
  fEta=-999; 
  fPhi=-999;
  fEt=-999;
  fPtCut=0.;
  fCEta=-999; 
  fCPhi=-999;
  fCEt=-999;
  fPLength=0.;
  fXAxis=0.;
  fYAxis=0.;
  fZAxis=0.;
  fPtLead=0.;
  fLeadPart=0;
  fNParts=0;
  fParts=0;
}

void AliTkConeJet::calculateFromParticles(Float_t &et, Float_t &eta, Float_t &phi, Float_t ptcut)
{
  Float_t px=0.,py=0.,pz=0.;
  TIterator *titer = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)titer->Next()) != NULL) {
    TClonesArray *particles = tower->getParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      px+=particle->Px();
      py+=particle->Py();
      pz+=particle->Pz();
    }
    delete part;
  }
  delete titer;

  et = TMath::Sqrt(px*px+py*py);
  Float_t p = TMath::Sqrt(px*px+py*py+pz*pz);
  Float_t theta = (pz==0)?TMath::PiOver2():TMath::ACos(pz/p);
  phi = TMath::Pi()+TMath::ATan2(-py,-px);
  eta = -TMath::Log(TMath::Tan(theta/2.));
}

/* or take e = p_t * cosh(eta) */
Float_t AliTkConeJet::getE(Float_t ptcut) const 
{
  Float_t e=0.;

  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      e+=particle->P();
    }
    delete part;
  }
  delete iter;
  return e;
}

/* or take e = p_t * cosh(eta) */
Float_t AliTkConeJet::getECharged(Float_t ptcut) const 
{
  Float_t e=0.;

  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getChargedParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      e+=particle->P();
    }
    particles->Delete();
    delete particles;
    delete part;
  }
  delete iter;
  return e;
}

Float_t AliTkConeJet::getEtCharged() const 
{
  Float_t EtCharged = 0;
  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    EtCharged += tower->getEtCharged();
  }
  delete iter;
  return EtCharged;
}

/* use to get fake rate if 
   in mixed case particles were marked */
Float_t AliTkConeJet::getEtChargedMarked(Float_t ptcut) const 
{
  Float_t et=0.;

  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getChargedParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->GetWeight()>-100) continue; //(-123)
      if(particle->Pt()<ptcut) continue;
      et+=particle->Pt();
    }
    particles->Delete();
    delete particles;
   delete part;
  }
  delete iter;
  return et;
}

/* use to get fake rate if 
   in mixed case particles were marked */
Float_t AliTkConeJet::getEtMarked(Float_t ptcut) const 
{
  Float_t et=0.;
  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->GetWeight()>-100) continue; //(-123)
      if(particle->Pt()<ptcut) continue;
      et+=particle->Pt();
    }
    delete part;
  }
  delete iter;
  return et;
}

Float_t AliTkConeJet::getEtMarkedFrac(Float_t ptcut) const 
{
  Float_t et=0.,etall=0.;
  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      etall+=particle->Pt();
      if(particle->GetWeight()>-100) continue; //(-123)
      et+=particle->Pt();
    }
    delete part;
  }
  delete iter;
  if(etall>0) return et/etall;
  return 1;
}

Float_t AliTkConeJet::getEtEM() const 
{
  Float_t EtEM = 0;
  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    EtEM += tower->getEtEM();
  }
  delete iter;
  return EtEM; 
}

TClonesArray *AliTkConeJet::getParticles(Float_t ptcut) const 
{
  TClonesArray *allParticles = new TClonesArray("TParticle",10000);
  Int_t nParticles = 0;
  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      new ((*allParticles)[nParticles]) TParticle(*particle);
      nParticles++;
    }
    delete part;
  }
  delete iter;
  return allParticles;
}

TClonesArray *AliTkConeJet::getChargedParticles(Float_t ptcut) const 
{
  TClonesArray *chargedParticles = new TClonesArray("TParticle",10000);
  Int_t nChargedParticles = 0;
  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getChargedParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      new ((*chargedParticles)[nChargedParticles]) TParticle(*particle);
      nChargedParticles++;
    }
    particles->Delete();
    delete particles;
    delete part;
  }
  delete iter;
  return chargedParticles;
}

TClonesArray *AliTkConeJet::getNeutralParticles(Float_t ptcut) const 
{
  TClonesArray *neutralParticles = new TClonesArray("TParticle",10000);
  Int_t nNeutralParticles = 0;
  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getNeutralParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      new ((*neutralParticles)[nNeutralParticles]) TParticle(*particle);
      nNeutralParticles++;
    }
    particles->Delete();
    delete particles;
    delete part;
  }
  delete iter;
  return neutralParticles;
}

Int_t AliTkConeJet::getNParticles() const 
{
  Int_t n = -1;
  TClonesArray *p = getParticles();
  n = p->GetEntries();
  p->Delete();
  delete p;
  return n;
}

Int_t AliTkConeJet::getNChargedParticles() const 
{
  Int_t n = -1;
  TClonesArray *p = getChargedParticles();
  n = p->GetEntries();
  p->Delete();
  delete p;
  return n;
}

Int_t AliTkConeJet::getNNeutralParticles() const 
{
  Int_t n = -1;
  TClonesArray *p = getNeutralParticles();
  n = p->GetEntries();
  p->Delete();
  delete p;
  return n;
}

Int_t AliTkConeJet::getNParticles(Float_t ptCut) const 
{
  Int_t n = 0;
  TClonesArray *particles = getParticles();
  TIterator *iter = particles->MakeIterator();
  TParticle *particle;
  while ((particle = (TParticle *)iter->Next()) != NULL) {
    if (particle->Pt() > ptCut) {
      n++;
    }
  }
  delete iter;
  particles->Delete();
  delete particles;
  return n;
}

Int_t AliTkConeJet::getNChargedParticles(Float_t ptCut) const 
{
  Int_t n = 0;
  TClonesArray *particles = getChargedParticles();
  TIterator *iter = particles->MakeIterator();
  TParticle *particle;
  while ((particle = (TParticle *)iter->Next()) != NULL) {
    if (particle->Pt() > ptCut) {
      n++;
    }
  }
  delete iter;
  particles->Delete();
  delete particles;
  return n;
}

Int_t AliTkConeJet::getNNeutralParticles(Float_t ptCut) const 
{
  Int_t n = 0;
  TClonesArray *particles = getNeutralParticles();
  TIterator *iter = particles->MakeIterator();
  TParticle *particle;
  while ((particle = (TParticle *)iter->Next()) != NULL) {
    if (particle->Pt() > ptCut) {
      n++;
    }
  }
  delete iter;
  particles->Delete();
  delete particles;
  return n;
}

void AliTkConeJet::getAxis(Float_t &x,Float_t &y,Float_t &z,Float_t ptcut) const
{
  x=0.;
  y=0.;
  z=0.;

  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      x+=particle->Px();
      y+=particle->Py();
      z+=particle->Pz();
    }
    delete part;
  }
  delete iter;
}

void AliTkConeJet::getChAxis(Float_t &x,Float_t &y,Float_t &z,Float_t ptcut) const
{
  x=0.;
  y=0.;
  z=0.;

  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getChargedParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      x+=particle->Px();
      y+=particle->Py();
      z+=particle->Pz();
    }
    particles->Delete();
    delete particles;
    delete part;
  }
  delete iter;
  Float_t length=TMath::Sqrt(x*x+y*y+z*z);
  if(length!=0){
    x/=length;
    y/=length;
    z/=length;
  }
}

TParticle* AliTkConeJet::getLeadingPart(Float_t ptcut) const
{
  TParticle* plead=new TParticle; 
  Float_t ptlead=0;

  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      if(particle->Pt()>ptlead){
	ptlead=particle->Pt();
	delete plead;
        plead=new TParticle(*particle);
      }
    }
    delete part;
  }
  delete iter;
  return plead;
}

TParticle* AliTkConeJet::getLeadingChPart(Float_t ptcut) const
{
  TParticle* plead=new TParticle; 
  Float_t ptlead=0;

  TIterator *iter = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)iter->Next()) != NULL) {
    TClonesArray *particles = tower->getChargedParticleList();
    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      if(particle->Pt()>ptlead){
	ptlead=particle->Pt();
	delete plead;
        plead=new TParticle(*particle);
      }
    }
    particles->Delete();
    delete particles;
    delete part;
  }
  delete iter;
  return plead;
}

void AliTkConeJet::calculateValues(Float_t ptcut)
{
  fPtCut=ptcut;
  if(fParts) delete fParts;
  fParts=new TClonesArray("TParticle",fNTowers*25);

  Float_t px=0.,py=0.,pz=0.;
  TParticle* plead=new TParticle; 
  Float_t ptlead=0;
  Int_t counter=0;
  TIterator *titer = fTowers->MakeIterator();
  AliTkTowerV2 *tower = NULL;
  while ((tower = (AliTkTowerV2 *)titer->Next()) != NULL) {
    TClonesArray *particles = tower->getParticleList();

    TParticle *particle = NULL;
    TIterator *part = particles->MakeIterator();
    while ((particle = (TParticle *)part->Next()) != NULL) {
      if(particle->Pt()<ptcut) continue;
      px+=particle->Px();
      py+=particle->Py();
      pz+=particle->Pz();
      new ((*fParts)[counter]) TParticle(*particle);
      counter++;
      if(particle->Pt()>ptlead){
	ptlead=particle->Pt();
	delete plead;
        plead=new TParticle(*particle);
      }
    }
    delete part;
  }
  delete titer;

  fCEt = TMath::Sqrt(px*px+py*py);
  Float_t p = TMath::Sqrt(px*px+py*py+pz*pz);
  fPLength = p;
  //Float_t theta = ((pz==0)||(p==0))?TMath::PiOver2():TMath::ACos(pz/p);
  fCPhi = TMath::Pi()+TMath::ATan2(-py,-px);
  //fCEta = (theta==0)?-100:TMath::Log(TMath::Tan(theta/2.));
  fCEta = (p==pz)?-100:0.5*TMath::Log((p+pz)/(p-pz));
  if(p){
    fXAxis= px/p;
    fYAxis= py/p;
    fZAxis= pz/p;
  }
  fPtLead=ptlead;
  if(fLeadPart) delete fLeadPart;
  fLeadPart=new TParticle(*plead);
  delete plead;
  fNParts=counter;
}

Int_t AliTkConeJet::Compare(const TObject *obj) const
{
  Double_t val=((AliTkConeJet*)obj)->getEt();

  if(fEt>val) return -1; //qsort is ascending
  else if (fEt<val) return 1;
  else return 0;
}

ostream& operator<<(ostream& s,const AliTkConeJet &p) 
{
  return s << "AliTkConeJet info: eta=" << p.getEta()
	   << " phi=" << p.getPhi() << " Et=" << p.getEt();
}

Float_t AliTkConeJet::Diff(const AliTkConeJet *jet1, const AliTkConeJet *jet2, Float_t &etdiff, Float_t &phidiff, Float_t &etadiff)
{
  Float_t ret=0;
  phidiff=TMath::Abs(jet1->getPhi()-jet2->getPhi());
  if(phidiff>TMath::Pi()) ret=TMath::TwoPi()-phidiff;
  etadiff=TMath::Abs(jet1->getEta()-jet2->getEta());
  ret=TMath::Sqrt(phidiff*phidiff+etadiff*etadiff);
  etdiff=jet1->getEt()-jet2->getEt();
  return ret;
}

ClassImp(AliTkConeJet)
