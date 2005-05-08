// $Id$

#include <Riostream.h>
#include <vector>

#include <TParticle.h>
#include <TClonesArray.h>
#include <TIterator.h>
#include <TMath.h>

#include "AliJFKtJet.h"
#include "AliJFPreCluster.h"
#include "AliJFCluster.h"

ClassImp(AliJFKtJet)

AliJFKtJet::AliJFKtJet(Int_t n) : AliJFJet(n)

{
}

AliJFKtJet::~AliJFKtJet()
{
  Clean();
}

void AliJFKtJet::AddJet(AliJFCluster *c)
{
  if(!c) return;
  if(!c->IsJet()){
    cerr << "AliJFKtJet Error: Something is strange here, supposed to add a jet!" << endl;
    return ;
  }

  vector<AliJFPreCluster*> clist=*c->GetClusterList();
  for(vector<AliJFPreCluster*>::iterator i=clist.begin();i!=clist.end();i++){

    //Get particles inside PreCluster
    TClonesArray *ps=(*i)->GetParticles();

    //Copy Particles inside PreCluster
    TIterator *iter=ps->MakeIterator();
    TParticle *p;
    while((p=(TParticle*)iter->Next()) != NULL){
      new (fParticles[fN]) TParticle(*p);
      
      fN++;
    }
  }

  fIsUpdated=kFALSE;
}

void AliJFKtJet::Update()
{
  if(fIsUpdated) return;

  fIsUpdated=kTRUE;

  fPhi=0;
  fEta=0;
  fY=0;
  fPt=0;
  fPx=0;
  fPy=0;
  fPz=0;
  fE=0;
  fE_=0;
  fPtSum=0;
  fPhiSum=0;
  fEtaSum=0;

  fPhiC=0;
  fEtaC=0;
  fYC=0;
  fPtC=0;
  fPxC=0;
  fPyC=0;
  fPzC=0;
  fEC=0;
  fE_C=0;
  fPtSumC=0;
  fPhiSumC=0;
  fEtaSumC=0;

  fPhiN=0;
  fEtaN=0;
  fYN=0;
  fPtN=0;
  fPxN=0;
  fPyN=0;
  fPzN=0;
  fEN=0;
  fE_N=0;
  fPtSumN=0;
  fPhiSumN=0;
  fEtaSumN=0;

  fPhiEM=0;
  fEtaEM=0;
  fYEM=0;
  fPtEM=0;
  fPxEM=0;
  fPyEM=0;
  fPzEM=0;
  fEEM=0;
  fE_EM=0;
  fPtSumEM=0;
  fPhiSumEM=0;
  fEtaSumEM=0;

  fNCharged=0;
  fNNeutral=0;
  fNEM=0;

  Float_t fMaxParticlePt=0;
  Float_t fMaxParticlePtC=0;
  Float_t fMaxParticlePtN=0;
  Float_t fMaxParticlePtEM=0;

  //loop over particles in jet
  Int_t N=0;
  TParticle *p;
  TIterator *iter=fParticles.MakeIterator();
  while((p=(TParticle*)iter->Next()) != NULL){
    N++; 

    Float_t px=p->Px();
    Float_t py=p->Py();
    Float_t pz=p->Pz();
    Float_t E=p->Energy();
    Float_t E_=TMath::Sqrt(px*px+py*py+pz*pz); //massless particles
    Float_t pt=p->Pt();
    Float_t phi=p->Phi();
    Float_t eta=p->Eta();

    fPx+=px;
    fPy+=py;
    fPz+=pz;
    fE+=E;
    fE_+=E_;

    fPtSum+=pt;
    fPhiSum+=phi*pt;
    fEtaSum+=eta*pt;

    if(pt>fMaxParticlePt){
      fMaxParticlePt=pt;
      fMaxParticle=*p;
    }

    Int_t pcode=p->GetPdgCode();
    if((pcode==11)||(pcode==-11)||(pcode==22)){ //EM Particle;
      fNEM++;

      fPxEM+=px;
      fPyEM+=py;
      fPzEM+=pz;
      fEEM+=E;
      fE_EM+=E_;

      fPtSumEM+=pt;
      fPhiSumEM+=phi*pt;
      fEtaSumEM+=eta*pt;

      if(pt>fMaxParticlePtEM){
	fMaxParticlePtEM=pt;
	fMaxParticleEM=*p;
      }
    } else if(p->GetPDG()->Charge()) {//Charged Particle
      fNCharged++;

      fPxC+=px;
      fPyC+=py;
      fPzC+=pz;
      fEC+=E;
      fE_C+=E_;
      fPtSumC+=pt;
      fPhiSumC+=phi*pt;
      fEtaSumC+=eta*pt;

      if(pt>fMaxParticlePtC){
	fMaxParticlePtC=pt;
	fMaxParticleC=*p;
      }
    }  else { //Neutral Particle
      fNNeutral++;

      fPxN+=px;
      fPyN+=py;
      fPzN+=pz;
      fEN+=E;
      fE_N+=E_;
      fPtSumN+=pt;
      fPhiSumN+=phi*pt;
      fEtaSumN+=eta*pt;

      if(pt>fMaxParticlePtN){
	fMaxParticlePtN=pt;
	fMaxParticleN=*p;
      }
    }    
  } //end loop 

  fPt=TMath::Sqrt(fPx*fPx+fPy*fPy);
  //fPhi=TMath::ATan(fPy/fPx);
  fPhi=TMath::Pi()+TMath::ATan2(-fPy,-fPx);
  fY=0.5*TMath::Log((fE+fPz)/(fE-fPz));
  fEta=0.5*TMath::Log((fE_+fPz)/(fE_-fPz));

  fPtC=TMath::Sqrt(fPxC*fPxC+fPyC*fPyC);
  //fPhiC=TMath::ATan(fPyC/fPxC);
  fPhiC=TMath::Pi()+TMath::ATan2(-fPyC,-fPxC);
  fYC=0.5*TMath::Log((fEC+fPzC)/(fEC-fPzC));
  fEtaC=0.5*TMath::Log((fE_C+fPzC)/(fE_C-fPzC));

  fPtN=TMath::Sqrt(fPxN*fPxN+fPyN*fPyN);
  //fPhiN=TMath::ATan(fPyN/fPxN);
  fPhiN=TMath::Pi()+TMath::ATan2(-fPyN,-fPxN);
  fYN=0.5*TMath::Log((fEN+fPz)/(fEN-fPz));
  fEtaN=0.5*TMath::Log((fE_N+fPz)/(fE_N-fPz));

  fPtEM=TMath::Sqrt(fPxEM*fPxEM+fPyEM*fPyEM);
  //fPhiEM=TMath::ATan(fPyEM/fPxEM);
  fPhiEM=TMath::Pi()+TMath::ATan2(-fPyEM,-fPxEM);
  fYEM=0.5*TMath::Log((fEEM+fPz)/(fEEM-fPz));
  fEtaEM=0.5*TMath::Log((fE_EM+fPz)/(fE_EM-fPz));

  fPhiSum/=fPtSum;
  fEtaSum/=fPtSum;
  fPhiSumEM/=fPtSumEM;
  fEtaSumEM/=fPtSumEM;
  fPhiSumC/=fPtSumC;
  fEtaSumC/=fPtSumC;
  fPhiSumN/=fPtSumN;
  fEtaSumN/=fPtSumN;

  //Do some simple checks
  if(N!=fN){
    cerr << "Fatal Error: Particles in Jets don't match!" << endl;
    exit(1);
  }
  if(fNCharged+fNNeutral+fNEM!=fN){
    cerr << "Fatal Error: Particles in Jets don't sum up!" << endl;
    exit(1);
  }
}

/*
void AliJFKtJet::Debug()
{
}

void AliJFKtJet::Clean()
{
  ::Clean();
}
*/
