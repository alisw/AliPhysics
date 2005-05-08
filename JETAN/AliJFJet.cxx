// $Id$

#include <Riostream.h>

#include "AliJFJet.h"

ClassImp(AliJFJet)

AliJFJet::AliJFJet(Int_t n) : TObject(),
		 fNJet(0),fN(0),fNCharged(0),fNNeutral(0),fNEM(0),
		 fPhi(0),fEta(0),fY(0),fPt(0),
		 fPx(0),fPy(0),fPz(0),fE(0),
		 fPtSum(0),fPhiSum(0),fEtaSum(0),
		 fPhiC(0),fEtaC(0),fYC(0),fPtC(0),
		 fPxC(0),fPyC(0),fPzC(0),fEC(0),
		 fPtSumC(0),fPhiSumC(0),fEtaSumC(0),
		 fPhiN(0),fEtaN(0),fYN(0),fPtN(0),
		 fPxN(0),fPyN(0),fPzN(0),fEN(0),
		 fPtSumN(0),fPhiSumN(0),fEtaSumN(0),
		 fPhiEM(0),fEtaEM(0),fYEM(0),fPtEM(0),
		 fPxEM(0),fPyEM(0),fPzEM(0),fEEM(0),
		 fPtSumEM(0),fPhiSumEM(0),fEtaSumEM(0),
                 fMaxParticle(),fMaxParticleC(),fMaxParticleN(),fMaxParticleEM(),
                 fParticles("TParticle",n),fIsUpdated(kTRUE)
{
}

AliJFJet::~AliJFJet()
{
}

Int_t AliJFJet::Compare(const TObject *obj) const
{
  Double_t val=((AliJFJet*)obj)->GetPt();

  if(fPt>val) return 1;
  else if (fPt<val) return -1;
  else return 0;
}

void AliJFJet::Clean()
{
  fParticles.Delete();
}

void AliJFJet::Debug()
{
  cout << "Jet: " << fNJet << endl;

  TIterator *iter=fParticles.MakeIterator();
  TParticle *p;
  Int_t i=0;
  while((p=(TParticle*)iter->Next()) != NULL){
    cout << i++ << ": " << p->Energy() << endl;
  }
}
