// $Id$

#include <Riostream.h>

#include <TClonesArray.h>
#include <TIterator.h>
#include <TParticle.h>
#include <TMath.h>
#include <TH2F.h>

#include "AliJFJetCalorimeterTrigger.h"
#include "AliJFJetCalorimeterTriggerResult.h"

ClassImp(AliJFJetCalorimeterTrigger)

AliJFJetCalorimeterTrigger::AliJFJetCalorimeterTrigger(Int_t n) : AliJFJetTrigger(n),
                                                            fRin(0.2),fRmid(0.7),fRout(1.0),
							    fPhiBins(25),fEtaBins(25)
{
  fAin=TMath::Pi()*fRin*fRin;
  fAout=TMath::Pi()*(fRout*fRout-fRin*fRin);
}

AliJFJetCalorimeterTrigger::~AliJFJetCalorimeterTrigger()
{
  Clean();
}

Int_t AliJFJetCalorimeterTrigger::Init(TClonesArray *particles)
{
  if(particles==NULL) return -1;

  TIterator *iter=particles->MakeIterator();
  TParticle *p;
  Int_t ret=0;

  while((p=(TParticle*)iter->Next()) != NULL){
    if(IsAcceptedParticle(p)){
      new((*fParticles)[ret]) TParticle(*p);
      ret++;
    }
  }
  delete iter;

  fSeedPlane=new TH2F("seedplane","seedplane",fEtaBins,fEtaMin,fEtaMax,fPhiBins,fPhiMin,fPhiMax);

  return ret;
}

Int_t AliJFJetCalorimeterTrigger::Run()
{
  TIterator *iter=fParticles->MakeIterator();
  TParticle *p;
  while((p=(TParticle*)iter->Next()) != NULL){
    fSeedPlane->Fill(p->Eta(),p->Phi(),p->Pt());
  }
  delete iter;

  fSeedPlane->Sumw2();
  Float_t thr=fSeedPlane->GetSumOfWeights()/fEtaBins/fPhiBins;

  //Find Maxima
  Int_t n=0;
  Float_t *mbin=new Float_t[Int_t(fEtaBins*fPhiBins)];
  Float_t *mphi=new Float_t[Int_t(fEtaBins*fPhiBins)];
  Float_t *meta=new Float_t[Int_t(fEtaBins*fPhiBins)];
  for(Int_t eta=1; eta<=fEtaBins; eta++) {
    for(Int_t phi=1; phi<=fPhiBins; phi++) {
      Int_t bin5=fSeedPlane->GetBin(eta,phi);

      Double_t val5=fSeedPlane->GetBinContent(bin5);
      Double_t val8=0,val2=0,val4=0,val6=0,val7=0,val9=0,val3=0,val1=0;

      Int_t phiplusone=phi+1;
      if(phiplusone>fPhiBins) phiplusone=0;
      Int_t phiminusone=phi-1;
      if(phiminusone<1) phiminusone=fPhiBins;

      val4=fSeedPlane->GetBinContent(fSeedPlane->GetBin(eta,phiminusone));
      val6=fSeedPlane->GetBinContent(fSeedPlane->GetBin(eta,phiplusone));

      if(eta>1){
	val2=fSeedPlane->GetBinContent(fSeedPlane->GetBin(eta-1,phi));
	val3=fSeedPlane->GetBinContent(fSeedPlane->GetBin(eta-1,phiplusone));
	val1=fSeedPlane->GetBinContent(fSeedPlane->GetBin(eta-1,phiminusone));
      }	
      if(eta<fEtaBins){
	val8=fSeedPlane->GetBinContent(fSeedPlane->GetBin(eta+1,phi));
	val7=fSeedPlane->GetBinContent(fSeedPlane->GetBin(eta+1,phiminusone));
	val9=fSeedPlane->GetBinContent(fSeedPlane->GetBin(eta+1,phiplusone));
      }

      if((val5>val1)&&(val5>val2)&&(val5>val3)&&(val5>val4)&&
	 (val5>val6)&&(val5>val7)&&(val5>val8)&&(val5>val9)&&
	 (val5>thr)){
	meta[n]=fSeedPlane->GetXaxis()->GetBinCenter(eta);
	mphi[n]=fSeedPlane->GetYaxis()->GetBinCenter(phi);
	mbin[n]=val5;
	n++;
      }
    }
  }

  Float_t *binweight=new Float_t[fNJetsMax];
  Float_t *phicenter=new Float_t[fNJetsMax];
  Float_t *etacenter=new Float_t[fNJetsMax];
  for(Int_t i=0;i<fNJetsMax;i++){
    Float_t tempmax=-1;
    Int_t tempind=-1;
    if(n==0) break;
    for(Int_t j=0;j<n;j++){
      if(tempmax<mbin[j]){
	tempmax=mbin[j];
	tempind=j;
      }
    }
    binweight[i]=tempmax;
    phicenter[i]=mphi[tempind];
    etacenter[i]=meta[tempind];
    mbin[tempind]=mbin[n]; //delete found maximum
    meta[tempind]=meta[n];
    mphi[tempind]=mphi[n];
    n--;
    fNJets++;
    //cout << i << " " << (Float_t)binweight[i] << " " << (Float_t)phicenter[i] << " " << (Float_t)etacenter[i] << endl;
  }
  delete[] mphi;
  delete[] meta;
  delete[] mbin;

  Float_t fRin2=fRin*fRin;
  Float_t fRmid2=fRmid*fRmid;
  Float_t fRout2=fRout*fRout;

  Float_t *inner=new Float_t[fNJets];
  Float_t *outer=new Float_t[fNJets];
  for(Int_t i=0;i<fNJets;i++){
    inner[i]=0;
    outer[i]=0;
  }
  iter=fParticles->MakeIterator();
  while((p=(TParticle*)iter->Next()) != NULL){
    Float_t pe=p->Eta();
    Float_t pp=p->Phi();
    Float_t pt=p->Pt();

    for(Int_t i=0;i<fNJets;i++){
      Float_t a=pe-etacenter[i];
      Float_t b=pp-phicenter[i];
      Float_t rad=a*a+b*b;
      if(rad<fRin2) inner[i]+=pt;
      else if((rad>fRmid2)&&(rad<fRout2)) outer[i]+=pt;
    }
  }
  delete iter;
  
  for(Int_t i=0;i<fNJets;i++){
    //inner[i]/=fAin;
    outer[i]=outer[i]*fAin/fAout;
    fJets.AddAt(new AliJFJetCalorimeterTriggerResult(fPtMin,fPtMax,fEtaMin,fEtaMax,
                                                  fPhiMin,fPhiMax,fNeutral,fCharged,fEM,
						  phicenter[i],etacenter[i],binweight[i],thr,
						  inner[i],outer[i],fRin,fRmid,fRout,
		                                  fPhiBins,fEtaBins),i);
    cout << i << " " << inner[i] << " " << outer[i] << endl;
  }

  delete[] binweight;
  delete[] phicenter;
  delete[] etacenter;
  delete[] inner;
  delete[] outer;
  return fNJets;
}

void AliJFJetCalorimeterTrigger::Clean()
{
  if(fSeedPlane) delete fSeedPlane;
  AliJFJetFinder::Clean();
}

void AliJFJetCalorimeterTrigger::SetRadius(Float_t in, Float_t mid, Float_t out)
{
  fRin=in;
  fRmid=mid;
  fRout=out;
  fAin=TMath::Pi()*fRin*fRin;
  fAout=TMath::Pi()*(fRout*fRout-fRmid*fRmid);
  cout << fAin << " " << fAout << endl;
}

void AliJFJetCalorimeterTrigger::SetBins(Int_t phis, Int_t etas)
{
  fPhiBins=phis;
  fEtaBins=etas;
}
