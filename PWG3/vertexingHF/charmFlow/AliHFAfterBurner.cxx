/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/////////////////////////////////////////////////////////////
//
// AliHFAfterBurner introduces a v2 modulation in MC for 
// the D mesons v2 analysis with event plane method
// Authors: Giacomo Ortona, ortona@to.infn.it
// 
/////////////////////////////////////////////////////////////

/* $Id$ */

#include <TDatabasePDG.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <AliAODEvent.h>
#include <AliAODMCParticle.h>
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliHFAfterBurner.h"

ClassImp(AliHFAfterBurner)

//________________________________________________
AliHFAfterBurner::AliHFAfterBurner():
  fSigv2(0),
  fBkgv2(0),
  fUseNewton(kTRUE),
  fPrecisionNewton(0),
  fDecChannel(0),
  fSignal(0)
{
  //empty constructor
}
//________________________________________________
AliHFAfterBurner::AliHFAfterBurner(Int_t decChannel):
  fSigv2(0.1),
  fBkgv2(0.2),
  fUseNewton(kTRUE),
  fPrecisionNewton(0.0005),
  fDecChannel(decChannel),
  fSignal(0)
{
  //default constructor
}
//______________________________________________________________________________
AliHFAfterBurner::AliHFAfterBurner(const AliHFAfterBurner &source):
  TObject(source),
  fSigv2(source.fSigv2),
  fBkgv2(source.fBkgv2),
  fUseNewton(source.fUseNewton),
  fPrecisionNewton(source.fPrecisionNewton),
  fDecChannel(source.fDecChannel),
  fSignal(source.fSignal)
{
  //copy constructor
}
//______________________________________________________________________________
AliHFAfterBurner &AliHFAfterBurner::operator=(const AliHFAfterBurner &source)
{
  //assignment operator
  if(&source == this) return *this;

  TObject::operator=(source);
  fSigv2 = source.fSigv2;
  fBkgv2=source.fBkgv2;
  fUseNewton=source.fUseNewton;
  fPrecisionNewton=source.fPrecisionNewton;
  fDecChannel=source.fDecChannel;
  fSignal=source.fSignal;
  return *this;
}
//______________________________________________________________________________
AliHFAfterBurner::~AliHFAfterBurner(){

}
//______________________________________________________________________________
Double_t AliHFAfterBurner::GetNewAngle(AliAODRecoDecayHF *d,TClonesArray *mcArray){
  Int_t lab=-1;
  Int_t *pdgdaughters;
  Int_t pdgmother=0;
  Int_t nProngs=0;
  switch(fDecChannel){
  case 0:
    //D+
    nProngs=3;
    pdgmother=411;
    pdgdaughters=new Int_t[3];
    pdgdaughters[0]=211;//pi
    pdgdaughters[1]=321;//K
    pdgdaughters[2]=211;//pi
    break;
  case 1:
    //D0
    nProngs=2;
    pdgmother=421;
    pdgdaughters=new Int_t[2];
    pdgdaughters[0]=211;//pi 
    pdgdaughters[1]=321;//K
    break;
  case 2:
    //D*
    nProngs=3;
    pdgmother=413;
    pdgdaughters=new Int_t[3];
    pdgdaughters[1]=211;//pi
    pdgdaughters[0]=321;//K
    pdgdaughters[2]=211;//pi (soft?)
    break;
  }
  lab = d->MatchToMC(pdgmother,mcArray,nProngs,pdgdaughters);
  Double_t phi=-999.;
  if(lab>=0){
    fSignal=kTRUE;
    phi=GetPhi(d->Phi(),fSigv2);
  }
  else {//phi=NewtonMethodv2(phi,fBkgv2,eventplane);
    //background
    fSignal=kFALSE;
    //const Int_t nProngs=d->GetNProngs();
    Float_t phidau[nProngs];
    for(Int_t ipr=0;ipr<nProngs;ipr++){
      phidau[ipr]=(Float_t)d->PhiProng(ipr);
      //AliAODTrack *trk = (AliAODTrack*)d->GetDaughter(ipr);
      Int_t labdau=(Int_t)d->GetProngID(ipr);
      if(labdau<0)continue;
      AliAODMCParticle *mcpart= (AliAODMCParticle*)mcArray->At(labdau);
      if(!mcpart)continue;
      Int_t laborig=TMath::Abs(CheckOrigin(mcpart,mcArray));
      if(laborig>=0){//charm
	mcpart= (AliAODMCParticle*)mcArray->At(laborig);
	if(mcpart)phidau[ipr]=GetPhi(mcpart->Phi(),fSigv2);
      }else{//not charm
	phidau[ipr]=GetPhi(phidau[ipr],fBkgv2);
      }
    }
    Float_t py=0,px=0;
    for(Int_t ipr=0;ipr<nProngs;ipr++){
      py+=d->PtProng(ipr)*TMath::Sin(phidau[ipr]);
      px+=d->PtProng(ipr)*TMath::Cos(phidau[ipr]);
    }
    phi=TMath::Pi()+TMath::ATan2(-py,-px);
  }
  return GetPhi02Pi(phi);
}
//______________________________________________________________________________
Double_t AliHFAfterBurner::GetPhi(Double_t phi,Float_t v2){
  Double_t evplane = fEventPlane;
  if(fUseNewton){
    return NewtonMethodv2(phi,v2);
  }else{
    return phi-v2*TMath::Sin(2*(phi-evplane));
  }
}
//______________________________________________________________________________
Double_t AliHFAfterBurner::NewtonMethodv2(Double_t phi,Double_t v2,Double_t phi0){
  Double_t eventplane = fEventPlane;
  Double_t phi1 = phi-(phi+v2*TMath::Sin(2.*(phi-eventplane))-phi0)/(1.+2.*v2*TMath::Cos(2.*(phi-eventplane)));
  if(TMath::Abs(phi/phi1-1.)<fPrecisionNewton){
    return phi1;
  }else {
    return NewtonMethodv2(phi1,v2,phi0);
  }
}
//______________________________________________________________________________
void AliHFAfterBurner::SetMCv2(Float_t v2sig,Float_t v2bkg){
  if(v2sig>=0)fSigv2=v2sig;
  if(v2bkg>=0)fBkgv2=v2bkg;
}
//______________________________________________________________________________
Int_t AliHFAfterBurner::CheckOrigin(const AliAODMCParticle* mcPart,TClonesArray *arrayMC)const{
  
  //
  // checking whether the mother of the particle is a charm
  //
  Float_t charmmother=kFALSE;
  Int_t pdgGranma = 0;
  Int_t mother = -1;
  mother = mcPart->GetMother();
  Int_t istep = 0;
  while (mother >=0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if(!mcGranma) break;
    pdgGranma = mcGranma->GetPdgCode();
    Int_t abspdgGranma = TMath::Abs(pdgGranma);
    if ((abspdgGranma > 400 && abspdgGranma < 500) || (abspdgGranma > 4000 && abspdgGranma < 5000)) {
      charmmother=kTRUE;
    }
    mother = mcGranma->GetMother();
  }
  if(!charmmother)mother=-1;
  return mother;
}
//________________________________________________________________________
Float_t AliHFAfterBurner::GetPhi02Pi(Float_t phi){
  Float_t result=phi;
  while(result<0){
    result=result+2.*TMath::Pi();
  }
  while(result>TMath::Pi()*2.){
    result=result-2.*TMath::Pi();
  }
  return result;
}
//________________________________________________________________________
void AliHFAfterBurner::SetDecChannel(Int_t decch){
  if(decch>2){
    AliWarning("Invalid decay channel");
    return;
  }
  fDecChannel=decch;
}
//________________________________________________________________________
void AliHFAfterBurner::SetEventPlaneMethod(Int_t method){
  //Only Random EP supported by now, feature to generate EP from histos to be added
  fMethodEP=method;
  return;
}
//________________________________________________________________________
void AliHFAfterBurner::SetEventPlane(){
  //Only Random EP supported by now, feature to generate EP from histos to be added
  TRandom3 *g = new TRandom3(0);
  fEventPlane=g->Rndm()*TMath::Pi();
  delete g;g=0x0;
  return;
}
