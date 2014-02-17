/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

#include <TMath.h>
#include <TRandom.h>
#include <TProfile.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF.h"
#include "AliVertexingHFUtils.h"

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class with functions useful for different D2H analyses        //
// - event plane resolution                                      //
// - <pt> calculation with side band subtraction                 //
// - tracklet multiplicity calculation                            //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

ClassImp(AliVertexingHFUtils)

//______________________________________________________________________
AliVertexingHFUtils::AliVertexingHFUtils():TObject(),
  fK(1),
  fSubRes(1.),
  fMinEtaForTracklets(-1.),
  fMaxEtaForTracklets(1.)
{
  // Default contructor
}


//______________________________________________________________________
AliVertexingHFUtils::AliVertexingHFUtils(Int_t k):
  TObject(),
  fK(k),
  fSubRes(1.),
  fMinEtaForTracklets(-1.),
  fMaxEtaForTracklets(1.)
{
  // Standard constructor
}


//______________________________________________________________________
void AliVertexingHFUtils::ComputeSignificance(Double_t signal, Double_t  errsignal, Double_t  background, Double_t  errbackground, Double_t &significance,Double_t &errsignificance){
  // calculate significance from S, B and errors


  Double_t errSigSq=errsignal*errsignal;
  Double_t errBkgSq=errbackground*errbackground;
  Double_t sigPlusBkg=signal+background;
  if(sigPlusBkg>0. && signal>0.){
    significance =  signal/TMath::Sqrt(signal+background);
    errsignificance = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
  }else{
    significance=0.;
    errsignificance=0.;
  }
  return;

}
//______________________________________________________________________
Double_t AliVertexingHFUtils::Pol(Double_t x, Int_t k){
  // compute chi from polynomial approximation
  Double_t c[5];
  if(k==1){ 
    c[0]=0.626657; c[1]=0.; c[2]=-0.09694; c[3]=0.02754; c[4]=-0.002283;
  }
  else if(k==2){
    c[0]=0.; c[1]=0.25; c[2]=-0.011414; c[3]=-0.034726; c[4]=0.006815;
  }
  else return -1;
  return c[0]*x+c[1]*x*x+c[2]*x*x*x+c[3]*x*x*x*x+c[4]*x*x*x*x*x;
}

//______________________________________________________________________
Double_t AliVertexingHFUtils:: ResolK1(Double_t x){
  return TMath::Sqrt(TMath::Pi()/8)*x*TMath::Exp(-x*x/4)*(TMath::BesselI0(x*x/4)+TMath::BesselI1(x*x/4));
}


//______________________________________________________________________
Double_t AliVertexingHFUtils::FindChi(Double_t res,  Int_t k){
  // compute chi variable (=v2*sqrt(N)) from external values

  Double_t x1=0;
  Double_t x2=15;
  Double_t y1,y2;
  if(k==1){
    y1=ResolK1(x1)-res;
    y2=ResolK1(x2)-res;
  }
  else if(k==2){
    y1=Pol(x1,2)-res;
    y2=Pol(x2,2)-res;
  }
  else return -1;

  if(y1*y2>0) return -1;
  if(y1==0) return y1;
  if(y2==0) return y2;
  Double_t xmed,ymed;
  Int_t jiter=0;
  while((x2-x1)>0.0001){
    xmed=0.5*(x1+x2);
    if(k==1){
      y1=ResolK1(x1)-res;
      y2=ResolK1(x2)-res;
      ymed=ResolK1(xmed)-res;
    }
    else if(k==2){
      y1=Pol(x1,2)-res;
      y2=Pol(x2,2)-res;
      ymed=Pol(xmed,2)-res;
    }
    else return -1;
    if((y1*ymed)<0) x2=xmed;
    if((y2*ymed)<0)x1=xmed;
    if(ymed==0) return xmed;
    jiter++;
  }
  return 0.5*(x1+x2);
}

//______________________________________________________________________
Double_t AliVertexingHFUtils::GetFullEvResol(Double_t resSub, Int_t k){
  // computes event plane resolution starting from sub event resolution
  Double_t chisub=FindChi(resSub,k);
  Double_t chifull=chisub*TMath::Sqrt(2);
  if(k==1) return ResolK1(chifull);
  else if(k==2) return Pol(chifull,2);
  else{
    printf("k should be <=2\n");
    return 1.;
  }
}

//______________________________________________________________________
Double_t AliVertexingHFUtils::GetFullEvResol(const TH1F* hSubEvCorr, Int_t k){
  // computes event plane resolution starting from sub event correlation histogram
  if(!hSubEvCorr) return 1.;
  Double_t resSub=GetSubEvResol(hSubEvCorr);
  return GetFullEvResol(resSub,k);
}
//______________________________________________________________________
Double_t AliVertexingHFUtils::GetFullEvResolLowLim(const TH1F* hSubEvCorr, Int_t k){
  // computes low limit event plane resolution starting from sub event correlation histogram
  if(!hSubEvCorr) return 1.;
  Double_t resSub=GetSubEvResolLowLim(hSubEvCorr);
  printf("%f\n",resSub);
  return GetFullEvResol(resSub,k);  
}
//______________________________________________________________________
Double_t AliVertexingHFUtils::GetFullEvResolHighLim(const TH1F* hSubEvCorr, Int_t k){
  // computes high limit event plane resolution starting from sub event correlation histogram
  if(!hSubEvCorr) return 1.;
  Double_t resSub=GetSubEvResolHighLim(hSubEvCorr);
  printf("%f\n",resSub);
  return GetFullEvResol(resSub,k);  
}
//______________________________________________________________________
Int_t AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(AliAODEvent* ev, Double_t mineta, Double_t maxeta){
  // counts tracklets in given eta range
  AliAODTracklets* tracklets=ev->GetTracklets();
  Int_t nTr=tracklets->GetNumberOfTracklets();
  Int_t count=0;
  for(Int_t iTr=0; iTr<nTr; iTr++){
    Double_t theta=tracklets->GetTheta(iTr);
    Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
    if(eta>mineta && eta<maxeta) count++;
  }
  return count;
}
//______________________________________________________________________
Int_t AliVertexingHFUtils::GetGeneratedMultiplicityInEtaRange(TClonesArray* arrayMC, Double_t mineta, Double_t maxeta){
  // counts generated particles in fgiven eta range

  Int_t nChargedMC=0;
  for(Int_t i=0;i<arrayMC->GetEntriesFast();i++){
    AliAODMCParticle *part=(AliAODMCParticle*)arrayMC->UncheckedAt(i);
    Int_t charge = part->Charge();
    Double_t eta = part->Eta();
    if(charge!=0 && eta>mineta && eta<maxeta) nChargedMC++;
  } 
  return nChargedMC;
}
//______________________________________________________________________
Int_t AliVertexingHFUtils::GetGeneratedPrimariesInEtaRange(TClonesArray* arrayMC, Double_t mineta, Double_t maxeta){
  // counts generated primary particles in given eta range

  Int_t nChargedMC=0;
  for(Int_t i=0;i<arrayMC->GetEntriesFast();i++){
    AliAODMCParticle *part=(AliAODMCParticle*)arrayMC->UncheckedAt(i);
    Int_t charge = part->Charge();
    Double_t eta = part->Eta();
    if(charge!=0 && eta>mineta && eta<maxeta){
      if(part->IsPrimary())nChargedMC++;
    } 
  }  
  return nChargedMC;
}
//______________________________________________________________________
Int_t AliVertexingHFUtils::GetGeneratedPhysicalPrimariesInEtaRange(TClonesArray* arrayMC, Double_t mineta, Double_t maxeta){
  // counts generated primary particles in given eta range

  Int_t nChargedMC=0;
  for(Int_t i=0;i<arrayMC->GetEntriesFast();i++){
    AliAODMCParticle *part=(AliAODMCParticle*)arrayMC->UncheckedAt(i);
    Int_t charge = part->Charge();
    Double_t eta = part->Eta();
    if(charge!=0 && eta>mineta && eta<maxeta){
      if(part->IsPhysicalPrimary())nChargedMC++;
    } 
  }
  return nChargedMC;
}
//______________________________________________________________________
void AliVertexingHFUtils::AveragePt(Float_t& averagePt, Float_t& errorPt,Float_t ptmin,Float_t ptmax, TH2F* hMassD, Float_t massFromFit, Float_t sigmaFromFit, TF1* funcB2, Float_t sigmaRangeForSig,Float_t sigmaRangeForBkg, Float_t minMass, Float_t maxMass, Int_t rebin){

  // Compute <pt> from 2D histogram M vs pt

  //Make 2D histos in the desired pt range
  Int_t start=hMassD->FindBin(ptmin);
  Int_t end=hMassD->FindBin(ptmax)-1;
  const Int_t nx=end-start;
  TH2F *hMassDpt=new TH2F("hptmass","hptmass",nx,ptmin,ptmax,hMassD->GetNbinsY(),hMassD->GetYaxis()->GetBinLowEdge(1),hMassD->GetYaxis()->GetBinLowEdge(hMassD->GetNbinsY())+hMassD->GetYaxis()->GetBinWidth(hMassD->GetNbinsY()));
  for(Int_t ix=start;ix<end;ix++){
    for(Int_t iy=1;iy<=hMassD->GetNbinsY();iy++){
      hMassDpt->SetBinContent(ix-start+1,iy,hMassD->GetBinContent(ix,iy));
      hMassDpt->SetBinError(ix-start+1,iy,hMassD->GetBinError(ix,iy));
    }
  }

  Double_t minMassSig=massFromFit-sigmaRangeForSig*sigmaFromFit;
  Double_t maxMassSig=massFromFit+sigmaRangeForSig*sigmaFromFit;
  Int_t minBinSig=hMassD->GetYaxis()->FindBin(minMassSig);
  Int_t maxBinSig=hMassD->GetYaxis()->FindBin(maxMassSig);
  Double_t minMassSigBin=hMassD->GetYaxis()->GetBinLowEdge(minBinSig);
  Double_t maxMassSigBin=hMassD->GetYaxis()->GetBinLowEdge(maxBinSig)+hMassD->GetYaxis()->GetBinWidth(maxBinSig);
  //  printf("Signal Fit Limits = %f %f\n",minMassSigBin,maxMassSigBin);

  Double_t maxMassBkgLow=massFromFit-sigmaRangeForBkg*sigmaFromFit;
  Int_t minBinBkgLow=TMath::Max(hMassD->GetYaxis()->FindBin(minMass),2);
  Int_t maxBinBkgLow=hMassD->GetYaxis()->FindBin(maxMassBkgLow);
  Double_t minMassBkgLowBin=hMassD->GetYaxis()->GetBinLowEdge(minBinBkgLow);
  Double_t maxMassBkgLowBin=hMassD->GetYaxis()->GetBinLowEdge(maxBinBkgLow)+hMassD->GetYaxis()->GetBinWidth(maxBinBkgLow);
  Double_t minMassBkgHi=massFromFit+sigmaRangeForBkg*sigmaFromFit;
  Int_t minBinBkgHi=hMassD->GetYaxis()->FindBin(minMassBkgHi);
  Int_t maxBinBkgHi=TMath::Min(hMassD->GetYaxis()->FindBin(maxMass),hMassD->GetNbinsY()-1);
  Double_t minMassBkgHiBin=hMassD->GetYaxis()->GetBinLowEdge(minBinBkgHi);
  Double_t maxMassBkgHiBin=hMassD->GetYaxis()->GetBinLowEdge(maxBinBkgHi)+hMassD->GetYaxis()->GetBinWidth(maxBinBkgHi);
  //  printf("BKG Fit Limits = %f %f  && %f %f\n",minMassBkgLowBin,maxMassBkgLowBin,minMassBkgHiBin,maxMassBkgHiBin);

  Double_t bkgSig=funcB2->Integral(minMassSigBin,maxMassSigBin);
  Double_t bkgLow=funcB2->Integral(minMassBkgLowBin,maxMassBkgLowBin);
  Double_t bkgHi=funcB2->Integral(minMassBkgHiBin,maxMassBkgHiBin);
  //  printf("Background integrals = %f %f %f\n",bkgLow,bkgSig,bkgHi);

  TH1F* hMptBkgLo=(TH1F*)hMassDpt->ProjectionX("hPtBkgLoBin",minBinBkgLow,maxBinBkgLow);
  TH1F* hMptBkgHi=(TH1F*)hMassDpt->ProjectionX("hPtBkgHiBin",minBinBkgHi,maxBinBkgHi);
  TH1F* hMptSigReg=(TH1F*)hMassDpt->ProjectionX("hCPtBkgSigBin",minBinSig,maxBinSig);

  hMptBkgLo->Rebin(rebin);
  hMptBkgHi->Rebin(rebin);
  hMptSigReg->Rebin(rebin);

  hMptBkgLo->Sumw2();
  hMptBkgHi->Sumw2();
  TH1F* hMptBkgLoScal=(TH1F*)hMptBkgLo->Clone("hPtBkgLoScalBin");
  hMptBkgLoScal->Scale(bkgSig/bkgLow);
  TH1F* hMptBkgHiScal=(TH1F*)hMptBkgHi->Clone("hPtBkgHiScalBin");
  hMptBkgHiScal->Scale(bkgSig/bkgHi);

  TH1F* hMptBkgAver=0x0;
  hMptBkgAver=(TH1F*)hMptBkgLoScal->Clone("hPtBkgAverBin");
  hMptBkgAver->Add(hMptBkgHiScal);
  hMptBkgAver->Scale(0.5);
  TH1F* hMptSig=(TH1F*)hMptSigReg->Clone("hCPtSigBin");
  hMptSig->Add(hMptBkgAver,-1.);   
 
  averagePt = hMptSig->GetMean();
  errorPt = hMptSig->GetMeanError();

  delete hMptBkgLo;
  delete hMptBkgHi;
  delete hMptSigReg;
  delete hMptBkgLoScal;
  delete hMptBkgHiScal;
  delete hMptBkgAver;
  delete hMassDpt;
  delete hMptSig;

}
//____________________________________________________________________________
Double_t AliVertexingHFUtils::GetTrueImpactParameterDzero(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partD) {
  // true impact parameter calculation for Dzero

  if(!partD || !arrayMC || !mcHeader) return 99999.;
  Int_t code=partD->GetPdgCode();
  if(TMath::Abs(code)!=421) return 99999.;

  Double_t vtxTrue[3];
  mcHeader->GetVertex(vtxTrue);
  Double_t origD[3];
  partD->XvYvZv(origD);
  Short_t charge=partD->Charge();
  Double_t pXdauTrue[2],pYdauTrue[2],pZdauTrue[2];
  for(Int_t iDau=0; iDau<2; iDau++){
    pXdauTrue[iDau]=0.;
    pYdauTrue[iDau]=0.;
    pZdauTrue[iDau]=0.;
  }

  Int_t nDau=partD->GetNDaughters();
  Int_t labelFirstDau = partD->GetDaughter(0); 
  if(nDau==2){
    for(Int_t iDau=0; iDau<2; iDau++){
      Int_t ind = labelFirstDau+iDau;
      AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(ind));
      if(!part){
	printf("Daughter particle not found in MC array");
	return 99999.;
      } 
      pXdauTrue[iDau]=part->Px();
      pYdauTrue[iDau]=part->Py();
      pZdauTrue[iDau]=part->Pz();
    }
  }else{
    return 99999.;
  }

  Double_t d0dummy[2]={0.,0.};
  AliAODRecoDecayHF aodDvsMC(vtxTrue,origD,2,charge,pXdauTrue,pYdauTrue,pZdauTrue,d0dummy);
  return aodDvsMC.ImpParXY();

}
//____________________________________________________________________________
Double_t AliVertexingHFUtils::GetTrueImpactParameterDplus(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partD) {
  // true impact parameter calculation for Dplus

  if(!partD || !arrayMC || !mcHeader) return 99999.;
  Int_t code=partD->GetPdgCode();
  if(TMath::Abs(code)!=411) return 99999.;

  Double_t vtxTrue[3];
  mcHeader->GetVertex(vtxTrue);
  Double_t origD[3];
  partD->XvYvZv(origD);
  Short_t charge=partD->Charge();
  Double_t pXdauTrue[3],pYdauTrue[3],pZdauTrue[3];
  for(Int_t iDau=0; iDau<3; iDau++){
    pXdauTrue[iDau]=0.;
    pYdauTrue[iDau]=0.;
    pZdauTrue[iDau]=0.;
  }

  Int_t nDau=partD->GetNDaughters();
  Int_t labelFirstDau = partD->GetDaughter(0); 
  if(nDau==3){
    for(Int_t iDau=0; iDau<3; iDau++){
      Int_t ind = labelFirstDau+iDau;
      AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(ind));
      if(!part){
	printf("Daughter particle not found in MC array");
	return 99999.;
      } 
      pXdauTrue[iDau]=part->Px();
      pYdauTrue[iDau]=part->Py();
      pZdauTrue[iDau]=part->Pz();
    }
  }else if(nDau==2){
    Int_t theDau=0;
    for(Int_t iDau=0; iDau<2; iDau++){
      Int_t ind = labelFirstDau+iDau;
      AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(arrayMC->At(ind));
      if(!part){
	printf("Daughter particle not found in MC array");
	return 99999.;
      } 
      Int_t pdgCode=TMath::Abs(part->GetPdgCode());
      if(pdgCode==211 || pdgCode==321){
	pXdauTrue[theDau]=part->Px();
	pYdauTrue[theDau]=part->Py();
	pZdauTrue[theDau]=part->Pz();
	++theDau;
      }else{
	Int_t nDauRes=part->GetNDaughters();
	if(nDauRes==2){
	  Int_t labelFirstDauRes = part->GetDaughter(0); 	
	  for(Int_t iDauRes=0; iDauRes<2; iDauRes++){
	    Int_t indDR = labelFirstDauRes+iDauRes;
	    AliAODMCParticle* partDR = dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDR));
	    if(!partDR){
	      printf("Daughter particle not found in MC array");
	      return 99999.;
	    } 
	    
	    Int_t pdgCodeDR=TMath::Abs(partDR->GetPdgCode());
 	    if(pdgCodeDR==211 || pdgCodeDR==321){
	      pXdauTrue[theDau]=partDR->Px();
	      pYdauTrue[theDau]=partDR->Py();
	      pZdauTrue[theDau]=partDR->Pz();
	      ++theDau;
	    }
	  }
	}
      }
    }
    if(theDau!=3){
      printf("Wrong number of decay prongs");
      return 99999.;
    }
  }

  Double_t d0dummy[3]={0.,0.,0.};
  AliAODRecoDecayHF aodDvsMC(vtxTrue,origD,3,charge,pXdauTrue,pYdauTrue,pZdauTrue,d0dummy);
  return aodDvsMC.ImpParXY();

}



//____________________________________________________________________________
Double_t AliVertexingHFUtils::GetCorrectedNtracklets(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult) {
  //
  // Correct the number of accepted tracklets based on a TProfile Hist
  //
  //

  if(TMath::Abs(vtxZ)>10.0){
    //    printf("ERROR: Z vertex out of range for correction of multiplicity\n");
    return uncorrectedNacc;
  }

  if(!estimatorAvg){
    printf("ERROR: Missing TProfile for correction of multiplicity\n");
    return uncorrectedNacc;
  }

  Double_t localAvg = estimatorAvg->GetBinContent(estimatorAvg->FindBin(vtxZ));

  Double_t deltaM = uncorrectedNacc*(refMult/localAvg - 1);

  Double_t correctedNacc = uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->Poisson(TMath::Abs(deltaM));

  return correctedNacc;
}
//______________________________________________________________________
TString AliVertexingHFUtils::GetGenerator(Int_t label, AliAODMCHeader* header){
  // get the name of the generator that produced a given particle

  Int_t nsumpart=0;
  TList *lh=header->GetCocktailHeaders();
  Int_t nh=lh->GetEntries();
  for(Int_t i=0;i<nh;i++){
    AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
    TString genname=gh->GetName();
    Int_t npart=gh->NProduced();
    if(label>=nsumpart && label<(nsumpart+npart)) return genname;
    nsumpart+=npart;
  }
  TString empty="";
  return empty;
}
//_____________________________________________________________________
void AliVertexingHFUtils::GetTrackPrimaryGenerator(AliAODTrack *track,AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen){

  // method to check if a track comes from a given generator

  Int_t lab=TMath::Abs(track->GetLabel());
  nameGen=GetGenerator(lab,header);
  
  //  Int_t countControl=0;
  
  while(nameGen.IsWhitespace()){
    AliAODMCParticle *mcpart= (AliAODMCParticle*)arrayMC->At(lab);
    if(!mcpart){
      printf("AliVertexingHFUtils::IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
      break;
    }
    Int_t mother = mcpart->GetMother();
    if(mother<0){
      printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
      break;
    }
    lab=mother;
    nameGen=GetGenerator(mother,header);
    // countControl++;
    // if(countControl>=10){ // 10 = arbitrary number; protection from infinite loops
    //   printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Protection from infinite loop active\n");
    //   break;
    // }
  }
  
  return;
}
//----------------------------------------------------------------------
Bool_t AliVertexingHFUtils::IsTrackInjected(AliAODTrack *track,AliAODMCHeader *header,TClonesArray *arrayMC){
  // method to check if a track comes from the signal event or from the underlying Hijing event
  TString nameGen;
  GetTrackPrimaryGenerator(track,header,arrayMC,nameGen);
  
  if(nameGen.IsWhitespace() || nameGen.Contains("ijing")) return kFALSE;
  
  return kTRUE;
}
//____________________________________________________________________________
Bool_t AliVertexingHFUtils::IsCandidateInjected(AliAODRecoDecayHF *cand, AliAODMCHeader *header,TClonesArray *arrayMC){
  // method to check if a D meson candidate comes from the signal event or from the underlying Hijing event

  Int_t nprongs=cand->GetNProngs();
  for(Int_t i=0;i<nprongs;i++){
    AliAODTrack *daugh=(AliAODTrack*)cand->GetDaughter(i);
    if(IsTrackInjected(daugh,header,arrayMC)) return kTRUE;
  }
  return kFALSE;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Bool_t searchUpToQuark){
  // checking whether the mother of the particles come from a charm or a bottom quark

  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPart->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;
  while (mother >0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetMother();
    }else{
      printf("AliVertexingHFUtils::CheckOrigin: Failed casting the mother particle!");
      break;
    }
  }
  if(searchUpToQuark && !isQuarkFound) return 0;
  if(isFromB) return 5;
  else return 4;

}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckD0Decay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  // Checks the D0 decay channel. Returns 1 for the D0->Kpi case, 2 for the D0->Kpipipi case, -1 in other cases


  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=421) return -1;

  Int_t nDau=mcPart->GetNDaughters();

  if(nDau==2){
    Int_t daughter0 = mcPart->GetDaughter(0);
    Int_t daughter1 = mcPart->GetDaughter(1);
    AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(daughter0));
    AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(arrayMC->At(daughter1));
    if(!mcPartDaughter0 || !mcPartDaughter1) return -1;
    arrayDauLab[0]=daughter0;
    arrayDauLab[1]=daughter1;
    Int_t pdgCode0=mcPartDaughter0->GetPdgCode();
    Int_t pdgCode1=mcPartDaughter1->GetPdgCode();
    if(!(TMath::Abs(pdgCode0)==321 && TMath::Abs(pdgCode1)==211) &&
       !(TMath::Abs(pdgCode0)==211 && TMath::Abs(pdgCode1)==321)){
      return -1;
    }
    if(TMath::Abs(pdgCode0)==321 && (pdgD*pdgCode0)>0) return -1;
    if(TMath::Abs(pdgCode1)==321 && (pdgD*pdgCode1)>0) return -1;
    if((pdgCode0*pdgCode1)>0) return -1;
    Double_t sumPxDau=mcPartDaughter0->Px()+mcPartDaughter1->Px();
    Double_t sumPyDau=mcPartDaughter0->Py()+mcPartDaughter1->Py();
    Double_t sumPzDau=mcPartDaughter0->Pz()+mcPartDaughter1->Pz();
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    return 1;
  }

  if(nDau==3 || nDau==4){
    Int_t nKaons=0;
    Int_t nPions=0;
    Double_t sumPxDau=0.;
    Double_t sumPyDau=0.;
    Double_t sumPzDau=0.;
    Int_t nFoundKpi=0;
    Int_t labelFirstDau = mcPart->GetDaughter(0);
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	if(pdgD>0 && pdgdau>0) return -1;
	if(pdgD<0 && pdgdau<0) return -1;
	nKaons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>4) return -1;
      }else if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>4) return -1;
      }else if(TMath::Abs(pdgdau)==113 || TMath::Abs(pdgdau)==313){
	Int_t nResDau=dau->GetNDaughters();
	if(nResDau!=2) return -1;
	Int_t indFirstResDau=dau->GetDaughter(0);
	for(Int_t resDau=0; resDau<2; resDau++){
	  Int_t indResDau=indFirstResDau+resDau;
	  if(indResDau<0) return -1;
	  AliAODMCParticle* resdau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indResDau));
	  if(!resdau) return -1;
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    if(pdgD>0 && pdgresdau>0) return -1;
	    if(pdgD<0 && pdgresdau<0) return -1;
	    nKaons++;
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    arrayDauLab[nFoundKpi++]=indResDau;
	    if(nFoundKpi>4) return -1;
	  }
	  if(TMath::Abs(pdgresdau)==211){
	    nPions++;
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    arrayDauLab[nFoundKpi++]=indResDau;
	    if(nFoundKpi>4) return -1;
	  }
	}
      }else{
	return -1;
      }
    }
    if(nPions!=3) return -1;
    if(nKaons!=1) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -1;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -1;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -1;
    return 2;
  }
  return -1;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDplusDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  // Checks the Dplus decay channel. Returns 1 for the non-resonant case, 2 for the resonant case, -1 in other cases

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=411) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughter(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;

  if(nDau==3 || nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	if(pdgD*pdgdau>0) return -1;
	nKaons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==211){
	if(pdgD*pdgdau<0) return -1;
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==313){
	Int_t nResDau=dau->GetNDaughters();
	if(nResDau!=2) return -1;
	Int_t indFirstResDau=dau->GetDaughter(0);
	for(Int_t resDau=0; resDau<2; resDau++){
	  Int_t indResDau=indFirstResDau+resDau;
	  if(indResDau<0) return -1;
	  AliAODMCParticle* resdau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indResDau));
	  if(!resdau) return -1;
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    if(pdgD*pdgresdau>0) return -1;
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nKaons++;
	    arrayDauLab[nFoundKpi++]=indResDau;
	    if(nFoundKpi>3) return -1;
	  }
	  if(TMath::Abs(pdgresdau)==211){
	    if(pdgD*pdgresdau<0) return -1;
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nPions++;
	    arrayDauLab[nFoundKpi++]=indResDau;
	    if(nFoundKpi>3) return -1;
	  }
	}
       }else{
	return -1;
      }
    }
    if(nPions!=2) return -1;
    if(nKaons!=1) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    if(nDau==3) return 1;
    else if(nDau==2) return 2;
  }

  return -1;

}
