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
#include <TParticle.h>
#include "AliMCEvent.h"
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

/// \cond CLASSIMP
ClassImp(AliVertexingHFUtils);
/// \endcond


//______________________________________________________________________
AliVertexingHFUtils::AliVertexingHFUtils():TObject(),
  fK(1),
  fSubRes(1.),
  fMinEtaForTracklets(-1.),
  fMaxEtaForTracklets(1.)
{
  /// Default contructor
}


//______________________________________________________________________
AliVertexingHFUtils::AliVertexingHFUtils(Int_t k):
  TObject(),
  fK(k),
  fSubRes(1.),
  fMinEtaForTracklets(-1.),
  fMaxEtaForTracklets(1.)
{
  /// Standard constructor
}


//______________________________________________________________________
void AliVertexingHFUtils::ComputeSignificance(Double_t signal, Double_t  errsignal, Double_t  background, Double_t  errbackground, Double_t &significance,Double_t &errsignificance){
  /// calculate significance from S, B and errors


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
  /// compute chi from polynomial approximation
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
  /// compute chi variable (=v2*sqrt(N)) from external values

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
  /// computes event plane resolution starting from sub event resolution
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
  /// computes event plane resolution starting from sub event correlation histogram
  if(!hSubEvCorr) return 1.;
  Double_t resSub=GetSubEvResol(hSubEvCorr);
  return GetFullEvResol(resSub,k);
}
//______________________________________________________________________
Double_t AliVertexingHFUtils::GetFullEvResolLowLim(const TH1F* hSubEvCorr, Int_t k){
  /// computes low limit event plane resolution starting from sub event correlation histogram
  if(!hSubEvCorr) return 1.;
  Double_t resSub=GetSubEvResolLowLim(hSubEvCorr);
  printf("%f\n",resSub);
  return GetFullEvResol(resSub,k);  
}
//______________________________________________________________________
Double_t AliVertexingHFUtils::GetFullEvResolHighLim(const TH1F* hSubEvCorr, Int_t k){
  /// computes high limit event plane resolution starting from sub event correlation histogram
  if(!hSubEvCorr) return 1.;
  Double_t resSub=GetSubEvResolHighLim(hSubEvCorr);
  printf("%f\n",resSub);
  return GetFullEvResol(resSub,k);  
}
//______________________________________________________________________
Int_t AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(AliAODEvent* ev, Double_t mineta, Double_t maxeta){
  /// counts tracklets in given eta range
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
  /// counts generated particles in fgiven eta range

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
  /// counts generated primary particles in given eta range

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
  /// counts generated primary particles in given eta range

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
Double_t AliVertexingHFUtils::GetVZEROAEqualizedMultiplicity(AliAODEvent* ev){
  //
  /// Method to get VZERO-A equalized multiplicty as done in AliCentralitySelectionTask
  ///  getting the equalized VZERO factors from the tender or AOD
  // http://git.cern.ch/pubweb/AliRoot.git/blob/HEAD:/ANALYSIS/AliCentralitySelectionTask.cxx#l1345

  Double_t multV0AEq=0;
  for(Int_t iCh = 32; iCh < 64; ++iCh) {
    Double_t mult = ev->GetVZEROEqMultiplicity(iCh);
    multV0AEq += mult;
  }
  return multV0AEq;
}

//______________________________________________________________________
Double_t AliVertexingHFUtils::GetVZEROCEqualizedMultiplicity(AliAODEvent* ev){
  /// Method to get VZERO-C equalized multiplicty as done in AliCentralitySelectionTask
  ///  getting the equalized VZERO factors from the tender or AOD
  /// http://git.cern.ch/pubweb/AliRoot.git/blob/HEAD:/ANALYSIS/AliCentralitySelectionTask.cxx#l1345

  Double_t multV0CEq=0;
  for(Int_t iCh = 0; iCh < 32; ++iCh) {
    Double_t mult = ev->GetVZEROEqMultiplicity(iCh);
    multV0CEq += mult;
  }
  return multV0CEq;
}

//______________________________________________________________________
void AliVertexingHFUtils::AveragePt(Float_t& averagePt, Float_t& errorPt,Float_t ptmin,Float_t ptmax, TH2F* hMassD, Float_t massFromFit, Float_t sigmaFromFit, TF1* funcB2, Float_t sigmaRangeForSig,Float_t sigmaRangeForBkg, Float_t minMass, Float_t maxMass, Int_t rebin){

  /// Compute <pt> from 2D histogram M vs pt

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
Bool_t AliVertexingHFUtils::CheckT0TriggerFired(AliAODEvent* aodEv){
  /// check if T0VTX trigger was fired, based on a workaround suggested by Alla
  const Double32_t *mean = aodEv->GetT0TOF();
  if(mean && mean[0]<9999.) return kTRUE;
  else return kFALSE;
}
//____________________________________________________________________________
Double_t AliVertexingHFUtils::GetTrueImpactParameterDzero(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partD) {
  /// true impact parameter calculation for Dzero

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
  /// true impact parameter calculation for Dplus

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
  /// get the name of the generator that produced a given particle

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

  /// method to check if a track comes from a given generator

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
  /// method to check if a track comes from the signal event or from the underlying Hijing event
  TString nameGen;
  GetTrackPrimaryGenerator(track,header,arrayMC,nameGen);
  
  if(nameGen.IsWhitespace() || nameGen.Contains("ijing")) return kFALSE;
  
  return kTRUE;
}
//____________________________________________________________________________
Bool_t AliVertexingHFUtils::IsCandidateInjected(AliAODRecoDecayHF *cand, AliAODMCHeader *header,TClonesArray *arrayMC){
  /// method to check if a D meson candidate comes from the signal event or from the underlying Hijing event

  Int_t nprongs=cand->GetNProngs();
  for(Int_t i=0;i<nprongs;i++){
    AliAODTrack *daugh=(AliAODTrack*)cand->GetDaughter(i);
    if(IsTrackInjected(daugh,header,arrayMC)) return kTRUE;
  }
  return kFALSE;
}
//____________________________________________________________________________
Bool_t AliVertexingHFUtils::HasCascadeCandidateAnyDaughInjected(AliAODRecoCascadeHF *cand, AliAODMCHeader *header,TClonesArray *arrayMC){
  /// method to check if a cascade candidate comes from the signal event or from the underlying Hijing event

  AliAODTrack* bach = cand->GetBachelor();
  if(IsTrackInjected(bach, header, arrayMC)) {
    AliDebug(2, "Bachelor is injected, the whole candidate is then injected");
    return kTRUE;
  }
  AliAODv0* v0 = cand->Getv0();
  Int_t nprongs = v0->GetNProngs();
  for(Int_t i = 0; i < nprongs; i++){
    AliAODTrack *daugh = (AliAODTrack*)v0->GetDaughter(i);
    if(IsTrackInjected(daugh,header,arrayMC)) {
      AliDebug(2, Form("V0 daughter number %d is injected, the whole candidate is then injected", i));
      return kTRUE;
    }
  }
  return kFALSE;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckOrigin(AliMCEvent* mcEvent, TParticle *mcPart, Bool_t searchUpToQuark){
  /// checking whether the mother of the particles come from a charm or a bottom quark

  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPart->GetFirstMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;
  while (mother >0 ){
    istep++;
    TParticle* mcGranma = mcEvent->Particle(mother);
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetFirstMother();
    }else{
      printf("CheckOrigin: Failed casting the mother particle!");
      break;
    }
  }
  if(searchUpToQuark && !isQuarkFound) return 0;
  if(isFromB) return 5;
  else return 4;

}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Bool_t searchUpToQuark){
  /// checking whether the mother of the particles come from a charm or a bottom quark

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
Double_t AliVertexingHFUtils::GetBeautyMotherPt(TClonesArray* arrayMC, AliAODMCParticle *mcPart){
  /// get the pt of the beauty hadron (feed-down case), returns negative value for prompt

  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPart->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  while (mother >0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	return mcGranma->Pt();
      }
      if(abspdgGranma==4) return -999.;
      if(abspdgGranma==5) return -1.;
      mother = mcGranma->GetMother();
    }else{
      printf("AliVertexingHFUtils::GetBeautyMotherPt: Failed casting the mother particle!");
      break;
    }
  }
  return -999.;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckD0Decay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the D0 decay channel. Returns 1 for the D0->Kpi case, 2 for the D0->Kpipipi case, -1 in other cases

  if(label<0) return -1;
  TParticle* mcPart = mcEvent->Particle(label);
  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=421) return -1;
  
  Int_t nDau=mcPart->GetNDaughters();
  
  if(nDau==2){
    Int_t daughter0 = mcPart->GetDaughter(0);
    Int_t daughter1 = mcPart->GetDaughter(1);
    TParticle* mcPartDaughter0 = mcEvent->Particle(daughter0);
    TParticle* mcPartDaughter1 = mcEvent->Particle(daughter1);
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
      TParticle* dau=mcEvent->Particle(indDau);
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
	  TParticle* resdau=mcEvent->Particle(indResDau);
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
Int_t AliVertexingHFUtils::CheckD0Decay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
 /// Checks the D0 decay channel. Returns 1 for the D0->Kpi case, 2 for the D0->Kpipipi case, -1 in other cases

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
Int_t AliVertexingHFUtils::CheckDplusDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Dplus decay channel. Returns 1 for the non-resonant case, 2 for the resonant case, -1 in other cases

  if(label<0) return -1;
  TParticle* mcPart = mcEvent->Particle(label);
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
      TParticle* dau=mcEvent->Particle(indDau);
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
	  TParticle* resdau=mcEvent->Particle(indResDau);
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

//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDplusDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  /// Checks the Dplus decay channel. Returns 1 for the non-resonant case, 2 for the resonant case, -1 in other cases

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
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDplusKKpiDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Ds decay channel. Returns 1 for Ds->phipi->KKpi, 2 for Ds->K0*K->KKpi, 3 for the non-resonant case

  if(label<0) return -1;
  TParticle* mcPart = mcEvent->Particle(label);
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
  Bool_t isPhi=kFALSE;
  Bool_t isk0st=kFALSE;
  
  if(nDau==3 || nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      TParticle* dau=mcEvent->Particle(indDau);
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	nKaons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==313 || TMath::Abs(pdgdau)==333){
	if(TMath::Abs(pdgdau)==313) isk0st=kTRUE;
	if(TMath::Abs(pdgdau)==333) isPhi=kTRUE;
	Int_t nResDau=dau->GetNDaughters();
	if(nResDau!=2) return -1;
	Int_t indFirstResDau=dau->GetDaughter(0);
	for(Int_t resDau=0; resDau<2; resDau++){
	  Int_t indResDau=indFirstResDau+resDau;
	  if(indResDau<0) return -1;
	  TParticle* resdau=mcEvent->Particle(indResDau);
	  if(!resdau) return -1;
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nKaons++;
	    arrayDauLab[nFoundKpi++]=indResDau;
	    if(nFoundKpi>3) return -1;
	  }
	  if(TMath::Abs(pdgresdau)==211){
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
    if(nPions!=1) return -1;
    if(nKaons!=2) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    if(nDau==3) return 3;
    else if(nDau==2){
      if(isk0st) return 2;
      if(isPhi) return 1;
    }
  }
  
  return -1;
  
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDplusKKpiDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  /// Checks the D+ decay channel. Returns 1 for D+->phipi->KKpi, 2 for D+->K0*K->KKpi, 3 for the non-resonant case
    
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
  Bool_t isPhi=kFALSE;
  Bool_t isk0st=kFALSE;
    
  if(nDau==3 || nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	nKaons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==313 || TMath::Abs(pdgdau)==333){
	if(TMath::Abs(pdgdau)==313) isk0st=kTRUE;
	if(TMath::Abs(pdgdau)==333) isPhi=kTRUE;
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
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nKaons++;
	    arrayDauLab[nFoundKpi++]=indResDau;
	    if(nFoundKpi>3) return -1;
	  }
	  if(TMath::Abs(pdgresdau)==211){
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
    if(nPions!=1) return -1;
    if(nKaons!=2) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    if(nDau==3) return 3;
    else if(nDau==2){
      if(isk0st) return 2;
      if(isPhi) return 1;
    }
  }
    
  return -1;
    
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDplusK0spiDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Dplus->V0+pion decay channel. Returns 1 if success, -1 otherwise

  if(label<0) return -1;
  TParticle* mcPart = mcEvent->Particle(label);
  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=411) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughter(0);
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundpi=0;

  Int_t codeV0=-1;
  if(nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      TParticle* dau=mcEvent->Particle(indDau);
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpi++]=indDau;
	if(nFoundpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==311){
	codeV0=TMath::Abs(pdgdau);
	TParticle* v0=dau;
	if(codeV0==311){
	  Int_t nK0Dau=dau->GetNDaughters();
	  if(nK0Dau!=1) return -1;
	  Int_t indK0s=dau->GetDaughter(0);
	  if(indK0s<0) return -1;
	  v0=mcEvent->Particle(indK0s);
	  if(!v0) return -1;
	  Int_t pdgK0sl=v0->GetPdgCode();
	  codeV0=TMath::Abs(pdgK0sl);
	}
	Int_t nV0Dau=v0->GetNDaughters();
	if(nV0Dau!=2) return -1;
	Int_t indFirstV0Dau=v0->GetDaughter(0);
	for(Int_t v0Dau=0; v0Dau<2; v0Dau++){
	  Int_t indV0Dau=indFirstV0Dau+v0Dau;
	  if(indV0Dau<0) return -1;
	  TParticle* v0dau=mcEvent->Particle(indV0Dau);
	  if(!v0dau) return -1;
	  Int_t pdgv0dau=v0dau->GetPdgCode();
	  if(TMath::Abs(pdgv0dau)==211){
	    sumPxDau+=v0dau->Px();
	    sumPyDau+=v0dau->Py();
	    sumPzDau+=v0dau->Pz();
	    nPions++;
	    arrayDauLab[nFoundpi++]=indV0Dau;
	    if(nFoundpi>3) return -1;
	  }
	}
      }else{
	return -1;
      }
    }
    if(nPions!=3) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    if(codeV0==310) return 1;
  }
  return -1;
  
}

 
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDsDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Ds decay channel. Returns 1 for Ds->phipi->KKpi, 2 for Ds->K0*K->KKpi, 3 for the non-resonant case

  if(label<0) return -1;
  TParticle* mcPart = mcEvent->Particle(label);
  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=431) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughter(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;
  Bool_t isPhi=kFALSE;
  Bool_t isk0st=kFALSE;
  
  if(nDau==3 || nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      TParticle* dau=mcEvent->Particle(indDau);
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	nKaons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==313 || TMath::Abs(pdgdau)==333){
	if(TMath::Abs(pdgdau)==313) isk0st=kTRUE;
	if(TMath::Abs(pdgdau)==333) isPhi=kTRUE;
	Int_t nResDau=dau->GetNDaughters();
	if(nResDau!=2) return -1;
	Int_t indFirstResDau=dau->GetDaughter(0);
	for(Int_t resDau=0; resDau<2; resDau++){
	  Int_t indResDau=indFirstResDau+resDau;
	  if(indResDau<0) return -1;
	  TParticle* resdau=mcEvent->Particle(indResDau);
	  if(!resdau) return -1;
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nKaons++;
	    arrayDauLab[nFoundKpi++]=indResDau;
	    if(nFoundKpi>3) return -1;
	  }
	  if(TMath::Abs(pdgresdau)==211){
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
    if(nPions!=1) return -1;
    if(nKaons!=2) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    if(nDau==3) return 3;
    else if(nDau==2){
      if(isk0st) return 2;
      if(isPhi) return 1;
    }
  }
  
  return -1;
  
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDsK0sKDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Ds->K0s+S decay channel. Returns 1 in case of success, otherwise -1

  if(label<0) return -1;
  TParticle* mcPart = mcEvent->Particle(label);
  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=431) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughter(0);
  Int_t nPions=0;
  Int_t nKaons=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;

  Int_t codeV0=-1;
  if(nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      TParticle* dau=mcEvent->Particle(indDau);
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==321){
	nKaons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==311){
	codeV0=TMath::Abs(pdgdau);
	TParticle* v0=dau;
	if(codeV0==311){
	  Int_t nK0Dau=dau->GetNDaughters();
	  if(nK0Dau!=1) return -1;
	  Int_t indK0s=dau->GetDaughter(0);
	  if(indK0s<0) return -1;
	  v0=mcEvent->Particle(indK0s);
	  if(!v0) return -1;
	  Int_t pdgK0sl=v0->GetPdgCode();
	  codeV0=TMath::Abs(pdgK0sl);
	}
	Int_t nV0Dau=v0->GetNDaughters();
	if(nV0Dau!=2) return -1;
	Int_t indFirstV0Dau=v0->GetDaughter(0);
	for(Int_t v0Dau=0; v0Dau<2; v0Dau++){
	  Int_t indV0Dau=indFirstV0Dau+v0Dau;
	  if(indV0Dau<0) return -1;
	  TParticle* v0dau=mcEvent->Particle(indV0Dau);
	  if(!v0dau) return -1;
	  Int_t pdgv0dau=v0dau->GetPdgCode();
	  if(TMath::Abs(pdgv0dau)==211){
	    sumPxDau+=v0dau->Px();
	    sumPyDau+=v0dau->Py();
	    sumPzDau+=v0dau->Pz();
	    nPions++;
	    arrayDauLab[nFoundKpi++]=indV0Dau;
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
    if(codeV0==310) return 1;
  }
  return -1;
  
}

 
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDsDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  /// Checks the Ds decay channel. Returns 1 for Ds->phipi->KKpi, 2 for Ds->K0*K->KKpi, 3 for the non-resonant case

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=431) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughter(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;
  Bool_t isPhi=kFALSE;
  Bool_t isk0st=kFALSE;
  
  if(nDau==3 || nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	nKaons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundKpi++]=indDau;
	if(nFoundKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==313 || TMath::Abs(pdgdau)==333){
	if(TMath::Abs(pdgdau)==313) isk0st=kTRUE;
	if(TMath::Abs(pdgdau)==333) isPhi=kTRUE;
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
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nKaons++;
	    arrayDauLab[nFoundKpi++]=indResDau;
	    if(nFoundKpi>3) return -1;
	  }
	  if(TMath::Abs(pdgresdau)==211){
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
    if(nPions!=1) return -1;
    if(nKaons!=2) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    if(nDau==3) return 3;
    else if(nDau==2){
      if(isk0st) return 2;
      if(isPhi) return 1;
    }
  }
  
  return -1;
  
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDstarDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Dstar decay channel. Returns 1 for D*->D0pi->Kpipi, -1 in other cases

  if(label<0) return -1;
  TParticle* mcPart = mcEvent->Particle(label);
  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=413) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughter(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;

  for(Int_t iDau=0; iDau<nDau; iDau++){
    Int_t indDau = labelFirstDau+iDau;
    if(indDau<0) return -1;
    TParticle* dau=mcEvent->Particle(indDau);
    if(!dau) return -1;
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==421){
      Int_t nResDau=dau->GetNDaughters();
      if(nResDau!=2) return -1;
      Int_t indFirstResDau=dau->GetDaughter(0);
      for(Int_t resDau=0; resDau<2; resDau++){
	Int_t indResDau=indFirstResDau+resDau;
	if(indResDau<0) return -1;
	TParticle* resdau=mcEvent->Particle(indResDau);
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
    }else if(TMath::Abs(pdgdau)==211){
      if(pdgD*pdgdau<0) return -1;
      nPions++;
      sumPxDau+=dau->Px();
      sumPyDau+=dau->Py();
      sumPzDau+=dau->Pz();
      arrayDauLab[nFoundKpi++]=indDau;
      if(nFoundKpi>3) return -1;
    }
  }

  if(nPions!=2) return -1;
  if(nKaons!=1) return -1;
  if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
  if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
  if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
  return 1;
  
}

//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDstarDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  /// Checks the Dstar decay channel. Returns 1 for D*->D0pi->Kpipi, -1 in other cases

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=413) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughter(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;

  for(Int_t iDau=0; iDau<nDau; iDau++){
    Int_t indDau = labelFirstDau+iDau;
    if(indDau<0) return -1;
    AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
    if(!dau) return -1;
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==421){
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
    }else if(TMath::Abs(pdgdau)==211){
      if(pdgD*pdgdau<0) return -1;
      nPions++;
      sumPxDau+=dau->Px();
      sumPyDau+=dau->Py();
      sumPzDau+=dau->Pz();
      arrayDauLab[nFoundKpi++]=indDau;
      if(nFoundKpi>3) return -1;
    }
  }

  if(nPions!=2) return -1;
  if(nKaons!=1) return -1;
  if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
  if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
  if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
  return 1;
  
}

//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckLcpKpiDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Lc->pKpi decay channel. Returns 1 for non-resonant decays and 2, 3 or 4 for resonant ones, -1 in other cases

  if(label<0) return -1;
  TParticle* mcPart = mcEvent->Particle(label);
  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=4122) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughter(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Int_t nProtons=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundpKpi=0;

  Int_t codeRes=-1;
  if(nDau==3 || nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      TParticle* dau=mcEvent->Particle(indDau);
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	nKaons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==2212){
	nProtons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==313 || TMath::Abs(pdgdau)==3124 || 
	       TMath::Abs(pdgdau)==2224){
	codeRes=TMath::Abs(pdgdau);
	Int_t nResDau=dau->GetNDaughters();
	if(nResDau!=2) return -1;
	Int_t indFirstResDau=dau->GetDaughter(0);
	for(Int_t resDau=0; resDau<2; resDau++){
	  Int_t indResDau=indFirstResDau+resDau;
	  if(indResDau<0) return -1;
	  TParticle* resdau=mcEvent->Particle(indResDau);
	  if(!resdau) return -1;
	  Int_t pdgresdau=resdau->GetPdgCode();
	  if(TMath::Abs(pdgresdau)==321){
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nKaons++;
	    arrayDauLab[nFoundpKpi++]=indResDau;
	    if(nFoundpKpi>3) return -1;
	  }else if(TMath::Abs(pdgresdau)==211){
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nPions++;
	    arrayDauLab[nFoundpKpi++]=indResDau;
	    if(nFoundpKpi>3) return -1;
	  }else if(TMath::Abs(pdgresdau)==2212){
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nProtons++;
	    arrayDauLab[nFoundpKpi++]=indResDau;
	    if(nFoundpKpi>3) return -1;
	  }
	}
      }else{
	return -1;
      }
    }
    if(nPions!=1) return -1;
    if(nKaons!=1) return -1;
    if(nProtons!=1) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    if(nDau==3) return 1;
    else if(nDau==2){
      if(codeRes==313) return 2;
      else if(codeRes==2224) return 3;
      else if(codeRes==3124) return 4;
    }  
  }
  return -1;
  
}

//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckLcpKpiDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  /// Checks the Lc->pKpi decay channel. Returns 1 for non-resonant decays and 2, 3 or 4 for resonant ones, -1 in other cases

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=4122) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughter(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Int_t nProtons=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundpKpi=0;

  Int_t codeRes=-1;
  if(nDau==3 || nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==321){
	nKaons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==2212){
	nProtons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundpKpi++]=indDau;
	if(nFoundpKpi>3) return -1;
      }else if(TMath::Abs(pdgdau)==313 || TMath::Abs(pdgdau)==3124 || 
	       TMath::Abs(pdgdau)==2224){
	codeRes=TMath::Abs(pdgdau);
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
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nKaons++;
	    arrayDauLab[nFoundpKpi++]=indResDau;
	    if(nFoundpKpi>3) return -1;
	  }else if(TMath::Abs(pdgresdau)==211){
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nPions++;
	    arrayDauLab[nFoundpKpi++]=indResDau;
	    if(nFoundpKpi>3) return -1;
	  }else if(TMath::Abs(pdgresdau)==2212){
	    sumPxDau+=resdau->Px();
	    sumPyDau+=resdau->Py();
	    sumPzDau+=resdau->Pz();
	    nProtons++;
	    arrayDauLab[nFoundpKpi++]=indResDau;
	    if(nFoundpKpi>3) return -1;
	  }
	}
      }else{
	return -1;
      }
    }
    if(nPions!=1) return -1;
    if(nKaons!=1) return -1;
    if(nProtons!=1) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    if(nDau==3) return 1;
    else if(nDau==2){
      if(codeRes==313) return 2;
      else if(codeRes==2224) return 3;
      else if(codeRes==3124) return 4;
    }  
  }
  return -1;
  
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckLcV0bachelorDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Lc->V0+bachelor decay channel. Returns 1 for pK0s, 2 for piLambda, 3 for pK0l -1 in other cases

  if(label<0) return -1;
  TParticle* mcPart = mcEvent->Particle(label);
  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=4122) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughter(0);
  Int_t nPions=0;
  Int_t nProtons=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundppi=0;

  Int_t codeV0=-1;
  if(nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      TParticle* dau=mcEvent->Particle(indDau);
      if(!dau) return -1;
      Int_t pdgdau=dau->GetPdgCode();
      if(TMath::Abs(pdgdau)==211){
	nPions++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundppi++]=indDau;
	if(nFoundppi>3) return -1;
      }else if(TMath::Abs(pdgdau)==2212){
	nProtons++;
	sumPxDau+=dau->Px();
	sumPyDau+=dau->Py();
	sumPzDau+=dau->Pz();
	arrayDauLab[nFoundppi++]=indDau;
	if(nFoundppi>3) return -1;
      }else if(TMath::Abs(pdgdau)==311 ||  TMath::Abs(pdgdau)==3122){
	codeV0=TMath::Abs(pdgdau);
	TParticle* v0=dau;
	if(codeV0==311){
	  Int_t nK0Dau=dau->GetNDaughters();
	  if(nK0Dau!=1) return -1;
	  Int_t indK0s=dau->GetDaughter(0);
	  if(indK0s<0) return -1;
	  v0=mcEvent->Particle(indK0s);
	  if(!v0) return -1;
	  Int_t pdgK0sl=v0->GetPdgCode();
	  codeV0=TMath::Abs(pdgK0sl);
	}
	Int_t nV0Dau=v0->GetNDaughters();
	if(nV0Dau!=2) return -1;
	Int_t indFirstV0Dau=v0->GetDaughter(0);
	for(Int_t v0Dau=0; v0Dau<2; v0Dau++){
	  Int_t indV0Dau=indFirstV0Dau+v0Dau;
	  if(indV0Dau<0) return -1;
	  TParticle* v0dau=mcEvent->Particle(indV0Dau);
	  if(!v0dau) return -1;
	  Int_t pdgv0dau=v0dau->GetPdgCode();
	  if(TMath::Abs(pdgv0dau)==211){
	    sumPxDau+=v0dau->Px();
	    sumPyDau+=v0dau->Py();
	    sumPzDau+=v0dau->Pz();
	    nPions++;
	    arrayDauLab[nFoundppi++]=indV0Dau;
	    if(nFoundppi>3) return -1;
	  }else if(TMath::Abs(pdgv0dau)==2212){
	    sumPxDau+=v0dau->Px();
	    sumPyDau+=v0dau->Py();
	    sumPzDau+=v0dau->Pz();
	    nProtons++;
	    arrayDauLab[nFoundppi++]=indV0Dau;
	    if(nFoundppi>3) return -1;
	  }
	}
      }else{
	return -1;
      }
    }
    if(nPions!=2) return -1;
    if(nProtons!=1) return -1;
    if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
    if(codeV0==310) return 1;
    else if(codeV0==3122) return 2;
  }
  return -1;
 
}

 
//__________________________________xic______________________________________
Int_t AliVertexingHFUtils::CheckXicXipipiDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Xic decay channel. Returns 1 for Xic->Xipipi, -1 in other cases

  if(label<0) return -1;
  TParticle* mcPart = mcEvent->Particle(label);
  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=4232) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=3) return -1;

  Int_t labelFirstDau = mcPart->GetDaughter(0);
  Int_t nXi=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundXi=0;

  for(Int_t iDau=0; iDau<nDau; iDau++){
    Int_t indDau = labelFirstDau+iDau;
    if(indDau<0) return -1;
    TParticle* dau=mcEvent->Particle(indDau);
    if(!dau) return -1;
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==3312){
      if(pdgD*pdgdau<0) return -1;
      sumPxDau+=dau->Px();
      sumPyDau+=dau->Py();
      sumPzDau+=dau->Pz();
      nXi++;
      arrayDauLab[nFoundXi++]=indDau;
      
    }
    if(TMath::Abs(pdgdau)==211){
      if(pdgD*pdgdau<0) return -1;
      nPions++;
      sumPxDau+=dau->Px();
      sumPyDau+=dau->Py();
      sumPzDau+=dau->Pz();
      arrayDauLab[nFoundXi++]=indDau;
      if(nFoundXi>3) return -1;
    }
  }
  
  if(nPions!=2) return -1;
  if(nXi!=1) return -1;
  if(TMath::Abs(mcPart->Px()-sumPxDau)>0.001) return -2;
  if(TMath::Abs(mcPart->Py()-sumPyDau)>0.001) return -2;
  if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.001) return -2;
  return 1;
  
}
//________________________________________________________________________
Double_t AliVertexingHFUtils::GetSphericity(AliAODEvent* aod, Double_t etaMin, Double_t etaMax, 
					    Double_t ptMin, Double_t ptMax,
					    Int_t filtbit1, Int_t filtbit2, 
					    Int_t minMult){
  /// compute sphericity

  Int_t nTracks=aod->GetNumberOfTracks();
  Int_t nSelTracks=0;

  Double_t sumpt=0.;
  Double_t s00=0.;
  Double_t s01=0.;
  Double_t s11=0.;
  if(ptMin<0.) ptMin=0.;
  
  for(Int_t it=0; it<nTracks; it++) {
    AliAODTrack *tr=dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
    if(!tr) continue;
    Float_t eta = tr->Eta();
    Float_t pt  = tr->Pt();
    Float_t phi  = tr->Phi();
    if(eta<etaMin || eta>etaMax) continue;
    if(pt<ptMin || pt>ptMax) continue;
    Bool_t fb1 = tr->TestFilterBit(filtbit1);
    Bool_t fb2 = tr->TestFilterBit(filtbit2);
    Bool_t tpcRefit=tr->GetStatus() & AliAODTrack::kTPCrefit;
    if(filtbit1==1 && !tpcRefit) fb1=kFALSE;
    if(filtbit2==1 && !tpcRefit) fb2=kFALSE;
    if( !(fb1 || fb2) ) continue;    
    Double_t px=pt*TMath::Cos(phi);
    Double_t py=pt*TMath::Sin(phi);
    s00 += (px * px)/pt;
    s01 += (py * px)/pt;
    s11 += (py * py)/pt;
    nSelTracks++;
    sumpt+=pt;
  }

  if(nSelTracks<minMult) return -0.5;

  if(sumpt>0.){
    s00/=sumpt;
    s01/=sumpt;
    s11/=sumpt;
  }else return -0.5;

  Double_t sphericity = -10;
  Double_t lambda1=((s00+s11)+TMath::Sqrt((s00+s11)*(s00+s11)-4*(s00*s11-s01*s01)))/2.;
  Double_t lambda2=((s00+s11)-TMath::Sqrt((s00+s11)*(s00+s11)-4*(s00*s11-s01*s01)))/2.;
  if(TMath::Abs(lambda2)<0.00001 && TMath::Abs(lambda1)<0.00001) sphericity=0;
  if(TMath::Abs(lambda1+lambda2)>0.000001) sphericity=2*TMath::Min(lambda1,lambda2)/(lambda1+lambda2);
  return sphericity;

}

//________________________________________________________________________
Double_t AliVertexingHFUtils::GetSpherocity(AliAODEvent* aod, 
					    Double_t etaMin, Double_t etaMax, 
					    Double_t ptMin, Double_t ptMax,
					    Int_t filtbit1, Int_t filtbit2, 
					    Int_t minMult, Double_t phiStepSizeDeg,
					    Int_t nTrksToSkip, Int_t* idToSkip
					    ){
  /// compute spherocity

  Int_t nTracks=aod->GetNumberOfTracks();
  Int_t nSelTracks=0;

  Double_t* ptArr=new Double_t[nTracks];
  Double_t* phiArr=new Double_t[nTracks];
  Double_t sumpt=0.;

  for(Int_t it=0; it<nTracks; it++) {
    AliAODTrack *tr=dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
    if(!tr) continue;
    Float_t eta = tr->Eta();
    Float_t pt  = tr->Pt();
    Float_t phi  = tr->Phi();
    if(eta<etaMin || eta>etaMax) continue;
    if(pt<ptMin || pt>ptMax) continue;    
    if(nTrksToSkip>0 && idToSkip){ 
      Int_t trid = (Int_t)tr->GetID();
      Bool_t keep=kTRUE;
      for(Int_t jt=0; jt<nTrksToSkip; jt++){ 
	if(trid==idToSkip[jt]) keep=kFALSE;
      }
      if(!keep) continue;
    }
    Bool_t fb1 = tr->TestFilterBit(filtbit1);
    Bool_t fb2 = tr->TestFilterBit(filtbit2);
    Bool_t tpcRefit=tr->GetStatus() & AliAODTrack::kTPCrefit;
    if(filtbit1==1 && !tpcRefit) fb1=kFALSE;
    if(filtbit2==1 && !tpcRefit) fb2=kFALSE;
    if( !(fb1 || fb2) ) continue;    
    ptArr[nSelTracks]=pt;
    phiArr[nSelTracks]=phi;
    nSelTracks++;
    sumpt+=pt;
  }

  if(nSelTracks<minMult) return -0.5;

  //Getting thrust
  Double_t spherocity=2.;
  for(Int_t i=0; i<360/phiStepSizeDeg; ++i){
    Double_t phistep=TMath::Pi()*(Double_t)i*phiStepSizeDeg/180.;
    Double_t nx=TMath::Cos(phistep);
    Double_t ny=TMath::Sin(phistep);
    Double_t numer=0.;
    for(Int_t j=0; j<nSelTracks; ++j){
      Double_t pxA=ptArr[j]*TMath::Cos(phiArr[j]);  // x component of an unitary vector n
      Double_t pyA=ptArr[j]*TMath::Sin(phiArr[j]);  // y component of an unitary vector n
      numer+=TMath::Abs(ny*pxA - nx*pyA);  
    }
    Double_t pFull=numer*numer/(sumpt*sumpt);
    if(pFull<spherocity) spherocity=pFull; // minimization;
  }

  delete [] ptArr;
  delete [] phiArr;

  spherocity*=(TMath::Pi()*TMath::Pi()/4.);
  return spherocity;

}
//________________________________________________________________________
Double_t AliVertexingHFUtils::GetGeneratedSpherocity(TClonesArray *arrayMC, 
						     Double_t etaMin, Double_t etaMax, 
						     Double_t ptMin, Double_t ptMax,
						     Int_t minMult, Double_t phiStepSizeDeg){

  /// compute generated spherocity

  Int_t nParticles=arrayMC->GetEntriesFast();
  Int_t nSelParticles=0;

  Double_t* ptArr=new Double_t[nParticles];
  Double_t* phiArr=new Double_t[nParticles];
  Double_t sumpt=0.;

  for(Int_t ip=0; ip<nParticles; ip++) {
    AliAODMCParticle *part=(AliAODMCParticle*)arrayMC->UncheckedAt(ip);
    if(!part) continue;
    Float_t eta = part->Eta();
    Float_t pt  = part->Pt();
    Float_t phi  = part->Phi();
    Int_t charge = part->Charge();
    Bool_t isPhysPrim = part->IsPhysicalPrimary();
    if(!isPhysPrim) continue;
    if(charge==0) continue;
    if(eta<etaMin || eta>etaMax) continue;
    if(pt<ptMin || pt>ptMax) continue;    

    ptArr[nSelParticles]=pt;
    phiArr[nSelParticles]=phi;
    nSelParticles++;
    sumpt+=pt;
  }

  if(nSelParticles<minMult) return -0.5;

  //Getting thrust
  Double_t spherocity=2.;
  for(Int_t i=0; i<360/phiStepSizeDeg; ++i){
    Double_t phistep=TMath::Pi()*(Double_t)i*phiStepSizeDeg/180.;
    Double_t nx=TMath::Cos(phistep);
    Double_t ny=TMath::Sin(phistep);
    Double_t numer=0.;
    for(Int_t j=0; j<nSelParticles; ++j){
      Double_t pxA=ptArr[j]*TMath::Cos(phiArr[j]);  // x component of an unitary vector n
      Double_t pyA=ptArr[j]*TMath::Sin(phiArr[j]);  // y component of an unitary vector n
      numer+=TMath::Abs(ny*pxA - nx*pyA);  
    }
    Double_t pFull=numer*numer/(sumpt*sumpt);
    if(pFull<spherocity) spherocity=pFull; // minimization;
  }

  delete [] ptArr;
  delete [] phiArr;

  spherocity*=(TMath::Pi()*TMath::Pi()/4.);
  return spherocity;

}
