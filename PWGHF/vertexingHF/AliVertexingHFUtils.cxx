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
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TParticle.h>
#include <TLorentzVector.h>
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF.h"
#include "AliVertexingHFUtils.h"

#ifndef HomogeneousField
#define HomogeneousField 
#endif

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class with functions useful for different D2H analyses        //
// - event plane resolution                                      //
// - <pt> calculation with side band subtraction                 //
// - tracklet multiplicity calculation                           //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
// - KF particle                                                 //
// Origin: Jianhui Zhu, <zjh@mail.ccnu.edu.cn>                   //
//                                                               //
///////////////////////////////////////////////////////////////////

/// \cond CLASSIMP
ClassImp(AliVertexingHFUtils);
/// \endcond

using namespace std;

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
void AliVertexingHFUtils::AveragePt(Float_t& averagePt, Float_t& errorPt,Float_t ptmin,Float_t ptmax, TH2F* hMassD, Float_t massFromFit, Float_t sigmaFromFit, 
                                    TF1* funcB2, Float_t sigmaRangeForSig,Float_t sigmaRangeForBkg, Float_t minMass, Float_t maxMass, Int_t rebin){

  /// Compute <pt> from 2D histogram M vs pt

  //Make 2D histos in the desired pt range
  Int_t start=hMassD->FindBin(ptmin);
  Int_t end=hMassD->FindBin(ptmax)-1;
  const Int_t nx=end-start;
  TH2F *hMassDpt=new TH2F("hptmass","hptmass",nx,ptmin,ptmax,hMassD->GetNbinsY(),hMassD->GetYaxis()->GetBinLowEdge(1),
                          hMassD->GetYaxis()->GetBinLowEdge(hMassD->GetNbinsY())+hMassD->GetYaxis()->GetBinWidth(hMassD->GetNbinsY()));
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
  Int_t labelFirstDau = partD->GetDaughterLabel(0);
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
  Int_t labelFirstDau = partD->GetDaughterLabel(0);
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
	  Int_t labelFirstDauRes = part->GetDaughterLabel(0);
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

  Double_t correctedNacc = std::max(0., uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->Poisson(TMath::Abs(deltaM)));

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
  GetTrackPrimaryGenerator(track->GetLabel(),header,arrayMC,nameGen);
}
//_____________________________________________________________________
void AliVertexingHFUtils::GetTrackPrimaryGenerator(Int_t label, AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen){

  /// method to check if a track comes from a given generator

  Int_t lab=TMath::Abs(label);
  nameGen=GetGenerator(lab,header);

  //  Int_t countControl=0;

  while(nameGen.IsWhitespace()){
    AliAODMCParticle *mcpart= (AliAODMCParticle*)arrayMC->At(lab);
    if(!mcpart) break;
    Int_t mother = mcpart->GetMother();
    if(mother<0) break;
    lab=mother;
    nameGen=GetGenerator(mother,header);
    // countControl++;
    // if(countControl>=10){ // 10 = arbitrary number; protection from infinite loops
    //   printf("AliVertexingHFUtils::GetTrackPrimaryGenerator - BREAK: Protection from infinite loop active\n");
    //   break;
    // }
  }

  return;
}
//----------------------------------------------------------------------
Bool_t AliVertexingHFUtils::IsTrackInjected(AliAODTrack *track,AliAODMCHeader *header,TClonesArray *arrayMC){
  /// method to check if a track comes from the signal event or from the underlying Hijing event

  return IsTrackInjected(track->GetLabel(),header,arrayMC);
}
//----------------------------------------------------------------------
Bool_t AliVertexingHFUtils::IsTrackInjected(Int_t label, AliAODMCHeader *header,TClonesArray *arrayMC){
  /// method to check if a track comes from the signal event or from the underlying Hijing event
  TString nameGen;
  Int_t lab=TMath::Abs(label);

  GetTrackPrimaryGenerator(lab,header,arrayMC,nameGen);

  if(nameGen.IsWhitespace() || nameGen.Contains("ijing")) return kFALSE;

  return kTRUE;
}
//____________________________________________________________________________
Bool_t AliVertexingHFUtils::IsCandidateInjected(AliAODRecoDecayHF *cand, AliAODMCHeader *header,TClonesArray *arrayMC){
  /// method to check if a D meson candidate comes from the signal event or from the underlying Hijing event
  /// works only for refilled candidates!
  Int_t nprongs=cand->GetNProngs();
  for(Int_t i=0;i<nprongs;i++){
    AliAODTrack *daugh=(AliAODTrack*)cand->GetDaughter(i);
    if(IsTrackInjected(daugh,header,arrayMC)) return kTRUE;
  }
  return kFALSE;
}
//____________________________________________________________________________
Bool_t AliVertexingHFUtils::IsCandidateInjected(AliAODRecoDecayHF *cand, AliAODEvent* aod, AliAODMCHeader *header,TClonesArray *arrayMC){
  /// method to check if a D meson candidate comes from the signal event or from the underlying Hijing event
  /// works also with not-refilled candidates of reduced AODs
  /// For refilled candidates other IsCandidateInjected method recommended

  Int_t nprongs=cand->GetNProngs();
  for(Int_t i=0;i<nprongs;i++){
    Int_t idDau=cand->GetProngID(i);
    for(Int_t i=0; i<aod->GetNumberOfTracks(); i++) {
      AliAODTrack *daugh=(AliAODTrack*)aod->GetTrack(i);
      if(daugh && (Int_t)daugh->GetID()==idDau){
        if(AliVertexingHFUtils::IsTrackInjected(daugh,header,arrayMC)) return kTRUE;
      }
    }
  }
  return kFALSE;
}
//____________________________________________________________________________
Bool_t AliVertexingHFUtils::HasCascadeCandidateAnyDaughInjected(AliAODRecoCascadeHF *cand, AliAODMCHeader *header,TClonesArray *arrayMC){
  /// method to check if a cascade candidate comes from the signal event or from the underlying Hijing event

  AliAODTrack* bach = cand->GetBachelor();
  if(IsTrackInjected(bach, header, arrayMC)) {
    //    printf("Bachelor is injected, the whole candidate is then injected\n");
    return kTRUE;
  }
  AliAODv0* v0 = cand->Getv0();
  Int_t nprongs = v0->GetNProngs();
  for(Int_t i = 0; i < nprongs; i++){
    AliAODTrack *daugh = (AliAODTrack*)v0->GetDaughter(i);
    if(IsTrackInjected(daugh,header,arrayMC)) {
      //      printf("V0 daughter number %d is injected, the whole candidate is then injected\n", i);
      return kTRUE;
    }
  }
  return kFALSE;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::PreSelectITSUpgrade(TClonesArray* arrayMC, AliAODMCHeader *header, TObjArray aodTracks, Int_t nDaug, Int_t pdgabs, const Int_t *pdgDg){
  /// Preselect function for ITS Upgrade MC's, to make a fast general preselection before filling the HF candidate
  /// Returns 0 when combination of injected+HIJING
  /// Returns 1 when purely injected + matched to MC
  /// Returns 2 when purely injected + not matched
  /// Returns 3 when purely HIJING

  Bool_t injected = kFALSE;
  Bool_t hijing = kFALSE;

  AliAODTrack *track[(const Int_t)nDaug];
  Int_t dgLabels[(const Int_t)nDaug];
  for(Int_t iD=0; iD<nDaug; iD++) {
    track[iD] = (AliAODTrack*)aodTracks.At(iD);
    dgLabels[iD] = track[iD]->GetLabel();

    Bool_t isTrInjected = IsTrackInjected(track[iD],header,arrayMC);
    if(isTrInjected) injected = kTRUE;
    else             hijing = kTRUE;

    //Combination of injected signal + HIJING background
    if(injected && hijing) return 0;
  }

  //Purely HIJING background, no need to MatchToMC
  if(hijing) return 3;

  //Purely injected signal, check by Matching to MC if it is considered as signal
  //  Mimic MatchToMC (can't access as it is protected), code below is the same as:
  //  Int_t MatchedToLabel = AliAODRecoDecay::MatchToMC(pdgabs, arrayMC, dgLabels, nDaug, 0, pdgDg);
  if(injected){
    //Label of mother particle, or -1 when some operation failed.
    Int_t MatchedToLabel = 0;

    Int_t labMom[10]={0,0,0,0,0,0,0,0,0,0};
    Int_t i,j,lab,labMother,pdgMother,pdgPart;
    AliAODMCParticle *part=0;
    AliAODMCParticle *mother=0;
    Bool_t pdgUsed[10]={kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};

    // loop on daughter labels
    for(i=0; i<nDaug; i++) {
      if(MatchedToLabel == -1) continue;

      labMom[i]=-1;
      lab = TMath::Abs(dgLabels[i]);
      if(lab<0){ MatchedToLabel = -1; continue; }

      part = (AliAODMCParticle*)arrayMC->At(lab);
      if(!part){ MatchedToLabel = -1; continue; }

      // check the PDG of the daughter, if requested
      pdgPart=TMath::Abs(part->GetPdgCode());
      for(j=0; j<nDaug; j++) {
        if(!pdgUsed[j] && pdgPart==pdgDg[j]) {
          pdgUsed[j]=kTRUE;
          break;
        }
      }

      mother = part;
      while(mother->GetMother()>=0) {
        labMother=mother->GetMother();
        mother = (AliAODMCParticle*)arrayMC->At(labMother);
        if(!mother) break;

        pdgMother = TMath::Abs(mother->GetPdgCode());
        if(pdgMother==pdgabs) {
          labMom[i]=labMother;
          break;
        } else if(pdgMother>pdgabs || pdgMother<10) {
          break;
        }
      }
      if(labMom[i]==-1) MatchedToLabel = -1; // mother PDG not ok for this daughter
    } // end loop on daughters

    if(MatchedToLabel == 0){
      // check if the candidate is signal
      labMother=labMom[0];
      // all labels have to be the same and !=-1
      for(i=0; i<nDaug; i++) {
        if(labMom[i]==-1)        MatchedToLabel = -1;
        if(labMom[i]!=labMother) MatchedToLabel = -1;
      }

      // check that all daughter PDGs are matched
      for(i=0; i<nDaug; i++) {
        if(pdgUsed[i]==kFALSE) MatchedToLabel = -1;
      }
    }

    if(MatchedToLabel==0) MatchedToLabel = labMother;

    if(MatchedToLabel != -1) return 1; //injected, matched to  MC
    else                     return 2; //injected, not matched to MC
  }

  //Should not reach this, if so check. Return 0 for compilation warning
  printf("AliVertexingHFUtils::PreSelectITSUpgrade: Neither injected, nor HIJING");
  return 0;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckOrigin(AliMCEvent* mcEvent, AliMCParticle *mcPart, Bool_t searchUpToQuark){
  /// checking whether the mother of the particles come from a charm or a bottom quark

  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPart->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  Bool_t isFromB=kFALSE;
  Bool_t isQuarkFound=kFALSE;
  while (mother >=0 ){
    istep++;
    AliMCParticle* mcGranma = (AliMCParticle*)mcEvent->GetTrack(mother);
    if (mcGranma){
      TParticle* partGranma = mcGranma->Particle();
      if(partGranma) pdgGranma = partGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
	isFromB=kTRUE;
      }
      if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
      mother = mcGranma->GetMother();
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
  while (mother >=0 ){
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
Bool_t AliVertexingHFUtils::IsTrackFromCharm(AliAODTrack* tr, TClonesArray* arrayMC){
  /// check if an AOD track originated from a charm hadron decay
  Int_t absLabel=TMath::Abs(tr->GetLabel());
  AliAODMCParticle* mcPart=dynamic_cast<AliAODMCParticle*>(arrayMC->At(absLabel));
  Int_t mother = mcPart->GetMother();
  Int_t istep = 0;
  while (mother >=0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      Int_t abspdgGranma = TMath::Abs(mcGranma->GetPdgCode());
      if ((abspdgGranma==4) ||(abspdgGranma>400 && abspdgGranma<500) || (abspdgGranma>4000 && abspdgGranma<5000)) return kTRUE;
      mother = mcGranma->GetMother();
    }else{
      printf("AliVertexingHFUtils::IsTrackFromCharm: Failed casting the mother particle!");
      break;
    }
  }
  return kFALSE;
}
//____________________________________________________________________________
Bool_t AliVertexingHFUtils::IsTrackFromBeauty(AliAODTrack* tr, TClonesArray* arrayMC){
  /// check if an AOD track originated from a charm hadron decay
  Int_t absLabel=TMath::Abs(tr->GetLabel());
  AliAODMCParticle* mcPart=dynamic_cast<AliAODMCParticle*>(arrayMC->At(absLabel));
  Int_t mother = mcPart->GetMother();
  Int_t istep = 0;
  while (mother >=0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      Int_t abspdgGranma = TMath::Abs(mcGranma->GetPdgCode());
      if ((abspdgGranma==5) ||(abspdgGranma>500 && abspdgGranma<600) || (abspdgGranma>5000 && abspdgGranma<6000)) return kTRUE;
      mother = mcGranma->GetMother();
    }else{
      printf("AliVertexingHFUtils::IsTrackFromBeauty: Failed casting the mother particle!");
      break;
    }
  }
  return kFALSE;
}
//____________________________________________________________________________
Bool_t AliVertexingHFUtils::IsTrackFromHadronDecay(Int_t pdgMoth, AliAODTrack* tr, TClonesArray* arrayMC){
  /// check if an AOD track originated from a charm hadron decay
  Int_t absLabel=TMath::Abs(tr->GetLabel());
  Int_t absPdgMoth=TMath::Abs(pdgMoth);
  AliAODMCParticle* mcPart=dynamic_cast<AliAODMCParticle*>(arrayMC->At(absLabel));
  Int_t mother = mcPart->GetMother();
  Int_t istep = 0;
  while (mother >=0 ){
    istep++;
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma){
      Int_t abspdgGranma = TMath::Abs(mcGranma->GetPdgCode());
      if (abspdgGranma==absPdgMoth) return kTRUE;
      mother = mcGranma->GetMother();
    }else{
      printf("AliVertexingHFUtils::IsTrackFromHadronDecay: Failed casting the mother particle!");
      break;
    }
  }
  return kFALSE;
}
//____________________________________________________________________________
Double_t AliVertexingHFUtils::GetBeautyMotherPt(AliMCEvent* mcEvent, AliMCParticle *mcPart){
  /// get the pt of the beauty hadron (feed-down case), returns negative value for prompt (ESD case)

  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPart->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  while (mother >=0 ){
    istep++;
    AliMCParticle* mcGranma = (AliMCParticle*)mcEvent->GetTrack(mother);
    if (mcGranma){
      TParticle* partGranma = mcGranma->Particle();
      if(partGranma) pdgGranma = partGranma->GetPdgCode();
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
Double_t AliVertexingHFUtils::GetBeautyMotherPt(TClonesArray* arrayMC, AliAODMCParticle *mcPart){
  /// get the pt of the beauty hadron (feed-down case), returns negative value for prompt

  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPart->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma =0;
  while (mother >=0 ){
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
Double_t AliVertexingHFUtils::GetBeautyMotherPtAndPDG(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t &pdgGranma){
  /// get the pt of the beauty hadron (feed-down case), returns negative value for prompt

  pdgGranma = 0;
  Int_t mother = 0;
  mother = mcPart->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma = 0;
  while (mother >=0 ){
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
  pdgGranma = 0;
  return -999.;
}

//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckD0Decay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the D0 decay channel. Returns 1 for the D0->Kpi case, 2 for the D0->Kpipipi case, -1 in other cases

  if(label<0) return -1;
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=421) return -1;

  Int_t nDau=mcPart->GetNDaughters();

  if(nDau==2){
    Int_t daughter0 = mcPart->GetDaughterFirst();
    Int_t daughter1 = mcPart->GetDaughterLast();
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
    if(TMath::Abs(part->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(part->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(part->Pz()-sumPzDau)>0.001) return -2;
    return 1;
  }

  if(nDau==3 || nDau==4){
    Int_t nKaons=0;
    Int_t nPions=0;
    Double_t sumPxDau=0.;
    Double_t sumPyDau=0.;
    Double_t sumPzDau=0.;
    Int_t nFoundKpi=0;
    Int_t labelFirstDau = mcPart->GetDaughterFirst();
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      AliMCParticle* mcDau=(AliMCParticle*)mcEvent->GetTrack(indDau);
      TParticle* dau=mcEvent->Particle(indDau);
      if(!mcDau || !dau) return -1;
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
	Int_t nResDau=mcDau->GetNDaughters();
	if(nResDau!=2) return -1;
	Int_t indFirstResDau=mcDau->GetDaughterFirst();
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
    if(TMath::Abs(part->Px()-sumPxDau)>0.001) return -1;
    if(TMath::Abs(part->Py()-sumPyDau)>0.001) return -1;
    if(TMath::Abs(part->Pz()-sumPzDau)>0.001) return -1;
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
    Int_t daughter0 = mcPart->GetDaughterLabel(0);
    Int_t daughter1 = mcPart->GetDaughterLabel(1);
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
    Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
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
	Int_t indFirstResDau=dau->GetDaughterLabel(0);
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
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=411) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughterFirst();
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
      AliMCParticle* mcDau=(AliMCParticle*)mcEvent->GetTrack(indDau);
      TParticle* dau=mcEvent->Particle(indDau);
      if(!mcDau || !dau) return -1;
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
	Int_t nResDau=mcDau->GetNDaughters();
	if(nResDau!=2) return -1;
	Int_t indFirstResDau=mcDau->GetDaughterFirst();
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
    if(TMath::Abs(part->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(part->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(part->Pz()-sumPzDau)>0.001) return -2;
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
  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
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
	Int_t indFirstResDau=dau->GetDaughterLabel(0);
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
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=411) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughterFirst();
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
      AliMCParticle* mcDau=(AliMCParticle*)mcEvent->GetTrack(indDau);
      TParticle* dau=mcEvent->Particle(indDau);
      if(!mcDau || !dau) return -1;
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
	Int_t nResDau=mcDau->GetNDaughters();
	if(nResDau!=2) return -1;
	Int_t indFirstResDau=mcDau->GetDaughterFirst();
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
    if(TMath::Abs(part->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(part->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(part->Pz()-sumPzDau)>0.001) return -2;
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
  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
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
	Int_t indFirstResDau=dau->GetDaughterLabel(0);
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
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=411) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughterFirst();
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
      AliMCParticle* mcDau=(AliMCParticle*)mcEvent->GetTrack(indDau);
      TParticle* dau=mcEvent->Particle(indDau);
      if(!mcDau || !dau) return -1;
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
	AliMCParticle* mcV0=mcDau;
	if(codeV0==311){
	  Int_t nK0Dau=mcDau->GetNDaughters();
	  if(nK0Dau!=1) return -1;
	  Int_t indK0s=mcDau->GetDaughterFirst();
	  if(indK0s<0) return -1;
	  v0=mcEvent->Particle(indK0s);
	  mcV0=(AliMCParticle*)mcEvent->GetTrack(indK0s);
	  if(!v0 || !mcV0) return -1;
	  Int_t pdgK0sl=v0->GetPdgCode();
	  codeV0=TMath::Abs(pdgK0sl);
	}
	Int_t nV0Dau=mcV0->GetNDaughters();
	if(nV0Dau!=2) return -1;
	Int_t indFirstV0Dau=mcV0->GetDaughterFirst();
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
    if(TMath::Abs(part->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(part->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(part->Pz()-sumPzDau)>0.001) return -2;
    if(codeV0==310) return 1;
  }
  return -1;

}


//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDsDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Ds decay channel. Returns 1 for Ds->phipi->KKpi, 2 for Ds->K0*K->KKpi, 3 for the non-resonant case, 4 for Ds->f0pi->KKpi

  if(label<0) return -1;
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=431) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughterFirst();
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;
  Bool_t isPhi=kFALSE;
  Bool_t isk0st=kFALSE;
  Bool_t isf0=kFALSE;

  if(nDau==3 || nDau==2){
    for(Int_t iDau=0; iDau<nDau; iDau++){
      Int_t indDau = labelFirstDau+iDau;
      if(indDau<0) return -1;
      AliMCParticle* mcDau=(AliMCParticle*)mcEvent->GetTrack(indDau);
      TParticle* dau=mcEvent->Particle(indDau);
      if(!mcDau || !dau) return -1;
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
      }else if(TMath::Abs(pdgdau)==313 || TMath::Abs(pdgdau)==333 || TMath::Abs(pdgdau)==9010221){
	if(TMath::Abs(pdgdau)==313) isk0st=kTRUE;
	if(TMath::Abs(pdgdau)==333) isPhi=kTRUE;
        if(TMath::Abs(pdgdau)==9010221) isf0=kTRUE;
	Int_t nResDau=mcDau->GetNDaughters();
	if(nResDau!=2) return -1;
	Int_t indFirstResDau=mcDau->GetDaughterFirst();
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
    if(TMath::Abs(part->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(part->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(part->Pz()-sumPzDau)>0.001) return -2;
    if(nDau==3) return 3;
    else if(nDau==2){
      if(isk0st) return 2;
      if(isPhi) return 1;
      if(isf0) return 4;
    }
  }

  return -1;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDsDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  /// Checks the Ds decay channel. Returns 1 for Ds->phipi->KKpi, 2 for Ds->K0*K->KKpi, 3 for the non-resonant case, 4 for Ds->f0pi->KKpi

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=431) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;
  Bool_t isPhi=kFALSE;
  Bool_t isk0st=kFALSE;
  Bool_t isf0=kFALSE;

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
      }else if(TMath::Abs(pdgdau)==313 || TMath::Abs(pdgdau)==333 || TMath::Abs(pdgdau)==9010221){
	if(TMath::Abs(pdgdau)==313) isk0st=kTRUE;
	if(TMath::Abs(pdgdau)==333) isPhi=kTRUE;
        if(TMath::Abs(pdgdau)==9010221) isf0=kTRUE;
	Int_t nResDau=dau->GetNDaughters();
	if(nResDau!=2) return -1;
	Int_t indFirstResDau=dau->GetDaughterLabel(0);
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
      if(isf0) return 4;
    }
  }

  return -1;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDsK0sKDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Ds->K0s+S decay channel. Returns 1 in case of success, otherwise -1

  if(label<0) return -1;
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=431) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughterFirst();
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
      AliMCParticle* mcDau=(AliMCParticle*)mcEvent->GetTrack(indDau);
      TParticle* dau=mcEvent->Particle(indDau);
      if(!mcDau || !dau) return -1;
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
	AliMCParticle* mcV0=mcDau;
	if(codeV0==311){
	  Int_t nK0Dau=mcDau->GetNDaughters();
	  if(nK0Dau!=1) return -1;
	  Int_t indK0s=mcDau->GetDaughterFirst();
	  if(indK0s<0) return -1;
	  v0=mcEvent->Particle(indK0s);
	  mcV0=(AliMCParticle*)mcEvent->GetTrack(indK0s);
	  if(!v0 || !mcV0) return -1;
	  Int_t pdgK0sl=v0->GetPdgCode();
	  codeV0=TMath::Abs(pdgK0sl);
	}
	Int_t nV0Dau=mcV0->GetNDaughters();
	if(nV0Dau!=2) return -1;
	Int_t indFirstV0Dau=mcV0->GetDaughterFirst();
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
    if(TMath::Abs(part->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(part->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(part->Pz()-sumPzDau)>0.001) return -2;
    if(codeV0==310) return 1;
  }
  return -1;

}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDstarDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Dstar decay channel. Returns 1 for D*->D0pi->Kpipi, -1 in other cases

  if(label<0) return -1;
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=413) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughterFirst();
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;

  for(Int_t iDau=0; iDau<nDau; iDau++){
    Int_t indDau = labelFirstDau+iDau;
    if(indDau<0) return -1;
    AliMCParticle* mcDau=(AliMCParticle*)mcEvent->GetTrack(indDau);
    TParticle* dau=mcEvent->Particle(indDau);
    if(!mcDau || !dau) return -1;
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==421){
      Int_t nResDau=mcDau->GetNDaughters();
      if(nResDau!=2) return -1;
      Int_t indFirstResDau=mcDau->GetDaughterFirst();
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
  if(TMath::Abs(part->Px()-sumPxDau)>0.001) return -2;
  if(TMath::Abs(part->Py()-sumPyDau)>0.001) return -2;
  if(TMath::Abs(part->Pz()-sumPzDau)>0.001) return -2;
  return 1;

}

//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckDstarDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  /// Checks the Dstar decay channel. Returns 1 for D*->D0pi->Kpipi, -1 in other cases

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=413) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
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
      Int_t indFirstResDau=dau->GetDaughterLabel(0);
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
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=4122) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughterFirst();
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
      AliMCParticle* mcDau=(AliMCParticle*)mcEvent->GetTrack(indDau);
      TParticle* dau=mcEvent->Particle(indDau);
      if(!mcDau || !dau) return -1;
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
	Int_t nResDau=mcDau->GetNDaughters();
	if(nResDau!=2) return -1;
	Int_t indFirstResDau=mcDau->GetDaughterFirst();
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
    if(TMath::Abs(part->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(part->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(part->Pz()-sumPzDau)>0.001) return -2;
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
  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
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
	Int_t indFirstResDau=dau->GetDaughterLabel(0);
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
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=4122) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughterFirst();
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
      AliMCParticle* mcDau=(AliMCParticle*)mcEvent->GetTrack(indDau);
      TParticle* dau=mcEvent->Particle(indDau);
      if(!mcDau || !dau) return -1;
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
	AliMCParticle* mcV0=mcDau;
	if(codeV0==311){
	  Int_t nK0Dau=mcDau->GetNDaughters();
	  if(nK0Dau!=1) return -1;
	  Int_t indK0s=mcDau->GetDaughterFirst();
	  if(indK0s<0) return -1;
	  v0=mcEvent->Particle(indK0s);
	  mcV0=(AliMCParticle*)mcEvent->GetTrack(indK0s);
	  if(!v0 || !mcV0) return -1;
	  Int_t pdgK0sl=v0->GetPdgCode();
	  codeV0=TMath::Abs(pdgK0sl);
	}
	Int_t nV0Dau=mcV0->GetNDaughters();
	if(nV0Dau!=2) return -1;
	Int_t indFirstV0Dau=mcV0->GetDaughterFirst();
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
    if(TMath::Abs(part->Px()-sumPxDau)>0.001) return -2;
    if(TMath::Abs(part->Py()-sumPyDau)>0.001) return -2;
    if(TMath::Abs(part->Pz()-sumPzDau)>0.001) return -2;
    if(codeV0==310) return 1;
    else if(codeV0==3122) return 2;
  }
  return -1;

}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckLcV0bachelorDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  /// Checks the Lc->V0+bachelor decay channel. Returns 1 for pK0s, 2 for piLambda, 3 for pK0l -1 in other cases

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=4122) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
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
      AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
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
  AliAODMCParticle* v0=dau;
	if(codeV0==311){
	  Int_t nK0Dau=dau->GetNDaughters();
	  if(nK0Dau!=1) return -1;
	  Int_t indK0s=dau->GetDaughterLabel(0);
	  if(indK0s<0) return -1;
	  v0=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indK0s));
	  if(!v0) return -1;
	  Int_t pdgK0sl=v0->GetPdgCode();
	  codeV0=TMath::Abs(pdgK0sl);
	}
	Int_t nV0Dau=v0->GetNDaughters();
	if(nV0Dau!=2) return -1;
	Int_t indFirstV0Dau=v0->GetDaughterLabel(0);
	for(Int_t v0Dau=0; v0Dau<2; v0Dau++){
	  Int_t indV0Dau=indFirstV0Dau+v0Dau;
	  if(indV0Dau<0) return -1;
    AliAODMCParticle* v0dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indV0Dau));
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
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=4232) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=3) return -1;

  Int_t labelFirstDau = mcPart->GetDaughterFirst();
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
  if(TMath::Abs(part->Px()-sumPxDau)>0.001) return -2;
  if(TMath::Abs(part->Py()-sumPyDau)>0.001) return -2;
  if(TMath::Abs(part->Pz()-sumPzDau)>0.001) return -2;
  return 1;

}

//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckBplusDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Bplus decay channel. Returns 1 for Bplus->D0pi->Kpipi, -1 in other cases
  /// If rejected by momentum conservation check, return (-1*decay - 1) (to allow checks at task level)
  ///
  /// NB: Loosened cut on mom. conserv. (needed because of small issue in ITS Upgrade productions)

  if(label<0) return -1;
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=521) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughterFirst();
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;

  for(Int_t iDau=0; iDau<nDau; iDau++){
    Int_t indDau = labelFirstDau+iDau;
    if(indDau<0) return -1;
    AliMCParticle* mcDau=(AliMCParticle*)mcEvent->GetTrack(indDau);
    TParticle* dau=mcEvent->Particle(indDau);
    if(!mcDau || !dau) return -1;
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==421){
      Int_t nResDau=mcDau->GetNDaughters();
      if(nResDau!=2) return -1;
      Int_t indFirstResDau=mcDau->GetDaughterFirst();
      for(Int_t resDau=0; resDau<2; resDau++){
        Int_t indResDau=indFirstResDau+resDau;
        if(indResDau<0) return -1;
        TParticle* resdau=mcEvent->Particle(indResDau);
        if(!resdau) return -1;
        Int_t pdgresdau=resdau->GetPdgCode();
        if(TMath::Abs(pdgresdau)==321){
          if(pdgD*pdgresdau<0) return -1;
          sumPxDau+=resdau->Px();
          sumPyDau+=resdau->Py();
          sumPzDau+=resdau->Pz();
          nKaons++;
          arrayDauLab[nFoundKpi++]=indResDau;
          if(nFoundKpi>3) return -1;
        }
        if(TMath::Abs(pdgresdau)==211){
          if(pdgD*pdgresdau>0) return -1;
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
  //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
  //Fix implemented, loosening cut from 0.001 to 0.1. If >0.1, (-1*decay - 1) is returned.
  if(TMath::Abs(part->Px()-sumPxDau)>0.1) return -2;
  if(TMath::Abs(part->Py()-sumPyDau)>0.1) return -2;
  if(TMath::Abs(part->Pz()-sumPzDau)>0.1) return -2;
  return 1;

}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckBplusDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  /// Checks the Bplus decay channel. Returns 1 for Bplus->D0pi->Kpipi, -1 in other cases
  /// If rejected by momentum conservation check, return (-1*decay - 1) (to allow checks at task level)
  ///
  /// NB: Loosened cut on mom. conserv. (needed because of small issue in ITS Upgrade productions)

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=521) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
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
      Int_t indFirstResDau=dau->GetDaughterLabel(0);
      for(Int_t resDau=0; resDau<2; resDau++){
        Int_t indResDau=indFirstResDau+resDau;
        if(indResDau<0) return -1;
        AliAODMCParticle* resdau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indResDau));
        if(!resdau) return -1;
        Int_t pdgresdau=resdau->GetPdgCode();
        if(TMath::Abs(pdgresdau)==321){
          if(pdgD*pdgresdau<0) return -1;
          sumPxDau+=resdau->Px();
          sumPyDau+=resdau->Py();
          sumPzDau+=resdau->Pz();
          nKaons++;
          arrayDauLab[nFoundKpi++]=indResDau;
          if(nFoundKpi>3) return -1;
        }
        if(TMath::Abs(pdgresdau)==211){
          if(pdgD*pdgresdau>0) return -1;
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
  //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
  //Fix implemented, loosening cut from 0.001 to 0.1. If >0.1, (-1*decay - 1) is returned.
  if(TMath::Abs(mcPart->Px()-sumPxDau)>0.1) return -2;
  if(TMath::Abs(mcPart->Py()-sumPyDau)>0.1) return -2;
  if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.1) return -2;
  return 1;

}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckB0toDminuspiDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Bs decay channel. Returns >= 1 for B0->Dminuspi->Kpipipi, <0 in other cases
  /// Returns 1 for the non-resonant case, 2 for the resonant case, -1 in other cases
  /// If rejected by momentum conservation check, return (-1*decay - 1) (to allow checks at task level)
  ///
  /// NB: Loosened cut on mom. conserv. (needed because of small issue in ITS Upgrade productions)

  if(label<0) return -1;
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=511) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughterFirst();
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;
  Int_t decayB0 = -1;

  for(Int_t iDau=0; iDau<nDau; iDau++){
    Int_t indDau = labelFirstDau+iDau;
    if(indDau<0) return -1;
    TParticle* dau=mcEvent->Particle(indDau);
    if(!dau) return -1;
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==411){
      /// Checks the Dplus decay channel. Returns 1 for the non-resonant case, 2 for the resonant case, -1 in other cases
      Int_t labDauDplus[3] = {-1,-1,-1};
      Int_t decayDplus = CheckDplusDecay(mcEvent, indDau, labDauDplus);

      //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
      //Fix implemented, to still select correct Dplus decay
      if(decayDplus==-2){
        AliMCParticle* mcDplus = (AliMCParticle*)mcEvent->GetTrack(indDau);
        Int_t labelFirstDauDplus = mcDplus->GetDaughterFirst();
        if(mcDplus->GetNDaughters() > 1){
          TParticle* dauDplus1 = mcEvent->Particle(labelFirstDauDplus);
          TParticle* dauDplus2 = mcEvent->Particle(labelFirstDauDplus+1);
          if(dauDplus1 && dauDplus2){
            Int_t pdgdauDplus1=dauDplus1->GetPdgCode();
            Int_t pdgdauDplus2=dauDplus2->GetPdgCode();
            if(TMath::Abs(pdgdauDplus1)==313 || TMath::Abs(pdgdauDplus2)==313) decayDplus=2;
          }
        }
      }

      if (decayDplus < 0 || labDauDplus[0] == -1) return -1;
      decayB0 = decayDplus;
      nPions+=2;
      nKaons++;
      for(Int_t iDplus = 0; iDplus < 3; iDplus++){
        TParticle* dauDplus=mcEvent->Particle(labDauDplus[iDplus]);
        sumPxDau+=dauDplus->Px();
        sumPyDau+=dauDplus->Py();
        sumPzDau+=dauDplus->Pz();
        arrayDauLab[nFoundKpi++]=labDauDplus[iDplus];
      }
      if(nFoundKpi>4) return -1;
    }else if(TMath::Abs(pdgdau)==211){
      //Temp fix for B0->D-+pi- and B0bar->D+pi- decays in ITS upgrade productions
      //if(pdgD*pdgdau>0) return -1;
      nPions++;
      sumPxDau+=dau->Px();
      sumPyDau+=dau->Py();
      sumPzDau+=dau->Pz();
      arrayDauLab[nFoundKpi++]=indDau;
      if(nFoundKpi>4) return -1;
    }
  }

  if(nPions!=3) return -1;
  if(nKaons!=1) return -1;
  //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
  //Fix implemented, loosening cut from 0.001 to 0.1. If >0.1, (-1*decay - 1) is returned.
  Int_t decayB0temp = decayB0;
  if(TMath::Abs(part->Px()-sumPxDau)>0.1) decayB0temp = -1*decayB0 - 1;
  if(TMath::Abs(part->Py()-sumPyDau)>0.1) decayB0temp = -1*decayB0 - 1;
  if(TMath::Abs(part->Pz()-sumPzDau)>0.1) decayB0temp = -1*decayB0 - 1;
  decayB0 = decayB0temp;
  return decayB0;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckB0toDminuspiDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  /// Checks the Bs decay channel. Returns >= 1 for B0->Dminuspi->Kpipipi, <0 in other cases
  /// Returns 1 for the non-resonant case, 2 for the resonant case, -1 in other cases
  /// If rejected by momentum conservation check, return (-1*decay - 1) (to allow checks at task level)
  ///
  /// NB: Loosened cut on mom. conserv. (needed because of small issue in ITS Upgrade productions)

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=511) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;
  Int_t decayB0 = -1;

  for(Int_t iDau=0; iDau<nDau; iDau++){
    Int_t indDau = labelFirstDau+iDau;
    if(indDau<0) return -1;
    AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
    if(!dau) return -1;
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==411){
      /// Checks the Dplus decay channel. Returns 1 for the non-resonant case, 2 for the resonant case, -1 in other cases
      Int_t labDauDplus[3] = {-1,-1,-1};
      Int_t decayDplus = CheckDplusDecay(arrayMC, dau, labDauDplus);

      //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
      //Fix implemented, to still select correct Dplus decay
      if(decayDplus==-2){
        Int_t labelFirstDauDplus = dau->GetDaughterLabel(0);
        if(dau->GetNDaughters() > 1){
          AliAODMCParticle* dauDplus1=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labelFirstDauDplus));
          AliAODMCParticle* dauDplus2=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labelFirstDauDplus+1));
          if(dauDplus1 && dauDplus2){
            Int_t pdgdauDplus1=dauDplus1->GetPdgCode();
            Int_t pdgdauDplus2=dauDplus2->GetPdgCode();
            if(TMath::Abs(pdgdauDplus1)==313 || TMath::Abs(pdgdauDplus2)==313) decayDplus=2;
          }
        }
      }

      if (decayDplus < 0 || labDauDplus[0] == -1) return -1;
      decayB0 = decayDplus;
      nPions+=2;
      nKaons++;
      for(Int_t iDplus = 0; iDplus < 3; iDplus++){
        AliAODMCParticle* dauDplus=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDauDplus[iDplus]));
        sumPxDau+=dauDplus->Px();
        sumPyDau+=dauDplus->Py();
        sumPzDau+=dauDplus->Pz();
        arrayDauLab[nFoundKpi++]=labDauDplus[iDplus];
      }
      if(nFoundKpi>4) return -1;
    }else if(TMath::Abs(pdgdau)==211){
      //Temp fix for B0->D-+pi- and B0bar->D+pi- decays in ITS upgrade productions
      //if(pdgD*pdgdau>0) return -1;
      nPions++;
      sumPxDau+=dau->Px();
      sumPyDau+=dau->Py();
      sumPzDau+=dau->Pz();
      arrayDauLab[nFoundKpi++]=indDau;
      if(nFoundKpi>4) return -1;
    }
  }

  if(nPions!=3) return -1;
  if(nKaons!=1) return -1;
  //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
  //Fix implemented, loosening cut from 0.001 to 0.1. If >0.1, (-1*decay - 1) is returned.
  Int_t decayB0temp = decayB0;
  if(TMath::Abs(mcPart->Px()-sumPxDau)>0.1) decayB0temp = -1*decayB0 - 1;
  if(TMath::Abs(mcPart->Py()-sumPyDau)>0.1) decayB0temp = -1*decayB0 - 1;
  if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.1) decayB0temp = -1*decayB0 - 1;
  decayB0 = decayB0temp;
  return decayB0;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckResoToDplusK0SDecay(TClonesArray *arrayMC, AliAODMCParticle *mcPart, Int_t *arrayDauLab) {
  /// Checks the D* -> D+(-> K-pi+pi+) K0S(->pi+pi-) decay channel.
  /// Returns pdg of the resonance, -1 in other cases

  Int_t pdgD = mcPart->GetPdgCode();
  if (TMath::Abs(pdgD) != 435 && TMath::Abs(pdgD) != 10433)
    return -1;

  Int_t nDau = mcPart->GetNDaughters();
  if (nDau != 2)
    return -1;

  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
  Int_t nKaons = 0;
  Int_t nPions = 0;
  Double_t sumPxDau = 0.;
  Double_t sumPyDau = 0.;
  Double_t sumPzDau = 0.;
  Int_t nFoundKpi = 0;

  for (Int_t iDau = 0; iDau < nDau; iDau++) {
    Int_t indDau = labelFirstDau + iDau;
    if (indDau < 0)
      return -1;
    AliAODMCParticle *dau = dynamic_cast<AliAODMCParticle *>(arrayMC->At(indDau));
    if (!dau)
      return -1;
    Int_t pdgdau = dau->GetPdgCode();
    if (TMath::Abs(pdgdau) == 411) {
      /// Checks the Dplus decay channel. Returns 1 for the non-resonant case, 2
      /// for the resonant case, -1 in other cases
      Int_t labDauDplus[3] = {-1, -1, -1};
      Int_t decayDplus = CheckDplusDecay(arrayMC, dau, labDauDplus);
      if (decayDplus < 0 || labDauDplus[0] == -1)
        return -1;

      nPions += 2;
      nKaons++;
      for (Int_t iDplus = 0; iDplus < 3; iDplus++) {
        AliAODMCParticle *dauDplus = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labDauDplus[iDplus]));
        sumPxDau += dauDplus->Px();
        sumPyDau += dauDplus->Py();
        sumPzDau += dauDplus->Pz();
        arrayDauLab[nFoundKpi++] = labDauDplus[iDplus];
      }
      if (nFoundKpi > 4)
        return -1;
    } else if (TMath::Abs(pdgdau) == 311) {
      AliAODMCParticle *v0 = dau;
      Int_t nK0Dau = dau->GetNDaughters();
      if (nK0Dau != 1)
        return -1;
      Int_t indK0s = dau->GetDaughterLabel(0);
      if (indK0s < 0)
        return -1;
      v0 = dynamic_cast<AliAODMCParticle *>(arrayMC->At(indK0s));
      if (!v0)
        return -1;
      Int_t nV0Dau = v0->GetNDaughters();
      if (nV0Dau != 2)
        return -1;

      Int_t indFirstV0Dau = v0->GetDaughterLabel(0);
      for (Int_t v0Dau = 0; v0Dau < 2; v0Dau++) {
        Int_t indV0Dau = indFirstV0Dau + v0Dau;
        if (indV0Dau < 0)
          return -1;
        AliAODMCParticle *v0dau = dynamic_cast<AliAODMCParticle *>(arrayMC->At(indV0Dau));
        if (!v0dau)
          return -1;
        Int_t pdgv0dau = v0dau->GetPdgCode();
        if (TMath::Abs(pdgv0dau) == 211) {
          sumPxDau += v0dau->Px();
          sumPyDau += v0dau->Py();
          sumPzDau += v0dau->Pz();
          nPions++;
          arrayDauLab[nFoundKpi++] = indV0Dau;
          if (nFoundKpi > 5)
            return -1;
        }
      }
    }
  }

  if (nPions != 4)
    return -1;
  if (nKaons != 1)
    return -1;

  if (TMath::Abs(mcPart->Px() - sumPxDau) > 0.1)
    return -1;
  if (TMath::Abs(mcPart->Py() - sumPyDau) > 0.1)
    return -1;
  if (TMath::Abs(mcPart->Pz() - sumPzDau) > 0.1)
    return -1;

  return TMath::Abs(pdgD);
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckResoToDstarK0SDecay(TClonesArray *arrayMC, AliAODMCParticle *mcPart, Int_t *arrayDauLab) {
  /// Checks the D* -> D*+(-> D0 pi+ -> K-pi+pi+) K0S(->pi+pi-) decay channel.
  /// Returns pdg of the resonance, -1 in other cases

  Int_t pdgD = mcPart->GetPdgCode();
  if (TMath::Abs(pdgD) != 10433 && TMath::Abs(pdgD) != 435)
    return -1;

  Int_t nDau = mcPart->GetNDaughters();
  if (nDau != 2)
    return -1;

  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
  Int_t nKaons = 0;
  Int_t nPions = 0;
  Double_t sumPxDau = 0.;
  Double_t sumPyDau = 0.;
  Double_t sumPzDau = 0.;
  Int_t nFoundKpi = 0;

  for (Int_t iDau = 0; iDau < nDau; iDau++) {
    Int_t indDau = labelFirstDau + iDau;
    if (indDau < 0)
      return -1;
    AliAODMCParticle *dau = dynamic_cast<AliAODMCParticle *>(arrayMC->At(indDau));
    if (!dau)
      return -1;
    Int_t pdgdau = dau->GetPdgCode();
    if (TMath::Abs(pdgdau) == 413) {
      /// Checks the Dstar decay channel. Returns 1 for the D0pi decay, -1 in other cases
      Int_t labDauDstar[3] = {-1, -1, -1};
      Int_t decayDstar = CheckDstarDecay(arrayMC, dau, labDauDstar);
      if (decayDstar < 0 || labDauDstar[0] == -1)
        return -1;

      nPions += 2;
      nKaons++;
      for (Int_t iDstar = 0; iDstar < 3; iDstar++) {
        AliAODMCParticle *dauDstar = dynamic_cast<AliAODMCParticle *>(arrayMC->At(labDauDstar[iDstar]));
        sumPxDau += dauDstar->Px();
        sumPyDau += dauDstar->Py();
        sumPzDau += dauDstar->Pz();
        arrayDauLab[nFoundKpi++] = labDauDstar[iDstar];
      }
      if (nFoundKpi > 4)
        return -1;
    } else if (TMath::Abs(pdgdau) == 311) {
      AliAODMCParticle *v0 = dau;
      Int_t nK0Dau = dau->GetNDaughters();
      if (nK0Dau != 1)
        return -1;
      Int_t indK0s = dau->GetDaughterLabel(0);
      if (indK0s < 0)
        return -1;
      v0 = dynamic_cast<AliAODMCParticle *>(arrayMC->At(indK0s));
      if (!v0)
        return -1;
      Int_t nV0Dau = v0->GetNDaughters();
      if (nV0Dau != 2)
        return -1;

      Int_t indFirstV0Dau = v0->GetDaughterLabel(0);
      for (Int_t v0Dau = 0; v0Dau < 2; v0Dau++) {
        Int_t indV0Dau = indFirstV0Dau + v0Dau;
        if (indV0Dau < 0)
          return -1;
        AliAODMCParticle *v0dau = dynamic_cast<AliAODMCParticle *>(arrayMC->At(indV0Dau));
        if (!v0dau)
          return -1;
        Int_t pdgv0dau = v0dau->GetPdgCode();
        if (TMath::Abs(pdgv0dau) == 211) {
          sumPxDau += v0dau->Px();
          sumPyDau += v0dau->Py();
          sumPzDau += v0dau->Pz();
          nPions++;
          arrayDauLab[nFoundKpi++] = indV0Dau;
          if (nFoundKpi > 5)
            return -1;
        }
      }
    }
  }

  if (nPions != 4)
    return -1;
  if (nKaons != 1)
    return -1;

  if (TMath::Abs(mcPart->Px() - sumPxDau) > 0.1)
    return -1;
  if (TMath::Abs(mcPart->Py() - sumPyDau) > 0.1)
    return -1;
  if (TMath::Abs(mcPart->Pz() - sumPzDau) > 0.1)
    return -1;

  return TMath::Abs(pdgD);
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckBsDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab, Bool_t ITS2UpgradeProd){
  /// Checks the Bs decay channel. Returns >= 1 for Bs->Dspi->KKpipi, <0 in other cases
  /// Returns 1 for Ds->phipi->KKpi, 2 for Ds->K0*K->KKpi, 3 for the non-resonant case, 4 for Ds->f0pi->KKpi
  /// If rejected by momentum conservation check, return (-1*decay - 1) (to allow checks at task level)
  ///
  /// NB: Loosened cut on mom. conserv. (needed because of small issue in ITS Upgrade productions)

  if(label<0) return -1;
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=531) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughterFirst();
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;
  Int_t decayBs = -1;

  for(Int_t iDau=0; iDau<nDau; iDau++){
    Int_t indDau = labelFirstDau+iDau;
    if(indDau<0) return -1;
    TParticle* dau=mcEvent->Particle(indDau);
    if(!dau) return -1;
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==431){
      //Returns 1 for Ds->phipi->KKpi, 2 for Ds->K0*K->KKpi, 3 for the non-resonant case, 4 for Ds->f0pi->KKpi
      Int_t labDauDs[3] = {-1,-1,-1};
      Int_t decayDs = CheckDsDecay(mcEvent, indDau, labDauDs);

      //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
      //Fix implemented, to still select correct Ds decay
      if(decayDs==-2){
        AliMCParticle* mcDs = (AliMCParticle*)mcEvent->GetTrack(indDau);
        Int_t labelFirstDauDs = mcDs->GetDaughterFirst();
        if(mcDs->GetNDaughters() > 1){
          TParticle* dauDs1 = mcEvent->Particle(labelFirstDauDs);
          TParticle* dauDs2 = mcEvent->Particle(labelFirstDauDs+1);
          if(dauDs1 && dauDs2){
            Int_t pdgdauDs1=dauDs1->GetPdgCode();
            Int_t pdgdauDs2=dauDs2->GetPdgCode();
            if(TMath::Abs(pdgdauDs1)==313 || TMath::Abs(pdgdauDs2)==313) decayDs=2;
            else if(TMath::Abs(pdgdauDs1)==333 || TMath::Abs(pdgdauDs2)==333) decayDs=1;
            else if(TMath::Abs(pdgdauDs1)==9010221 || TMath::Abs(pdgdauDs2)==9010221) decayDs=4;
            else decayDs=3;
          }
        }
      }

      //In ITS2 Upgrade production, phi not "stored", so decay read as non-resonant
      if(decayDs==3 && ITS2UpgradeProd){
        TParticle* dauK1 = mcEvent->Particle(labDauDs[0]);
        TParticle* dauK2 = mcEvent->Particle(labDauDs[1]);
        TParticle* dauK3 = mcEvent->Particle(labDauDs[2]);

        TLorentzVector vK1, vK2, vKK;
        if(TMath::Abs(dauK1->GetPdgCode())==321 && TMath::Abs(dauK2->GetPdgCode())==321){
          dauK1->Momentum(vK1);
          dauK2->Momentum(vK2);
        } else if(TMath::Abs(dauK1->GetPdgCode())==321 && TMath::Abs(dauK3->GetPdgCode())==321){
          dauK1->Momentum(vK1);
          dauK3->Momentum(vK2);
        } else {
          dauK2->Momentum(vK1);
          dauK3->Momentum(vK2);
        }
        vKK = vK1 + vK2;
        //Small window around phi-mass, tag as Ds->phipi->KKpi if inside
        if(vKK.M() > 1.00 && vKK.M() < 1.04) decayDs = 1;
      }

      if (decayDs < 0 || labDauDs[0] == -1) return -1;
      decayBs = decayDs;
      nPions++;
      nKaons+=2;
      for(Int_t iDs = 0; iDs < 3; iDs++){
        TParticle* dauDs=mcEvent->Particle(labDauDs[iDs]);
        sumPxDau+=dauDs->Px();
        sumPyDau+=dauDs->Py();
        sumPzDau+=dauDs->Pz();
        arrayDauLab[nFoundKpi++]=labDauDs[iDs];
      }
      if(nFoundKpi>4) return -1;
    }else if(TMath::Abs(pdgdau)==211){
      //Temp fix for Bs->Ds+pi- and Bs->Ds-pi+ decays in ITS upgrade productions
      //if(pdgD*pdgdau>0) return -1;
      nPions++;
      sumPxDau+=dau->Px();
      sumPyDau+=dau->Py();
      sumPzDau+=dau->Pz();
      arrayDauLab[nFoundKpi++]=indDau;
      if(nFoundKpi>4) return -1;
    }
  }

  if(nPions!=2) return -1;
  if(nKaons!=2) return -1;
  //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
  //Fix implemented, loosening cut from 0.001 to 0.1. If >0.1, (-1*decay - 1) is returned.
  Int_t decayBstemp = decayBs;
  if(TMath::Abs(part->Px()-sumPxDau)>0.1) decayBstemp = -1*decayBs - 1;
  if(TMath::Abs(part->Py()-sumPyDau)>0.1) decayBstemp = -1*decayBs - 1;
  if(TMath::Abs(part->Pz()-sumPzDau)>0.1) decayBstemp = -1*decayBs - 1;
  decayBs = decayBstemp;
  return decayBs;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckBsDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab, Bool_t ITS2UpgradeProd){
  /// Checks the Bs decay channel. Returns >= 1 for Bs->Dspi->KKpipi, <0 in other cases
  /// Returns 1 for Ds->phipi->KKpi, 2 for Ds->K0*K->KKpi, 3 for the non-resonant case, 4 for Ds->f0pi->KKpi
  /// If rejected by momentum conservation check, return (-1*decay - 1) (to allow checks at task level)
  ///
  /// NB: Loosened cut on mom. conserv. (needed because of small issue in ITS Upgrade productions)

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=531) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundKpi=0;
  Int_t decayBs = -1;

  for(Int_t iDau=0; iDau<nDau; iDau++){
    Int_t indDau = labelFirstDau+iDau;
    if(indDau<0) return -1;
    AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
    if(!dau) return -1;
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==431){
      //Returns 1 for Ds->phipi->KKpi, 2 for Ds->K0*K->KKpi, 3 for the non-resonant case, 4 for Ds->f0pi->KKpi
      Int_t labDauDs[3] = {-1,-1,-1};
      Int_t decayDs = CheckDsDecay(arrayMC, dau, labDauDs);

      //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
      //Fix implemented, to still select correct Ds decay
      if(decayDs==-2){
        Int_t labelFirstDauDs = dau->GetDaughterLabel(0);
        if(dau->GetNDaughters() > 1){
          AliAODMCParticle* dauDs1=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labelFirstDauDs));
          AliAODMCParticle* dauDs2=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labelFirstDauDs+1));
          if(dauDs1 && dauDs2){
            Int_t pdgdauDs1=dauDs1->GetPdgCode();
            Int_t pdgdauDs2=dauDs2->GetPdgCode();
            if(TMath::Abs(pdgdauDs1)==313 || TMath::Abs(pdgdauDs2)==313) decayDs=2;
            else if(TMath::Abs(pdgdauDs1)==333 || TMath::Abs(pdgdauDs2)==333) decayDs=1;
            else if(TMath::Abs(pdgdauDs1)==9010221 || TMath::Abs(pdgdauDs2)==9010221) decayDs=4;
            else decayDs=3;
          }
        }
      }

      //In ITS2 Upgrade production, phi not "stored", so decay read as non-resonant
      if(decayDs==3 && ITS2UpgradeProd){
        AliAODMCParticle* dauK1=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDauDs[0]));
        AliAODMCParticle* dauK2=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDauDs[1]));
        AliAODMCParticle* dauK3=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDauDs[2]));

        TLorentzVector vK1, vK2, vKK;
        if(TMath::Abs(dauK1->GetPdgCode())==321 && TMath::Abs(dauK2->GetPdgCode())==321){
          dauK1->Momentum(vK1);
          dauK2->Momentum(vK2);
        } else if(TMath::Abs(dauK1->GetPdgCode())==321 && TMath::Abs(dauK3->GetPdgCode())==321){
          dauK1->Momentum(vK1);
          dauK3->Momentum(vK2);
        } else {
          dauK2->Momentum(vK1);
          dauK3->Momentum(vK2);
        }
        vKK = vK1 + vK2;
        //Small window around phi-mass, tag as Ds->phipi->KKpi if inside
        if(vKK.M() > 1.00 && vKK.M() < 1.04) decayDs = 1;
      }

      if (decayDs < 0 || labDauDs[0] == -1) return -1;
      decayBs = decayDs;
      nPions++;
      nKaons+=2;
      for(Int_t iDs = 0; iDs < 3; iDs++){
        AliAODMCParticle* dauDs=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDauDs[iDs]));
        sumPxDau+=dauDs->Px();
        sumPyDau+=dauDs->Py();
        sumPzDau+=dauDs->Pz();
        arrayDauLab[nFoundKpi++]=labDauDs[iDs];
      }
      if(nFoundKpi>4) return -1;
    }else if(TMath::Abs(pdgdau)==211){
      //Temp fix for Bs->Ds+pi- and Bs->Ds-pi+ decays in ITS upgrade productions
      //if(pdgD*pdgdau>0) return -1;
      nPions++;
      sumPxDau+=dau->Px();
      sumPyDau+=dau->Py();
      sumPzDau+=dau->Pz();
      arrayDauLab[nFoundKpi++]=indDau;
      if(nFoundKpi>4) return -1;
    }
  }

  if(nPions!=2) return -1;
  if(nKaons!=2) return -1;
  //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
  //Fix implemented, loosening cut from 0.001 to 0.1. If >0.1, (-1*decay - 1) is returned.
  Int_t decayBstemp = decayBs;
  if(TMath::Abs(mcPart->Px()-sumPxDau)>0.1) decayBstemp = -1*decayBs - 1;
  if(TMath::Abs(mcPart->Py()-sumPyDau)>0.1) decayBstemp = -1*decayBs - 1;
  if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.1) decayBstemp = -1*decayBs - 1;
  decayBs = decayBstemp;
  return decayBs;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckLbDecay(AliMCEvent* mcEvent, Int_t label, Int_t* arrayDauLab){
  /// Checks the Lb decay channel. Returns >= 1 for Lb->Lcpi->pKpipi, <0 in other cases
  /// Returns 1 for non-resonant Lc decays and 2, 3 or 4 for resonant ones, -1 in other cases
  /// If rejected by momentum conservation check, return (-1*decay - 1) (to allow checks at task level)
  ///
  /// NB: Loosened cut on mom. conserv. (needed because of small issue in ITS Upgrade productions)

  if(label<0) return -1;
  AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(label);
  TParticle* part = mcEvent->Particle(label);
  if(!part || !mcPart) return -1;
  Int_t pdgD=part->GetPdgCode();
  if(TMath::Abs(pdgD)!=5122) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughterFirst();
  Int_t nKaons=0;
  Int_t nPions=0;
  Int_t nProtons=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundpKpi=0;
  Int_t decayLb = -1;

  for(Int_t iDau=0; iDau<nDau; iDau++){
    Int_t indDau = labelFirstDau+iDau;
    if(indDau<0) return -1;
    TParticle* dau=mcEvent->Particle(indDau);
    if(!dau) return -1;
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==4122){
      Int_t labDauLc[3] = {-1,-1,-1};
      //Returns 1 for non-resonant decays and 2, 3 or 4 for resonant ones, -1 in other cases
      Int_t decayLc = CheckLcpKpiDecay(mcEvent, indDau, labDauLc);

      //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
      //Fix implemented, to still select correct Lc decay
      if(decayLc==-2){
        AliMCParticle* mcLc = (AliMCParticle*)mcEvent->GetTrack(indDau);
        Int_t labelFirstDauLc = mcLc->GetDaughterFirst();
        if(mcLc->GetNDaughters() > 1){
          TParticle* dauLc1 = mcEvent->Particle(labelFirstDauLc);
          TParticle* dauLc2 = mcEvent->Particle(labelFirstDauLc+1);
          if(dauLc1 && dauLc2){
            Int_t pdgdauLc1=dauLc1->GetPdgCode();
            Int_t pdgdauLc2=dauLc2->GetPdgCode();
            if(TMath::Abs(pdgdauLc1)==313 || TMath::Abs(pdgdauLc2)==313) decayLc=2;
            else if(TMath::Abs(pdgdauLc1)==2224 || TMath::Abs(pdgdauLc2)==2224) decayLc=3;
            else if(TMath::Abs(pdgdauLc1)==3124 || TMath::Abs(pdgdauLc2)==3124) decayLc=4;
            else decayLc=1;
          }
        }
      }

      if (decayLc < 0 || labDauLc[0] == -1) return -1;
      decayLb = decayLc;
      nProtons++;
      nKaons++;
      nPions++;
      for(Int_t iLc = 0; iLc < 3; iLc++){
        TParticle* dauLc=mcEvent->Particle(labDauLc[iLc]);
        sumPxDau+=dauLc->Px();
        sumPyDau+=dauLc->Py();
        sumPzDau+=dauLc->Pz();
        arrayDauLab[nFoundpKpi++]=labDauLc[iLc];
      }
      if(nFoundpKpi>4) return -1;
    }else if(TMath::Abs(pdgdau)==211){
      //Temp fix for Lb->Lc+pi- and Lb->Lc-pi+ decays in ITS upgrade productions
      //if(pdgD*pdgdau>0) return -1;
      nPions++;
      sumPxDau+=dau->Px();
      sumPyDau+=dau->Py();
      sumPzDau+=dau->Pz();
      arrayDauLab[nFoundpKpi++]=indDau;
      if(nFoundpKpi>4) return -1;
    }
  }

  if(nProtons!=1) return -1;
  if(nKaons!=1) return -1;
  if(nPions!=2) return -1;
  //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
  //Fix implemented, loosening cut from 0.001 to 0.1. If >0.1, (-1*decay - 1) is returned.
  Int_t decayLbtemp = decayLb;
  if(TMath::Abs(part->Px()-sumPxDau)>0.1) decayLbtemp = -1*decayLb - 1;
  if(TMath::Abs(part->Py()-sumPyDau)>0.1) decayLbtemp = -1*decayLb - 1;
  if(TMath::Abs(part->Pz()-sumPzDau)>0.1) decayLbtemp = -1*decayLb - 1;
  decayLb = decayLbtemp;
  return decayLb;
}
//____________________________________________________________________________
Int_t AliVertexingHFUtils::CheckLbDecay(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab){
  /// Checks the Lb decay channel. Returns >= 1 for Lb->Lcpi->pKpipi, <0 in other cases
  /// Returns 1 for non-resonant Lc decays and 2, 3 or 4 for resonant ones, -1 in other cases
  /// If rejected by momentum conservation check, return (-1*decay - 1) (to allow checks at task level)
  ///
  /// NB: Loosened cut on mom. conserv. (needed because of small issue in ITS Upgrade productions)

  Int_t pdgD=mcPart->GetPdgCode();
  if(TMath::Abs(pdgD)!=5122) return -1;

  Int_t nDau=mcPart->GetNDaughters();
  if(nDau!=2) return -1;

  Int_t labelFirstDau = mcPart->GetDaughterLabel(0);
  Int_t nKaons=0;
  Int_t nPions=0;
  Int_t nProtons=0;
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;
  Int_t nFoundpKpi=0;
  Int_t decayLb = -1;

  for(Int_t iDau=0; iDau<nDau; iDau++){
    Int_t indDau = labelFirstDau+iDau;
    if(indDau<0) return -1;
    AliAODMCParticle* dau=dynamic_cast<AliAODMCParticle*>(arrayMC->At(indDau));
    if(!dau) return -1;
    Int_t pdgdau=dau->GetPdgCode();
    if(TMath::Abs(pdgdau)==4122){
      Int_t labDauLc[3] = {-1,-1,-1};
      //Returns 1 for non-resonant decays and 2, 3 or 4 for resonant ones, -1 in other cases
      Int_t decayLc = CheckLcpKpiDecay(arrayMC, dau, labDauLc);

      //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
      //Fix implemented, to still select correct Lc decay
      if(decayLc==-2){
        Int_t labelFirstDauLc = dau->GetDaughterLabel(0);
        if(dau->GetNDaughters() > 1){
          AliAODMCParticle* dauLc1=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labelFirstDauLc));
          AliAODMCParticle* dauLc2=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labelFirstDauLc+1));
          if(dauLc1 && dauLc2){
            Int_t pdgdauLc1=dauLc1->GetPdgCode();
            Int_t pdgdauLc2=dauLc2->GetPdgCode();
            if(TMath::Abs(pdgdauLc1)==313 || TMath::Abs(pdgdauLc2)==313) decayLc=2;
            else if(TMath::Abs(pdgdauLc1)==2224 || TMath::Abs(pdgdauLc2)==2224) decayLc=3;
            else if(TMath::Abs(pdgdauLc1)==3124 || TMath::Abs(pdgdauLc2)==3124) decayLc=4;
            else decayLc=1;
          }
        }
      }

      if (decayLc < 0 || labDauLc[0] == -1) return -1;
      decayLb = decayLc;
      nProtons++;
      nKaons++;
      nPions++;
      for(Int_t iLc = 0; iLc < 3; iLc++){
        AliAODMCParticle* dauLc=dynamic_cast<AliAODMCParticle*>(arrayMC->At(labDauLc[iLc]));
        sumPxDau+=dauLc->Px();
        sumPyDau+=dauLc->Py();
        sumPzDau+=dauLc->Pz();
        arrayDauLab[nFoundpKpi++]=labDauLc[iLc];
      }
      if(nFoundpKpi>4) return -1;
    }else if(TMath::Abs(pdgdau)==211){
      //Temp fix for Lb->Lc+pi- and Lb->Lc-pi+ decays in ITS upgrade productions
      //if(pdgD*pdgdau>0) return -1;
      nPions++;
      sumPxDau+=dau->Px();
      sumPyDau+=dau->Py();
      sumPzDau+=dau->Pz();
      arrayDauLab[nFoundpKpi++]=indDau;
      if(nFoundpKpi>4) return -1;
    }
  }

  if(nProtons!=1) return -1;
  if(nKaons!=1) return -1;
  if(nPions!=2) return -1;
  //Momentum conservation for several beauty decays not satisfied at gen. level in Upgrade MC's.
  //Fix implemented, loosening cut from 0.001 to 0.1. If >0.1, (-1*decay - 1) is returned.
  Int_t decayLbtemp = decayLb;
  if(TMath::Abs(mcPart->Px()-sumPxDau)>0.1) decayLbtemp = -1*decayLb - 1;
  if(TMath::Abs(mcPart->Py()-sumPyDau)>0.1) decayLbtemp = -1*decayLb - 1;
  if(TMath::Abs(mcPart->Pz()-sumPzDau)>0.1) decayLbtemp = -1*decayLb - 1;
  decayLb = decayLbtemp;
  return decayLb;
}
//________________________________________________________________________
Double_t AliVertexingHFUtils::GetSphericity(AliAODEvent* aod,
                                            Double_t etaMin, Double_t etaMax,
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
void AliVertexingHFUtils::GetSpherocity(AliAODEvent* aod,
                                        Double_t &spherocity, Double_t &phiRef,
                                        Double_t etaMin, Double_t etaMax,
                                        Double_t ptMin, Double_t ptMax,
                                        Int_t filtbit1, Int_t filtbit2,
                                        Int_t minMult, Double_t phiStepSizeDeg,
                                        Int_t nTrksToSkip, Int_t* idToSkip, 
                                        Bool_t unweight){
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

  if(nSelTracks<minMult){spherocity = -0.5; return;}

  //Getting thrust
  spherocity=2.;
  for(Int_t i=0; i<360/phiStepSizeDeg; ++i){
    Double_t phistep=TMath::Pi()*(Double_t)i*phiStepSizeDeg/180.;
    Double_t nx=TMath::Cos(phistep);
    Double_t ny=TMath::Sin(phistep);
    Double_t numer=0.;
    for(Int_t j=0; j<nSelTracks; ++j){
      Double_t pxA = 1.*TMath::Cos(phiArr[j]);  // x component of an unitary vector n for unweighted spherocity definition
      Double_t pyA = 1.*TMath::Sin(phiArr[j]);  // y component of an unitary vector n for unweighted spherocity definition 
        if(unweight == kFALSE){
            pxA=ptArr[j]*TMath::Cos(phiArr[j]);  // x component of an unitary vector n
            pyA=ptArr[j]*TMath::Sin(phiArr[j]);  // y component of an unitary vector n
        }
      numer+=TMath::Abs(ny*pxA - nx*pyA);
    }
    
    Double_t pFull = 0.;    
    if(unweight == kTRUE)
      pFull=numer*numer/(nSelTracks*nSelTracks);  
    else
      pFull=numer*numer/(sumpt*sumpt);

    if(pFull<spherocity){
        spherocity=pFull; // minimization;
        phiRef=phistep;
    }
  }

  delete [] ptArr;
  delete [] phiArr;

  spherocity*=(TMath::Pi()*TMath::Pi()/4.);
  return;
}

//________________________________________________________________________

void AliVertexingHFUtils::GetGeneratedSpherocity(TClonesArray *arrayMC,
                                                 Double_t &spherocity, Double_t &phiRef,
                                                 Double_t etaMin, Double_t etaMax,
                                                 Double_t ptMin, Double_t ptMax,
                                                 Int_t minMult, Double_t phiStepSizeDeg,
                                                 Bool_t unweight){

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

  if(nSelParticles<minMult){spherocity = -0.5; return;}

  //Getting thrust
  spherocity=2.;
  for(Int_t i=0; i<360/phiStepSizeDeg; ++i){
    Double_t phistep=TMath::Pi()*(Double_t)i*phiStepSizeDeg/180.;
    Double_t nx=TMath::Cos(phistep);
    Double_t ny=TMath::Sin(phistep);
    Double_t numer=0.;
    for(Int_t j=0; j<nSelParticles; ++j){
      Double_t pxA = 1.*TMath::Cos(phiArr[j]);  // x component of an unitary vector n for unweighted spherocity definition
      Double_t pyA = 1.*TMath::Sin(phiArr[j]);  // y component of an unitary vector n for unweighted spherocity definition 
        if(unweight == kFALSE){
            pxA=ptArr[j]*TMath::Cos(phiArr[j]);  // x component of an unitary vector n
            pyA=ptArr[j]*TMath::Sin(phiArr[j]);  // y component of an unitary vector n
        }
      numer+=TMath::Abs(ny*pxA - nx*pyA);
    }
    
    Double_t pFull = 0.;    
    if(unweight == kTRUE)
      pFull=numer*numer/(nSelParticles*nSelParticles);  
    else
      pFull=numer*numer/(sumpt*sumpt);

    if(pFull<spherocity){
        spherocity=pFull; // minimization;
        phiRef=phistep;
    }
  }

  delete [] ptArr;
  delete [] phiArr;

  spherocity*=(TMath::Pi()*TMath::Pi()/4.);
  return;
}

//________________________________________________________________________
TH1D* AliVertexingHFUtils::RebinHisto(TH1* hOrig, Int_t reb, Int_t firstUse){
  /// Rebin histogram, from bin firstUse to lastUse
  /// Use all bins if firstUse=-1
  /// If ngroup is not an exact divider of the number of bins,
  ///  the bin width is kept as reb*original width
  ///  and the range of rebinned histogram is adapted

  Int_t nBinOrig=hOrig->GetNbinsX();
  Int_t firstBinOrig=1;
  Int_t lastBinOrig=nBinOrig;
  Int_t nBinOrigUsed=nBinOrig;
  Int_t nBinFinal=nBinOrig/reb;
  if(firstUse>=1){
    firstBinOrig=firstUse;
    nBinFinal=(nBinOrig-firstUse+1)/reb;
    nBinOrigUsed=nBinFinal*reb;
    lastBinOrig=firstBinOrig+nBinOrigUsed-1;
  }else{
    Int_t exc=nBinOrigUsed%reb;
    if(exc!=0){
      nBinOrigUsed-=exc;
      lastBinOrig=firstBinOrig+nBinOrigUsed-1;
    }
  }

  printf("Rebin from %d bins to %d bins -- Used bins=%d in range %d-%d\n",nBinOrig,nBinFinal,nBinOrigUsed,firstBinOrig,lastBinOrig);
  Float_t lowLim=hOrig->GetXaxis()->GetBinLowEdge(firstBinOrig);
  Float_t hiLim=hOrig->GetXaxis()->GetBinUpEdge(lastBinOrig);
  TH1D* hRebin=new TH1D(Form("%s-rebin",hOrig->GetName()),hOrig->GetTitle(),nBinFinal,lowLim,hiLim);
  Int_t lastSummed=firstBinOrig-1;
  for(Int_t iBin=1;iBin<=nBinFinal; iBin++){
    Float_t sum=0.;
    Float_t sume2=0.;
    for(Int_t iOrigBin=0;iOrigBin<reb;iOrigBin++){
      sum+=hOrig->GetBinContent(lastSummed+1);
      sume2+=(hOrig->GetBinError(lastSummed+1)*hOrig->GetBinError(lastSummed+1));
      lastSummed++;
    }
    hRebin->SetBinContent(iBin,sum);
    hRebin->SetBinError(iBin,TMath::Sqrt(sume2));
  }
  return hRebin;
}
//________________________________________________________________________
TH1* AliVertexingHFUtils::AdaptTemplateRangeAndBinning(const TH1 *hMC,TH1 *hData, Double_t minFit, Double_t maxFit){
  /// Adapt the MC histograms (for signal and reflections) to the binning of the data histogram

  Int_t binmin=TMath::Max(1,hData->FindBin(hMC->GetXaxis()->GetXmin()));
  Bool_t found=kFALSE;
  Int_t binminD=-1;
  Int_t binminMC=-1;
  for(Int_t j=binmin; j<hData->GetNbinsX(); j++){
    if(found) break;
    for(Int_t k=1; k<hMC->GetNbinsX(); k++){
      Double_t delta=TMath::Abs(hMC->GetBinLowEdge(k)-hData->GetBinLowEdge(j));
      if(delta<0.0001){
	found=kTRUE;
	binminMC=k;
	binminD=j;
      }
      if(found) break;
    }
  }
  Int_t binmax=TMath::Min(hData->GetNbinsX(),hData->FindBin(hMC->GetXaxis()->GetXmax()*0.99999));
  found=kFALSE;
  Int_t binmaxD=-1;
  Int_t binmaxMC=-1;
  for(Int_t j=binmax; j>1; j--){
    if(found) break;
    for(Int_t k=hMC->GetNbinsX(); k>400; k--){
      Double_t delta=TMath::Abs(hMC->GetBinLowEdge(k+1)-hData->GetBinLowEdge(j+1));
      if(delta<0.0001){
	found=kTRUE;
	binmaxMC=k;
	binmaxD=j;
      }
      if(found) break;
    }
  }

  Double_t min=hData->GetBinLowEdge(binminD);
  Double_t max=hData->GetBinLowEdge(binmaxD)+hData->GetBinWidth(binmaxD);
  Double_t minMC=hMC->GetBinLowEdge(binminMC);
  Double_t maxMC=hMC->GetBinLowEdge(binmaxMC)+hMC->GetBinWidth(binmaxMC);
  Double_t width=hData->GetBinWidth(binminD);
  Double_t widthMC=hMC->GetBinWidth(binminMC);

  if(TMath::Abs(minMC-min)>0.0001*min || TMath::Abs(maxMC-max)>0.0001*max){
    printf("Cannot adapt range and rebin histo:\n");
    printf("Range for data histo: %f-%f GeV/c2    bins %d-%d width=%f\n",min,max,binminD,binmaxD,width);
    printf("Range for reflection histo: %f-%f GeV/c2    bins %d-%d width=%f\n",minMC,maxMC,binminMC,binmaxMC,widthMC);
    return 0x0;
  }

  Double_t rebin=width/widthMC;
  if(TMath::Abs(rebin-TMath::Nint(rebin))>0.001){
    printf("Cannot adapt histo: rebin %f issue, width MC = %f, width hData=%f (check=%f)\n",rebin,widthMC,width,TMath::Abs(rebin-TMath::Nint(rebin)));
    return 0x0;
  }

  Int_t nBinsNew=binmaxD-binminD+1;
  TH1 *hOut;
  TString stype=hMC->ClassName();
  if(stype.Contains("TH1F")){
    hOut=new TH1F(Form("%s-rebinned",hMC->GetName()),hMC->GetTitle(),nBinsNew,min,max);
  }else if(stype.Contains("TH1D")){
    hOut=new TH1D(Form("%s-rebinned",hMC->GetName()),hMC->GetTitle(),nBinsNew,min,max);
  }else{
    printf("Wrong type %s\n",stype.Data());
    return 0x0;
  }

  for(Int_t j=1; j<=hMC->GetNbinsX(); j++){
    Double_t m=hMC->GetBinCenter(j);
    Int_t binFin=hOut->FindBin(m);
    if(binFin>=1 && binFin<=nBinsNew){
      hOut->AddBinContent(binFin,hMC->GetBinContent(j));
    }
  }
  return hOut;
}

//___________________________________________________________________________________//
//method that performs simultaneus fit of in-plane and out-of-plane inv-mass spectra
ROOT::Fit::FitResult AliVertexingHFUtils::DoInPlaneOutOfPlaneSimultaneusFit(AliHFInvMassFitter *&massfitterInPlane, AliHFInvMassFitter *&massfitterOutOfPlane, 
                                                                            TH1F* hMassInPlane, TH1F* hMassOutOfPlane, Double_t MinMass, Double_t MaxMass, 
                                                                            Double_t massD, vector<UInt_t> commonpars) {

  cout << "\nIn-plane - out-of-plane simultaneus fit" << endl;
  cout << "\nIndependent prefits" << endl;

  //prefits to initialise parameters
  massfitterInPlane->SetUseLikelihoodFit();
  massfitterInPlane->SetInitialGaussianMean(massD);
  massfitterInPlane->SetInitialGaussianSigma(0.010);
  massfitterInPlane->MassFitter(kFALSE);

  massfitterOutOfPlane->SetUseLikelihoodFit();
  massfitterOutOfPlane->SetInitialGaussianMean(massD);
  massfitterOutOfPlane->SetInitialGaussianSigma(0.010);
  massfitterOutOfPlane->MassFitter(kFALSE);

  ROOT::Math::WrappedMultiTF1 wfInPlane(*(massfitterInPlane->GetMassFunc()),1);
  ROOT::Math::WrappedMultiTF1 wfOutOfPlane(*(massfitterOutOfPlane->GetMassFunc()),1);

  // set data options and ranges
  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange rangeMass; //same range for two functions

  rangeMass.SetRange(MinMass,MaxMass);
  ROOT::Fit::BinData dataInPlane(opt,rangeMass);
  ROOT::Fit::FillData(dataInPlane, hMassInPlane);
  ROOT::Fit::BinData dataOutOfPlane(opt,rangeMass);
  ROOT::Fit::FillData(dataOutOfPlane, hMassOutOfPlane);

  //define the 2 chi squares
  ROOT::Fit::Chi2Function chi2InPlane(dataInPlane, wfInPlane);
  ROOT::Fit::Chi2Function chi2OutOfPlane(dataOutOfPlane, wfOutOfPlane);

  //define the global chi square and get initial parameters from prefits
  const Int_t npars = massfitterInPlane->GetMassFunc()->GetNpar();
  const UInt_t ncommonpars = commonpars.size();
  Int_t nparsBkg = massfitterInPlane->GetBackgroundRecalcFunc()->GetNpar();
  Int_t nparsSgn = massfitterInPlane->GetSignalFunc()->GetNpar();
  Int_t nparsSecPeak = 0, nparsRefl = 0;
  if(massfitterInPlane->GetSecondPeakFunc())
    nparsSecPeak = massfitterInPlane->GetSecondPeakFunc()->GetNpar();
  if(massfitterInPlane->GetReflFunc())
    nparsRefl = massfitterInPlane->GetReflFunc()->GetNpar();

  GlobalInOutOfPlaneChi2 globalChi2(chi2InPlane, chi2OutOfPlane, npars, commonpars);

  vector<UInt_t>::iterator iter;
  vector<Double_t> initpars;
  for(Int_t iPar=0; iPar<2*npars; iPar++) {
    if(iPar<npars) { //in-plane
      initpars.push_back(massfitterInPlane->GetMassFunc()->GetParameter(iPar));
    }
    else { //out-of-plane
      iter = find(commonpars.begin(),commonpars.end(),iPar-npars);
      if(iter!=commonpars.end()) continue;
      else {
        initpars.push_back(massfitterOutOfPlane->GetMassFunc()->GetParameter(iPar-npars));
      }
    }
  }

  //define fitter and fit
  ROOT::Fit::Fitter simulfitter;
  simulfitter.Config().SetParamsSettings(npars*2-ncommonpars,initpars.data()); //set initial parameters from prefits
  if(nparsRefl>0) { //fix S/R
    simulfitter.Config().ParSettings(nparsBkg+nparsSgn+nparsSecPeak).Fix();

    iter = find(commonpars.begin(),commonpars.end(),npars+nparsBkg+nparsSgn+nparsSecPeak);
    if(iter==commonpars.end()) //if not included in common pars, need to be fixed also for out-of-plane func
      simulfitter.Config().ParSettings(npars+nparsBkg+nparsSgn+nparsSecPeak-ncommonpars).Fix();
  }

  simulfitter.Config().MinimizerOptions().SetPrintLevel(0);
  simulfitter.Config().SetMinimizer("Minuit2","Migrad");
  simulfitter.FitFCN(npars*2-ncommonpars,globalChi2,0,dataInPlane.Size()+dataOutOfPlane.Size(),kFALSE);
  ROOT::Fit::FitResult result = simulfitter.Result();
  cout << "\nSimultaneus fit" << endl;
  result.Print(cout);

  //Set new parameters to functions of mass fitters
  Int_t ncommonparsused = 0;
  for(Int_t iPar=0; iPar<npars; iPar++) {
    massfitterInPlane->GetMassFunc()->SetParameter(iPar,result.Parameter(iPar));
    massfitterInPlane->GetMassFunc()->SetParError(iPar,result.ParError(iPar));
    iter = find(commonpars.begin(),commonpars.end(),iPar);
    if(iter!=commonpars.end()) { //is common parameter
      massfitterOutOfPlane->GetMassFunc()->SetParameter(iPar,result.Parameter(iPar));
      massfitterOutOfPlane->GetMassFunc()->SetParError(iPar,result.ParError(iPar));
      ncommonparsused++;
    }
    else {
      massfitterOutOfPlane->GetMassFunc()->SetParameter(iPar,result.Parameter(iPar+npars-ncommonparsused));
      massfitterOutOfPlane->GetMassFunc()->SetParError(iPar,result.ParError(iPar+npars-ncommonparsused));
    }

    if(iPar < nparsBkg) { //bkg
      massfitterInPlane->GetBackgroundRecalcFunc()->SetParameter(iPar,massfitterInPlane->GetMassFunc()->GetParameter(iPar));
      massfitterInPlane->GetBackgroundRecalcFunc()->SetParError(iPar,massfitterInPlane->GetMassFunc()->GetParError(iPar));
      massfitterOutOfPlane->GetBackgroundRecalcFunc()->SetParameter(iPar,massfitterOutOfPlane->GetMassFunc()->GetParameter(iPar));
      massfitterOutOfPlane->GetBackgroundRecalcFunc()->SetParError(iPar,massfitterOutOfPlane->GetMassFunc()->GetParError(iPar));
    }
    else if(iPar >= nparsBkg && iPar < nparsBkg+nparsSgn){ //signal
      massfitterInPlane->GetSignalFunc()->SetParameter(iPar-nparsBkg,massfitterInPlane->GetMassFunc()->GetParameter(iPar));
      massfitterInPlane->GetSignalFunc()->SetParError(iPar-nparsBkg,massfitterInPlane->GetMassFunc()->GetParError(iPar));
      massfitterOutOfPlane->GetSignalFunc()->SetParameter(iPar-nparsBkg,massfitterOutOfPlane->GetMassFunc()->GetParameter(iPar));
      massfitterOutOfPlane->GetSignalFunc()->SetParError(iPar-nparsBkg,massfitterOutOfPlane->GetMassFunc()->GetParError(iPar));
    }
    else if(iPar >= nparsBkg && iPar < nparsBkg+nparsSgn){ //signal
      massfitterInPlane->GetSignalFunc()->SetParameter(iPar-nparsBkg,massfitterInPlane->GetMassFunc()->GetParameter(iPar));
      massfitterInPlane->GetSignalFunc()->SetParError(iPar-nparsBkg,massfitterInPlane->GetMassFunc()->GetParError(iPar));
      massfitterOutOfPlane->GetSignalFunc()->SetParameter(iPar-nparsBkg,massfitterOutOfPlane->GetMassFunc()->GetParameter(iPar));
      massfitterOutOfPlane->GetSignalFunc()->SetParError(iPar-nparsBkg,massfitterOutOfPlane->GetMassFunc()->GetParError(iPar));
    }
  }

  cout << "\n" << endl;
  return result;
}

//________________________________________________________________________
Double_t AliVertexingHFUtils::ComputeMaxd0MeasMinusExp(AliAODRecoDecayHF *cand, Double_t bfield) {
  // Compute maximum difference between observed and expected impact parameter of the candidate prongs
  Double_t dd0max = 0;
  UInt_t fNProngsCand = static_cast<UInt_t>(cand->GetNProngs());
  for (UInt_t iProng = 0; iProng < fNProngsCand; iProng++) {
    Double_t d0diff, errd0diff;
    cand->Getd0MeasMinusExpProng(iProng, bfield, d0diff, errd0diff);
    Double_t normdd0 = d0diff / errd0diff;
    if (iProng == 0 || TMath::Abs(normdd0) > TMath::Abs(dd0max))
        dd0max = normdd0;
  }
  return dd0max;
}

//______________________________________________________________________
Double_t AliVertexingHFUtils::CombineNsigmaTPCTOF(Double_t nsigmaTPC, Double_t nsigmaTOF) {
  // Combine TPC and TOF nsigma (sum in quadrature + special cases)
  if (nsigmaTPC > -998. && nsigmaTOF > -998.)
      return TMath::Sqrt((nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) / 2);
  else if (nsigmaTPC > -998. && nsigmaTOF < -998.)
      return TMath::Abs(nsigmaTPC);
  else if (nsigmaTPC < -998. && nsigmaTOF > -998.)
      return TMath::Abs(nsigmaTOF);
  else
      return -999.;
}

//______________________________________________________________________
Double_t AliVertexingHFUtils::CosPointingAngleFromKF(KFParticle kfp, KFParticle kfpmother)
{
  Double_t v[3];
  v[0] = kfp.GetX() - kfpmother.GetX();
  v[1] = kfp.GetY() - kfpmother.GetY();
  v[2] = kfp.GetZ() - kfpmother.GetZ();

  Double_t p[3];
  p[0] = kfp.GetPx();
  p[1] = kfp.GetPy();
  p[2] = kfp.GetPz();

  Double_t ptimesv2 = (p[0]*p[0]+p[1]*p[1]+p[2]*p[2])*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  if ( ptimesv2<=0 ) return 0.0;
  else {
    Double_t cos = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2]) / TMath::Sqrt(ptimesv2);
    if(cos >  1.0) cos =  1.0;
    if(cos < -1.0) cos = -1.0;
    return cos;
  }
}

//______________________________________________________________________
Double_t AliVertexingHFUtils::CosPointingAngleXYFromKF(KFParticle kfp, KFParticle kfpmother)
{
  Double_t v[2];
  v[0] = kfp.GetX() - kfpmother.GetX();
  v[1] = kfp.GetY() - kfpmother.GetY();

  Double_t p[2];
  p[0] = kfp.GetPx();
  p[1] = kfp.GetPy();

  Double_t ptimesv2 = (p[0]*p[0]+p[1]*p[1])*(v[0]*v[0]+v[1]*v[1]);

  if ( ptimesv2<=0 ) return 0.0;
  else {
    Double_t cos = (v[0]*p[0]+v[1]*p[1]) / TMath::Sqrt(ptimesv2);
    if(cos >  1.0) cos =  1.0;
    if(cos < -1.0) cos = -1.0;
    return cos;
  }
}

//______________________________________________________________________
Double_t AliVertexingHFUtils::CosThetaStarFromKF(Int_t ip, UInt_t pdgvtx, UInt_t pdgprong0, UInt_t pdgprong1, KFParticle kfpvtx, KFParticle kfpprong0, KFParticle kfpprong1)
{
  Double_t massvtx = TDatabasePDG::Instance()->GetParticle(pdgvtx)->Mass();
  Double_t massp[2];

  massp[0] = TDatabasePDG::Instance()->GetParticle(pdgprong0)->Mass();
  massp[1] = TDatabasePDG::Instance()->GetParticle(pdgprong1)->Mass();

  Double_t pStar = TMath::Sqrt((massvtx*massvtx-massp[0]*massp[0]-massp[1]*massp[1])*(massvtx*massvtx-massp[0]*massp[0]-massp[1]*massp[1])-4.*massp[0]*massp[0]*massp[1]*massp[1])/(2.*massvtx);

  Double_t e = kfpvtx.GetE();
  Double_t beta = kfpvtx.GetP()/e;
  Double_t gamma = e/massvtx;

  TVector3 mom;
  TVector3 momTot(kfpvtx.GetPx(), kfpvtx.GetPy(), kfpvtx.GetPz());

  if (ip==0) {
    mom.SetXYZ(kfpprong0.GetPx(), kfpprong0.GetPy(), kfpprong0.GetPz());
  }
  if (ip==1) {
    mom.SetXYZ(kfpprong1.GetPx(), kfpprong1.GetPy(), kfpprong1.GetPz());
  }

  Double_t cts = ( (mom.Dot(momTot)/momTot.Mag()) /gamma-beta*TMath::Sqrt(pStar*pStar+massp[ip]*massp[ip]) ) / pStar;

  return cts;
}

//______________________________________________________________________
Bool_t AliVertexingHFUtils::CheckAODvertexCov(AliAODVertex *vtx)                                  
{
  Double_t covMatrix[6];
  vtx->GetCovarianceMatrix(covMatrix);
  Double_t cov[3][3]={0.};
  cov[0][0] = covMatrix[0];
  cov[1][0] = covMatrix[1];
  cov[1][1] = covMatrix[2];
  cov[2][0] = covMatrix[3];
  cov[2][1] = covMatrix[4];
  cov[2][2] = covMatrix[5];
  if ( cov[0][0]<0 || cov[1][1]<0 || cov[2][2]<0 ) return kFALSE;
  for (Int_t i=0; i<3; i++) {
    for (Int_t j=0; j<3; j++) {
      if (i<=j) continue;
      if ( fabs(cov[i][j]) > TMath::Sqrt(cov[i][i]*cov[j][j]) ) return kFALSE;
    }
  }

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliVertexingHFUtils::CheckAODtrackCov(AliAODTrack *track)
{
  Double_t covMatrix[21];
  track->GetCovarianceXYZPxPyPz(covMatrix);
  Double_t cov[6][6]={0.};
  cov[0][0] = covMatrix[0];
  cov[1][0] = covMatrix[1];
  cov[1][1] = covMatrix[2];
  cov[2][0] = covMatrix[3];
  cov[2][1] = covMatrix[4];
  cov[2][2] = covMatrix[5];
  cov[3][0] = covMatrix[6];
  cov[3][1] = covMatrix[7];
  cov[3][2] = covMatrix[8];
  cov[3][3] = covMatrix[9];
  cov[4][0] = covMatrix[10];
  cov[4][1] = covMatrix[11];
  cov[4][2] = covMatrix[12];
  cov[4][3] = covMatrix[13];
  cov[4][4] = covMatrix[14];
  cov[5][0] = covMatrix[15];
  cov[5][1] = covMatrix[16];
  cov[5][2] = covMatrix[17];
  cov[5][3] = covMatrix[18];
  cov[5][4] = covMatrix[19];
  cov[5][5] = covMatrix[20];
  if ( cov[0][0]<0 || cov[1][1]<0 || cov[2][2]<0 || cov[3][3]<0 || cov[4][4]<0 || cov[5][5]<0 ) return kFALSE;
  for (Int_t i=0; i<6; i++) {
    for (Int_t j=0; j<6; j++) {
      if (i<=j) continue;
      if ( fabs(cov[i][j]) > TMath::Sqrt(cov[i][i]*cov[j][j]) ) return kFALSE;
    }
  }

  return kTRUE;
}

//______________________________________________________________________
Bool_t AliVertexingHFUtils::CheckKFParticleCov(KFParticle kfp)                                 
{
  if ( kfp.GetCovariance(0,0)<0 || kfp.GetCovariance(1,1)<0 || kfp.GetCovariance(2,2)<0 || kfp.GetCovariance(3,3)<0 || kfp.GetCovariance(4,4)<0 || kfp.GetCovariance(5,5)<0 ) return kFALSE;
  for (Int_t i=0; i<6; i++) {
    for (Int_t j=0; j<6; j++) {
      if (i<=j) continue;
      if ( fabs(kfp.GetCovariance(i,j)) > TMath::Sqrt(kfp.GetCovariance(i,i)*kfp.GetCovariance(j,j)) ) return kFALSE;
    }
  }

  return kTRUE;
}

//______________________________________________________________________
KFVertex AliVertexingHFUtils::CreateKFVertex(Double_t *param, Double_t *cov)                   
{
  KFPVertex kfpVtx;
  // Set the values
  Float_t paramF[3] = {(Float_t) param[0],(Float_t) param[1],(Float_t) param[2]};
  kfpVtx.SetXYZ(paramF);
  Float_t covF[6] = {(Float_t) cov[0],(Float_t) cov[1],(Float_t) cov[2],
                     (Float_t) cov[3],(Float_t) cov[4],(Float_t) cov[5]};
  kfpVtx.SetCovarianceMatrix(covF);

  return kfpVtx;
}

//______________________________________________________________________
KFVertex AliVertexingHFUtils::CreateKFVertexFromAODvertex(AliAODVertex *vtx)
{
  Double_t pos[3], cov[6];
  vtx->GetXYZ(pos);
  vtx->GetCovarianceMatrix(cov);
  Float_t covF[6];
  for (Int_t i=0; i<6; i++) { covF[i] = (Float_t)cov[i]; }

  KFPVertex kfpVertex;
  kfpVertex.SetXYZ((Float_t)pos[0], (Float_t)pos[1], (Float_t)pos[2]);
  kfpVertex.SetCovarianceMatrix(covF);
  kfpVertex.SetChi2(vtx->GetChi2());
  kfpVertex.SetNDF(vtx->GetNDF());
  kfpVertex.SetNContributors(vtx->GetNContributors());

  return kfpVertex;
}

//______________________________________________________________________
KFParticle AliVertexingHFUtils::CreateKFParticle(Double_t *param, Double_t *cov, Float_t Chi2perNDF, Int_t charge, Int_t pdg)
{
  // Interface to KFParticle
  KFPTrack kfpTrk;
  // Set the values
  kfpTrk.SetParameters((Float_t) param[0],(Float_t) param[1],(Float_t) param[2],
                       (Float_t) param[3],(Float_t) param[4],(Float_t) param[5]);
  kfpTrk.SetCharge(charge);
  Float_t covF[21];
  for (Int_t i = 0; i<21;i++) { covF[i] = (Float_t) cov[i]; }
  kfpTrk.SetCovarianceMatrix(covF);
  kfpTrk.SetNDF(1);
  kfpTrk.SetChi2(Chi2perNDF);

  // Build KFParticle
  KFParticle kfp(kfpTrk, pdg);
  return kfp;
}

//______________________________________________________________________
KFParticle AliVertexingHFUtils::CreateKFParticleFromAODtrack(AliAODTrack *track, Int_t pdg)
{
  Double_t trackParam[6];
  Double_t covMatrix[21];

  Bool_t IsDCA = track->GetXYZ(trackParam);
//  if (IsDCA) cout << "track position is at DCA" << endl;
//  if (!IsDCA) cout << "track position is at first point" << endl;
//  track->XvYvZv(trackParam);
  track->GetPxPyPz(&trackParam[3]);
  track->GetCovarianceXYZPxPyPz(covMatrix);

  KFPTrack kfpt;
  kfpt.SetParameters(trackParam);
  kfpt.SetCovarianceMatrix(covMatrix);
  kfpt.SetCharge(track->Charge());
  kfpt.SetNDF(1);
  kfpt.SetChi2(track->Chi2perNDF());

  KFParticle kfp(kfpt, pdg);

  return kfp;
}

//______________________________________________________________________
KFParticle AliVertexingHFUtils::CreateKFParticleV0(AliAODTrack *track1, AliAODTrack *track2, Int_t pdg1, Int_t pdg2)
{
  Double_t trackParam[6];
  Double_t covMatrix[21];

  track1->GetXYZ(trackParam);
  track1->GetPxPyPz(&trackParam[3]);
  track1->GetCovarianceXYZPxPyPz(covMatrix);

  KFPTrack kfpt1;
  kfpt1.SetParameters(trackParam);
  kfpt1.SetCovarianceMatrix(covMatrix);
  kfpt1.SetCharge(track1->Charge());
  kfpt1.SetNDF(1);
  kfpt1.SetChi2(track1->Chi2perNDF());

  track2->GetXYZ(trackParam);
  track2->GetPxPyPz(&trackParam[3]);
  track2->GetCovarianceXYZPxPyPz(covMatrix);

  KFPTrack kfpt2;
  kfpt2.SetParameters(trackParam);
  kfpt2.SetCovarianceMatrix(covMatrix);
  kfpt2.SetCharge(track2->Charge());
  kfpt2.SetNDF(1);
  kfpt2.SetChi2(track2->Chi2perNDF());

  // now we have all info to create the KFParticle version of the daughters
  KFParticle kfpDaughter1(kfpt1, pdg1);
  KFParticle kfpDaughter2(kfpt2, pdg2);
                       
  KFParticle kfpMother(kfpDaughter1, kfpDaughter2);

  return kfpMother;
}

//______________________________________________________________________
KFParticle AliVertexingHFUtils::CreateKFParticleCasc(KFParticle kfpV0, AliAODTrack *btrack, Int_t pdg_V0, Int_t pdg_btrack)
{
  Double_t trackParam[6];
  Double_t covMatrix[21];

  btrack->GetXYZ(trackParam);
//  btrack->XvYvZv(trackParam);
  btrack->GetPxPyPz(&trackParam[3]);
  btrack->GetCovarianceXYZPxPyPz(covMatrix);

  KFPTrack kfpbtrack;
  kfpbtrack.SetParameters(trackParam);
  kfpbtrack.SetCovarianceMatrix(covMatrix);
  kfpbtrack.SetCharge(btrack->Charge());
  kfpbtrack.SetNDF(1);
  kfpbtrack.SetChi2(btrack->Chi2perNDF());

  KFParticle kfpDaughter2(kfpbtrack, pdg_btrack);

  kfpV0.SetPDG(pdg_V0);
//  KFParticle kfpCasc(kfpV0, kfpDaughter2);
  KFParticle kfpCasc;
//  const KFParticle *vDaughters[2] = {&kfpV0, &kfpDaughter2};
  const KFParticle *vDaughters[2] = {&kfpDaughter2, &kfpV0}; // the order is important
  kfpCasc.Construct(vDaughters, 2);

  return kfpCasc;
}

//______________________________________________________________________
Double_t AliVertexingHFUtils::DecayLengthFromKF(KFParticle kfpParticle, KFParticle PV)
{
  Double_t dx_particle = PV.GetX()-kfpParticle.GetX();
  Double_t dy_particle = PV.GetY()-kfpParticle.GetY();
  Double_t dz_particle = PV.GetZ()-kfpParticle.GetZ();
  Double_t l_particle = TMath::Sqrt(dx_particle*dx_particle + dy_particle*dy_particle + dz_particle*dz_particle);
  return l_particle;
}

//______________________________________________________________________
Double_t AliVertexingHFUtils::DecayLengthXYFromKF(KFParticle kfpParticle, KFParticle PV)
{
  Double_t dx_particle = PV.GetX()-kfpParticle.GetX();
  Double_t dy_particle = PV.GetY()-kfpParticle.GetY();
  Double_t l_particle = TMath::Sqrt(dx_particle*dx_particle + dy_particle*dy_particle);
  return l_particle;
}

//______________________________________________________________________
Double_t AliVertexingHFUtils::ldlFromKF(KFParticle kfpParticle, KFParticle PV)
{
  Double_t dx_particle = PV.GetX()-kfpParticle.GetX();
  Double_t dy_particle = PV.GetY()-kfpParticle.GetY();
  Double_t dz_particle = PV.GetZ()-kfpParticle.GetZ();
  Double_t l_particle = TMath::Sqrt(dx_particle*dx_particle + dy_particle*dy_particle + dz_particle*dz_particle);
  Double_t dl_particle = (PV.GetCovariance(0)+kfpParticle.GetCovariance(0))*dx_particle*dx_particle + (PV.GetCovariance(2)+kfpParticle.GetCovariance(2))*dy_particle*dy_particle + (PV.GetCovariance(5)+kfpParticle.GetCovariance(5))*dz_particle*dz_particle + 2*( (PV.GetCovariance(1)+kfpParticle.GetCovariance(1))*dx_particle*dy_particle + (PV.GetCovariance(3)+kfpParticle.GetCovariance(3))*dx_particle*dz_particle + (PV.GetCovariance(4)+kfpParticle.GetCovariance(4))*dy_particle*dz_particle );
  if ( fabs(l_particle)<1.e-8f ) l_particle = 1.e-8f;
  dl_particle = dl_particle<0. ? 1.e8f : sqrt(dl_particle)/l_particle;
  if ( dl_particle==0. ) return 9999.;
  return l_particle/dl_particle;
}

//______________________________________________________________________
Double_t AliVertexingHFUtils::ldlXYFromKF(KFParticle kfpParticle, KFParticle PV)
{
  Double_t dx_particle = PV.GetX()-kfpParticle.GetX();
  Double_t dy_particle = PV.GetY()-kfpParticle.GetY();
  Double_t l_particle = TMath::Sqrt(dx_particle*dx_particle + dy_particle*dy_particle);
  Double_t dl_particle = (PV.GetCovariance(0)+kfpParticle.GetCovariance(0))*dx_particle*dx_particle + (PV.GetCovariance(2)+kfpParticle.GetCovariance(2))*dy_particle*dy_particle + 2*( (PV.GetCovariance(1)+kfpParticle.GetCovariance(1))*dx_particle*dy_particle );
  if ( fabs(l_particle)<1.e-8f ) l_particle = 1.e-8f;
  dl_particle = dl_particle<0. ? 1.e8f : sqrt(dl_particle)/l_particle;
  if ( dl_particle==0. ) return 9999.;
  return l_particle/dl_particle;
}
//______________________________________________________________________
TH1D* AliVertexingHFUtils::ComputeGenAccOverGenLimAcc(TFile* fileToyMCoutput,
						      Int_t nPtBins, Double_t* binLims,
						      Bool_t useSimpleFormula){
  // method to propagate the uncertainty on the ratio GenAcc/GenLimAcc
  // starting from ToyMC output file
  
  TH2D* hPtVsYGenAccToy=(TH2D*)fileToyMCoutput->Get("hPtVsYGenAcc");
  TH2D* hPtVsYGenLimAccToy=(TH2D*)fileToyMCoutput->Get("hPtVsYGenLimAcc");
  ComputeGenAccOverGenLimAcc(hPtVsYGenAccToy,hPtVsYGenLimAccToy,nPtBins,binLims,useSimpleFormula);
}
//______________________________________________________________________
TH1D* AliVertexingHFUtils::ComputeGenAccOverGenLimAcc(TH2D* hPtVsYGenAccToy, TH2D* hPtVsYGenLimAccToy,
						      Int_t nPtBins, Double_t* binLims,
						      Bool_t useSimpleFormula){

  // method to propagate the uncertainty on the ratio GenAcc/GenLimAcc
  // starting from the 2D histos (pt,y) of the ToyMC
  
  Int_t iybin1=hPtVsYGenAccToy->GetYaxis()->FindBin(-0.499999);
  Int_t iybin2=hPtVsYGenAccToy->GetYaxis()->FindBin(0.499999);
  TH1D* hptga=(TH1D*)hPtVsYGenAccToy->ProjectionX("hptga");
  TH1D* hptgay05=(TH1D*)hPtVsYGenAccToy->ProjectionX("hptga05",iybin1,iybin2);
  TH1D* hptgla=(TH1D*)hPtVsYGenLimAccToy->ProjectionX("hptgla");
  TH1D* hptgaR=0x0;
  TH1D* hptglaR=0x0;
  TH1D* hptgay05R=0x0;
  Bool_t shouldDelete=kFALSE;
  if(nPtBins>0 && binLims){
    hptgaR=(TH1D*)hptga->Rebin(nPtBins,"hptgaR",binLims);
    hptglaR=(TH1D*)hptgla->Rebin(nPtBins,"hptglaR",binLims);
    hptgay05R=(TH1D*)hptgay05->Rebin(nPtBins,"hptgay05R",binLims);
  }else{
    hptgaR=(TH1D*)hptga->Clone("hptgaR");
    hptglaR=(TH1D*)hptgla->Clone("hptglaR");
    hptgay05R=(TH1D*)hptgay05->Clone("hptgay05R");
    nPtBins=hptgay05R->GetNbinsX();
    binLims=new Double_t[nPtBins+1];
    for(Int_t ib=0; ib<nPtBins; ib++) binLims[ib]=hptgay05R->GetBinLowEdge(ib+1);
    binLims[nPtBins]=hptgay05R->GetBinLowEdge(nPtBins+1);
    shouldDelete=kTRUE;
  }
  TH1D* hGenAccOverGenLimAcc = new TH1D("hGenAccOverGenLimAcc",";p_{T} (GeV/c); GenAcc/GenLimAcc",nPtBins,binLims);
  for(Int_t ib=1; ib<=nPtBins; ib++){
    Double_t countGenAcc=hptgaR->GetBinContent(ib);
    Double_t countGenAccY05=hptgay05R->GetBinContent(ib);
    Double_t countGenLimAcc=hptglaR->GetBinContent(ib);
    Double_t yfid=0.8;
    Double_t ptcent=hptgaR->GetBinCenter(ib);
    if(ptcent<5) yfid=-0.2/15*ptcent*ptcent+1.9/15*ptcent+0.5;
    Double_t acc=0;
    Double_t erracc=0;
    if(countGenAcc>0 && countGenLimAcc>0){
      acc=countGenAcc/countGenLimAcc;
      Double_t rho=TMath::Sqrt(countGenAccY05*countGenAccY05/(countGenAcc*countGenLimAcc));
      if(useSimpleFormula) rho=TMath::Sqrt(0.5/yfid*acc/1.6);
      erracc=acc*TMath::Sqrt(1./countGenAcc+1./countGenLimAcc-2*rho*1/TMath::Sqrt(countGenAcc*countGenLimAcc));
    }
    hGenAccOverGenLimAcc->SetBinContent(ib,acc);
    hGenAccOverGenLimAcc->SetBinError(ib,erracc);
  }
  delete hptga;
  delete hptgay05;
  delete hptgla;
  delete hptgaR;
  delete hptgay05R;
  delete hptglaR;
  if(shouldDelete) delete binLims;
  return hGenAccOverGenLimAcc;
}
