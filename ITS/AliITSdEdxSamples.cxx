/**************************************************************************
 * Copyright(c) 2009-2012, ALICE Experiment at CERN, All rights reserved. *
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
/* $Id$ */



///////////////////////////////////////////////////////////////////
//                                                               //
// Class to store information for PID with ITS                   //
// and truncated mean computation methods                        //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSPidParams.h"
#include "AliITSdEdxSamples.h"
#include "AliLog.h"
#include <TMath.h>

ClassImp(AliITSdEdxSamples)

//______________________________________________________________________
AliITSdEdxSamples::AliITSdEdxSamples():TObject(),
  fNSamples(0),
  fP(0.),
  fParticleSpecie(0)
{
  // Default constructor

  for(Int_t i=0; i<fgkMaxSamples; i++) {
    fdEdxSamples[i]=0.;
    fPAtSample[i]=0.;
  }
}

//______________________________________________________________________
AliITSdEdxSamples::AliITSdEdxSamples(Int_t nSamples,Double_t* samples, Double_t mom, Int_t specie) :
  TObject(),
  fNSamples(nSamples),
  fP(mom),
  fParticleSpecie(specie)
{
  // Standard constructor
  for(Int_t i=0; i<fgkMaxSamples; i++)fPAtSample[i]=0.;
  SetSamples(nSamples,samples);
}
//______________________________________________________________________
void AliITSdEdxSamples::SetSamples(Int_t nSamples, Double_t* samples){
  // Set the samples

  if(nSamples>fgkMaxSamples){
    AliWarning(Form("Too many dE/dx samples,only first %d will be used",fgkMaxSamples));
    fNSamples=fgkMaxSamples;
  }else{
    fNSamples=nSamples;
  }
  for(Int_t i=0; i<fNSamples; i++) fdEdxSamples[i]=samples[i];
  for(Int_t i=fNSamples; i<fgkMaxSamples; i++) fdEdxSamples[i]=0.;
  return;
}
//______________________________________________________________________
void AliITSdEdxSamples::SetSamplesAndMomenta(Int_t nSamples, Double_t* samples, Double_t* mom){
  // Set the samples
  SetSamples(nSamples,samples);
  for(Int_t i=0; i<fNSamples; i++) fPAtSample[i]=mom[i];
  for(Int_t i=fNSamples; i<fgkMaxSamples; i++) fPAtSample[i]=0.;
  return;
}
//______________________________________________________________________
Double_t AliITSdEdxSamples::GetTruncatedMean(Double_t frac, Double_t mindedx) const {
  // compute truncated mean 

  Int_t nc=0;
  Double_t dedx[fgkMaxSamples];
  for (Int_t il=0; il<fNSamples; il++) { // count good (>0) dE/dx values
    if(fdEdxSamples[il]>mindedx){
      dedx[nc]= fdEdxSamples[il];
      nc++;
    }
  }
  if(nc<1) return 0.;
  
  Int_t swap; // sort in ascending order
  do {
    swap=0;
    for (Int_t i=0; i<nc-1; i++) {
      if (dedx[i]<=dedx[i+1]) continue;
      Double_t tmp=dedx[i];
      dedx[i]=dedx[i+1]; 
      dedx[i+1]=tmp;
      swap++;
    }
  } while (swap);

  Double_t sumamp=0,sumweight=0;
  Double_t weight[fgkMaxSamples];
  for(Int_t iw=0; iw<fgkMaxSamples; iw++) weight[iw]=0.;
  Int_t lastUsed=(Int_t)(frac*nc+0.00001);
  if(lastUsed==0) lastUsed=1;
  if(lastUsed>fgkMaxSamples) lastUsed=fgkMaxSamples;
  for(Int_t iw=0; iw<lastUsed; iw++) weight[iw]=1.;
  if((frac*nc-lastUsed)>0.4 && lastUsed<fgkMaxSamples) weight[lastUsed]=0.5;
  for (Int_t i=0; i<nc; i++) {
    AliDebug(5,Form("dE/dx %f   weight %f",dedx[i],weight[i]));
    sumamp+= dedx[i]*weight[i];
    sumweight+=weight[i];
  }
  if(sumweight>0.) return sumamp/sumweight;
  else return 0.;

}
//______________________________________________________________________
Double_t AliITSdEdxSamples::GetWeightedMean(Double_t mindedx) const {
  // compute generalized mean with k=-2 (used by CMS)

  Int_t nc=0;
  Double_t dedx[fgkMaxSamples];
  for (Int_t il=0; il<fNSamples; il++) { // count good (>0) dE/dx values
    if(fdEdxSamples[il]>mindedx){
      dedx[nc]= fdEdxSamples[il];
      nc++;
    }
  }
  if(nc<1) return 0.;

  Double_t weiSum = 0.;
  for (Int_t i=0; i<nc; i++) {
    weiSum+=TMath::Power(dedx[i],-2);
  }
  Double_t wMean=0.;
  if(weiSum>0.) wMean= TMath::Power(weiSum/nc,-0.5);
  return wMean;

}
//______________________________________________________________________
void  AliITSdEdxSamples::GetConditionalProbabilities(AliITSPidParams* pars, Double_t condprob[AliPID::kSPECIES], Double_t mindedx) const {
  // compute conditional probablilities
  const Int_t nPart = 3;
  Double_t itsProb[nPart] = {1,1,1}; // p, K, pi

  for(Int_t iS=0; iS<fNSamples; iS++){
    Int_t iLayer=iS+3; // to match with present ITS
    if(iLayer>6) iLayer=6; // all extra points are treated as SSD
    if(fdEdxSamples[iS]<mindedx) continue;
    Float_t dedx = fdEdxSamples[iS];
    Float_t layProb = pars->GetLandauGausNorm(dedx,AliPID::kProton,fP,iLayer);
    itsProb[0] *= layProb;

    layProb = pars->GetLandauGausNorm(dedx,AliPID::kKaon,fP,iLayer);
    if (fP < 0.16) layProb=0.00001;
    itsProb[1] *= layProb;
    
    layProb = pars->GetLandauGausNorm(dedx,AliPID::kPion,fP,iLayer);
    itsProb[2] *= layProb;
  }

  // Normalise probabilities
  Double_t sumProb = 0;
  for (Int_t iPart = 0; iPart < nPart; iPart++) {
    sumProb += itsProb[iPart];
  }
  sumProb += 2*itsProb[2]; // muon and electron cannot be distinguished from pions

  for (Int_t iPart = 0; iPart < nPart; iPart++) {
    itsProb[iPart]/=sumProb;
  }
  
  condprob[AliPID::kElectron] = itsProb[2];
  condprob[AliPID::kMuon] = itsProb[2];
  condprob[AliPID::kPion] = itsProb[2];
  condprob[AliPID::kKaon] = itsProb[1];
  condprob[AliPID::kProton] = itsProb[0];
  return;

}
