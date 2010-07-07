/**************************************************************************
 * Copyright(c) 2005-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
// ITS PID method # 1
//           Implementation of the ITS PID class
// Very naive one... Should be made better by the detector experts...
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------
#include "TMath.h"
#include "AliITSPIDResponse.h"
#include "AliITSPidParams.h"
#include "AliExternalTrackParam.h"

ClassImp(AliITSPIDResponse)

AliITSPIDResponse::AliITSPIDResponse(Bool_t isMC): 
  fRes(0.13),
  fKp1(15.77),
  fKp2(4.95),
  fKp3(0.312),
  fKp4(2.14),
  fKp5(0.82)
{
  if(!isMC){
    fBBtpcits[0]=0.73;
    fBBtpcits[1]=14.68;
    fBBtpcits[2]=0.905;
    fBBtpcits[3]=1.2;
    fBBtpcits[4]=6.6;
    fBBsa[0]=5.33458E4;
    fBBsa[1]=16.5303;
    fBBsa[2]=2.60065E-3;
    fBBsa[3]=3.59533E-4;
    fBBsa[4]=7.51168E-5;  
    for(Int_t i=0; i<5;i++) fResolSA[i]=0.15;
    for(Int_t i=0; i<5;i++) fResolTPCITS[i]=0.13;
  }else{
    fBBtpcits[0]=0.73;
    fBBtpcits[1]=14.68;
    fBBtpcits[2]=0.905;
    fBBtpcits[3]=0.2;
    fBBtpcits[4]=6.6;
    fBBsa[0]=139.1;
    fBBsa[1]=23.36;
    fBBsa[2]=0.06052;
    fBBsa[3]=0.2043;
    fBBsa[4]=-0.0004999;
    for(Int_t i=0; i<5;i++) fResolSA[i]=0.15;
    for(Int_t i=0; i<5;i++) fResolTPCITS[i]=0.13;
  }
}

//_________________________________________________________________________
AliITSPIDResponse::AliITSPIDResponse(Double_t *param): 
  fRes(param[0]),
  fKp1(15.77),
  fKp2(4.95),
  fKp3(0.312),
  fKp4(2.14),
  fKp5(0.82)
{
  //
  //  The main constructor
  //
}


Double_t AliITSPIDResponse::BetheAleph(Double_t p, Double_t mass) const {
  //
  // returns AliExternalTrackParam::BetheBloch normalized to 
  // fgMIP at the minimum
  //
  
  Double_t bb=
    AliExternalTrackParam::BetheBlochAleph(p/mass,fKp1,fKp2,fKp3,fKp4,fKp5);
  return bb;
}

Double_t AliITSPIDResponse::Bethe(Double_t p, Double_t mass, Bool_t isSA) const {
  //
  // returns AliExternalTrackParam::BetheBloch normalized to 
  // fgMIP at the minimum
  //

  Double_t bg=p/mass;
  Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);
  Double_t gamma=bg/beta;
  Double_t par[5];
  if(isSA){
    for(Int_t ip=0; ip<5;ip++) par[ip]=fBBsa[ip];
  }else{
    for(Int_t ip=0; ip<5;ip++) par[ip]=fBBtpcits[ip];
  }
  Double_t eff=1.0;
  if(bg<par[2])
    eff=(bg-par[3])*(bg-par[3])+par[4];
  else
    eff=(par[2]-par[3])*(par[2]-par[3])+par[4];

  Double_t bb=0.;
  if(gamma>=0. && beta>0.){
    bb=(par[1]+2.0*TMath::Log(gamma)-beta*beta)*(par[0]/(beta*beta))*eff;
  }
  return bb;
}

Double_t AliITSPIDResponse::GetResolution(Double_t bethe, 
					  Int_t nPtsForPid, 
					  Bool_t isSA) const {
  // 
  // Calculate expected resolution for truncated mean
  //
  Float_t r;
  if(isSA) r=fResolSA[nPtsForPid];
  else r=fResolTPCITS[nPtsForPid];
  return r*bethe;
}




void AliITSPIDResponse::GetITSProbabilities(Float_t mom, Double_t qclu[4], Double_t condprobfun[AliPID::kSPECIES]) const {
  //
  // Method to calculate PID probabilities for a single track
  // using the likelihood method
  //
  const Int_t nLay = 4;
  const Int_t nPart = 3;

  static AliITSPidParams pars;  // Pid parametrisation parameters
  
  Double_t itsProb[nPart] = {1,1,1}; // p, K, pi

  for (Int_t iLay = 0; iLay < nLay; iLay++) {
    if (qclu[iLay] <= 0)
      continue;

    Float_t dedx = qclu[iLay];
    Float_t layProb = pars.GetLandauGausNorm(dedx,AliPID::kProton,mom,iLay+3);
    itsProb[0] *= layProb;
    
    layProb = pars.GetLandauGausNorm(dedx,AliPID::kKaon,mom,iLay+3);
    if (mom < 0.16) layProb=0.00001;
    itsProb[1] *= layProb;
    
    layProb = pars.GetLandauGausNorm(dedx,AliPID::kPion,mom,iLay+3);
    itsProb[2] *= layProb;
  }

  // Normalise probabilities
  Double_t sumProb = 0;
  for (Int_t iPart = 0; iPart < nPart; iPart++) {
    sumProb += itsProb[iPart];
  }

  for (Int_t iPart = 0; iPart < nPart; iPart++) {
    itsProb[iPart]/=sumProb;
  }
  
  condprobfun[AliPID::kElectron] = itsProb[2]/3.;
  condprobfun[AliPID::kMuon] = itsProb[2]/3.;
  condprobfun[AliPID::kPion] = itsProb[2]/3.;
  condprobfun[AliPID::kKaon] = itsProb[1];
  condprobfun[AliPID::kProton] = itsProb[0];
  return;
}

