/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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


//-----------------------------------------------------------------
//        Base class for combining PID of different detectors    //
//        (user selected) and compute Bayesian probabilities     //
//                                                               //
//                                                               //
//   Origin: Pietro Antonioli, INFN-BO Pietro.Antonioli@cern.ch  //
//                                                               //
//-----------------------------------------------------------------

#include <TH1.h>

#include <AliVTrack.h>
#include <AliLog.h>
#include <AliPID.h>
#include <AliPIDResponse.h>

#include "AliPIDCombined.h"

AliPIDCombined::AliPIDCombined():
	TNamed("CombinedPID","CombinedPID"),
	fDetectorMask(0),
	fRejectMismatchMask(0x7F),
	fEnablePriors(kTRUE),
	fSelectedSpecies(AliPID::kSPECIES)
{	
  //
  // default constructor
  //	
  for (Int_t i=0;i<AliPID::kSPECIES+AliPID::kSPECIESLN;i++) fPriorsDistributions[i]=NULL;
  AliLog::SetClassDebugLevel("AliPIDCombined",10);
}

//-------------------------------------------------------------------------------------------------	
AliPIDCombined::AliPIDCombined(const TString& name,const TString& title):
	TNamed(name,title),
	fDetectorMask(0),
	fRejectMismatchMask(0x7F),
	fEnablePriors(kTRUE),
	fSelectedSpecies(AliPID::kSPECIES)
{
  //
  // default constructor with name and title
  //
  for (Int_t i=0;i<AliPID::kSPECIES+AliPID::kSPECIESLN;i++) fPriorsDistributions[i]=NULL;
  AliLog::SetClassDebugLevel("AliPIDCombined",10);

}

//-------------------------------------------------------------------------------------------------	
AliPIDCombined::~AliPIDCombined() {

  for(Int_t i=0;i < AliPID::kSPECIES+AliPID::kSPECIESLN;i++){
    if(fPriorsDistributions[i])
      delete fPriorsDistributions[i];
  }
}

//-------------------------------------------------------------------------------------------------	
void AliPIDCombined::SetPriorDistribution(AliPID::EParticleType type,TH1F *prior) {
  if ( (type < 0) || ( type >= (AliPID::kSPECIESN+AliPID::kSPECIESLN) ) ){
    AliError(Form("Invalid EParticleType setting prior (offending type: %d)",type));
    return;
  }
  if(prior) {
    Int_t i = type;
    if (type >= (AliPID::EParticleType)AliPID::kSPECIES ) {
      if (type < AliPID::kDeuteron) {
	AliError(Form("Invalid EParticleType setting prior. Type: %d (neutral) not supported",type));
	return;
      } else i = (Int_t)type - (AliPID::kSPECIESN-AliPID::kSPECIES);
    }
    if (fPriorsDistributions[i]) {
      delete fPriorsDistributions[i]; 
    }
    fPriorsDistributions[i]=new TH1F(*prior);
  }
}

//-------------------------------------------------------------------------------------------------	
UInt_t AliPIDCombined::ComputeProbabilities(const AliVTrack *track, const AliPIDResponse *response, Double_t* bayesProbabilities) const {
	Double_t pITS[fSelectedSpecies];
	Double_t pTPC[fSelectedSpecies];
	Double_t pTRD[fSelectedSpecies];
	Double_t pTOF[fSelectedSpecies];
	Double_t pHMPID[fSelectedSpecies];
	Double_t pEMCAL[fSelectedSpecies];
	Double_t pPHOS[fSelectedSpecies];
	for (Int_t i=0;i<fSelectedSpecies;i++) {
	 pITS[i]=1.;
	 pTPC[i]=1.;
	 pTRD[i]=1.;
	 pTOF[i]=1.;
	 pHMPID[i]=1.;
	 pEMCAL[i]=1.;
	 pPHOS[i]=1.;
	}
	UInt_t usedMask=0;
	AliPIDResponse::EDetPidStatus status=AliPIDResponse::kDetNoSignal;
	Double_t p[fSelectedSpecies];

	// getting probability distributions for selected detectors only
	if (fDetectorMask & AliPIDResponse::kDetITS) {
	  status = response->ComputeITSProbability(track, fSelectedSpecies, pITS);
	  SetCombinedStatus(status,&usedMask,AliPIDResponse::kDetITS,pITS);
	}

	if (fDetectorMask & AliPIDResponse::kDetTPC) { 
	  status = response->ComputeTPCProbability(track, fSelectedSpecies, pTPC);
	  SetCombinedStatus(status,&usedMask,AliPIDResponse::kDetTPC,pTPC);
	}


	if (fDetectorMask & AliPIDResponse::kDetTRD) { 
	  status = response->ComputeTRDProbability(track, fSelectedSpecies, pTRD);
	  SetCombinedStatus(status,&usedMask,AliPIDResponse::kDetTRD,pTRD);
	}

	if (fDetectorMask &  AliPIDResponse::kDetTOF) { 
	  status = response->ComputeTOFProbability(track, fSelectedSpecies, pTOF);
	  SetCombinedStatus(status,&usedMask,AliPIDResponse::kDetTOF,pTOF);
	}

	if (fDetectorMask & AliPIDResponse::kDetHMPID) { 
	  status = response->ComputeHMPIDProbability(track, fSelectedSpecies, pHMPID);
	  SetCombinedStatus(status,&usedMask,AliPIDResponse::kDetHMPID,pHMPID);
	}


	if (fDetectorMask & AliPIDResponse::kDetEMCAL) { 
	  status = response->ComputeEMCALProbability(track, fSelectedSpecies, pEMCAL);
	  SetCombinedStatus(status,&usedMask,AliPIDResponse::kDetEMCAL,pEMCAL);
	}


	if (fDetectorMask & AliPIDResponse::kDetPHOS) { 
	  status = response->ComputePHOSProbability(track, fSelectedSpecies, pPHOS);
	  SetCombinedStatus(status,&usedMask,AliPIDResponse::kDetPHOS,pPHOS);
	}


	for (Int_t i =0; i<fSelectedSpecies; i++){
	  p[i] = pITS[i]*pTPC[i]*pTRD[i]*pTOF[i]*pHMPID[i]*pEMCAL[i]*pPHOS[i];
	}
	Double_t priors[fSelectedSpecies];
	if (fEnablePriors) GetPriors(track,priors);
	else { for (Int_t i=0;i<fSelectedSpecies;i++) priors[i]=1.;}
	ComputeBayesProbabilities(bayesProbabilities,p,priors);
	return usedMask;
}


//-------------------------------------------------------------------------------------------------
void AliPIDCombined::GetPriors(const AliVTrack *track, Double_t* p) const {
	
	//
	// get priors from distributions
	//
	
	Double_t pt=TMath::Abs(track->Pt());
	Double_t sumPriors = 0;
	for (Int_t i=0;i<fSelectedSpecies;++i){
		if (fPriorsDistributions[i]){
			p[i]=fPriorsDistributions[i]->Interpolate(pt);
		}
		else {
			p[i]=0.;
		}
		sumPriors+=p[i];		
	}

	// normalizing priors
	if (sumPriors == 0) return;
	for (Int_t i=0;i<fSelectedSpecies;++i){
	   p[i]=p[i]/sumPriors;
	}
	return;
}

//-------------------------------------------------------------------------------------------------	
void AliPIDCombined::ComputeBayesProbabilities(Double_t* probabilities, const Double_t* probDensity, const Double_t* prior) const {


  //
  // calculate Bayesian probabilities
  //
  Double_t sum = 0.;
  for (Int_t i = 0; i < fSelectedSpecies; i++) {
    sum += probDensity[i] * prior[i];
  }
  if (sum <= 0) {
    AliError("Invalid probability densities or priors");
    for (Int_t i = 0; i < fSelectedSpecies; i++) probabilities[i] = -1;
    return;
  }
  for (Int_t i = 0; i < fSelectedSpecies; i++) {
    probabilities[i] = probDensity[i] * prior[i] / sum;
  }


}

//----------------------------------------------------------------------------------------
void AliPIDCombined::SetCombinedStatus(AliPIDResponse::EDetPidStatus status, UInt_t *maskDetIn, AliPIDResponse::EDetCode bit, Double_t* p) const {
  switch (status) {
  case AliPIDResponse::kDetNoSignal:
    break;
  case AliPIDResponse::kDetPidOk:
    *maskDetIn = *maskDetIn | bit;
    break;
  case AliPIDResponse::kDetMismatch:
    if ( fRejectMismatchMask & bit) for (Int_t j=0;j<fSelectedSpecies;j++) p[j]=1./fSelectedSpecies;
    else *maskDetIn = *maskDetIn | bit;
    break;
  }

}

ClassImp(AliPIDCombined);

