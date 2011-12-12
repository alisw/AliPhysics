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

#include "TMath.h"
#include "TFile.h"

#include "AliOADBContainer.h"

// initialize static members
TH2F* AliPIDCombined::fDefaultPriorsTPC[5];

ClassImp(AliPIDCombined);

AliPIDCombined::AliPIDCombined():
	TNamed("CombinedPID","CombinedPID"),
	fDetectorMask(0),
	fRejectMismatchMask(0x7F),
	fEnablePriors(kTRUE),
	fSelectedSpecies(AliPID::kSPECIES),
	fUseDefaultTPCPriors(kFALSE)
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
	fSelectedSpecies(AliPID::kSPECIES),
	fUseDefaultTPCPriors(kFALSE)
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
    Int_t i = (Int_t)type;
    if (i >= AliPID::kSPECIES ) {
      if (i < (Int_t)AliPID::kDeuteron) {
	AliError(Form("Invalid EParticleType setting prior. Type: %d (neutral) not supported",i));
	return;
      } else i -= (AliPID::kSPECIESN-AliPID::kSPECIES);
    }
    if (i>(AliPID::kSPECIES+AliPID::kSPECIESLN)) {             // coverity fix (to put it mildly....)
      AliError(Form("Unexpected EParticleType setting prior. Type: %d (neutral) not supported",i));
      return;
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
	memset(p,0,sizeof(Double_t)*fSelectedSpecies);

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
	memset(priors,0,fSelectedSpecies*sizeof(Double_t));
	if (fEnablePriors){
	  GetPriors(track,priors,response->GetCurrentCentrality());
	  
	  // apply tof matching efficiency to priors if TOF joined PID for this track
	  if(fUseDefaultTPCPriors && (usedMask & AliPIDResponse::kDetTOF)){
	    Double_t pt=TMath::Abs(track->Pt());
	    Float_t kaonTOFfactor = 0.1;
	    if(pt > 0.35) kaonTOFfactor = 1 - TMath::Exp(-TMath::Power(pt,4.19618E-07)/5.68017E-01)*TMath::Power(pt,-1.50705);
	    Float_t protonTOFfactor = 0.1;
	    if(pt > 0.4) protonTOFfactor = 1 - TMath::Exp(-TMath::Power(pt,3.30978)/8.57396E-02)*TMath::Power(pt,-4.42661E-01);
	    
	    if(track->Charge() < 0){ // for negative tracks
	      kaonTOFfactor *= 1 - TMath::Exp(-TMath::Power(pt,4.87912E-07)/3.26431E-01)*TMath::Power(pt,-1.22893);
	      protonTOFfactor *= 1 - TMath::Exp(-TMath::Power(pt,2.00575E-07)/4.95605E-01)*TMath::Power(pt,-6.71305E-01);
	    }
	    
	    p[3] *= kaonTOFfactor;
	    p[4] *= protonTOFfactor;
	  }
	}
	else { for (Int_t i=0;i<fSelectedSpecies;i++) priors[i]=1.;}
	ComputeBayesProbabilities(bayesProbabilities,p,priors);
	return usedMask;
}


//-------------------------------------------------------------------------------------------------
void AliPIDCombined::GetPriors(const AliVTrack *track, Double_t* p,Float_t centrality) const {
	
	//
	// get priors from distributions
	//
	
	Double_t pt=TMath::Abs(track->Pt());

        if(fUseDefaultTPCPriors){ // use default priors if requested
	  Float_t usedCentr = centrality+5*(centrality>0);
	  if(usedCentr < -0.99) usedCentr = -0.99;
	  else if(usedCentr > 99) usedCentr = 99;
	  if(pt > 9.99) pt = 9.99;
	  else if(pt < 0.01)  pt = 0.01;
	  
	  for(Int_t i=0;i<5;i++) p[i] = fDefaultPriorsTPC[i]->Interpolate(usedCentr,pt);

	  return;
	}
	
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
void AliPIDCombined::GetPriors(const AliVTrack *track,Double_t* p,const AliPIDResponse *response,UInt_t detUsed) const{
	
	//
	// get priors from distributions
	//
        Double_t pt=TMath::Abs(track->Pt());
	
        if(fUseDefaultTPCPriors){ // use default priors if requested
	  Float_t usedCentr = response->GetCurrentCentrality()+5*(response->GetCurrentCentrality()>0);
	  if(usedCentr < -0.99) usedCentr = -0.99;
	  else if(usedCentr > 99) usedCentr = 99;
	  if(pt > 9.99) pt = 9.99;
	  else if(pt < 0.01)  pt = 0.01;
	  
	  for(Int_t i=0;i<5;i++) p[i] = fDefaultPriorsTPC[i]->Interpolate(usedCentr,pt);

	  // Extra factor if TOF matching was required
	  if(detUsed & AliPIDResponse::kDetTOF){
	    Float_t kaonTOFfactor = 0.1;
	    if(pt > 0.35) kaonTOFfactor = 1 - TMath::Exp(-TMath::Power(pt,4.19618E-07)/5.68017E-01)*TMath::Power(pt,-1.50705);
	    Float_t protonTOFfactor = 0.1;
	    if(pt > 0.4) protonTOFfactor = 1 - TMath::Exp(-TMath::Power(pt,3.30978)/8.57396E-02)*TMath::Power(pt,-4.42661E-01);
	    
	    if(track->Charge() < 0){ // for negative tracks
	      kaonTOFfactor *= 1 - TMath::Exp(-TMath::Power(pt,4.87912E-07)/3.26431E-01)*TMath::Power(pt,-1.22893);
	      protonTOFfactor *= 1 - TMath::Exp(-TMath::Power(pt,2.00575E-07)/4.95605E-01)*TMath::Power(pt,-6.71305E-01);
	    }
	    
	    p[3] *= kaonTOFfactor;
	    p[4] *= protonTOFfactor;
	  }

	  return;
	}


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

//----------------------------------------------------------------------------------------
void AliPIDCombined::SetDefaultTPCPriors(){
  fEnablePriors=kTRUE;
  fUseDefaultTPCPriors = kTRUE;

  TString oadbfilename("$ALICE_ROOT/OADB/COMMON/PID/data/");
  oadbfilename += "/PIDdefaultPriors.root";
  TFile * foadb = TFile::Open(oadbfilename.Data());
  if(!foadb->IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbfilename.Data()));
  
  AliOADBContainer * priorsContainer = (AliOADBContainer*) foadb->Get("priorsTPC");
  if (!priorsContainer) AliFatal("Cannot fetch OADB container for Priors");
  
  const char *namespecies[5] = {"El","Mu","Pi","Ka","Pr"};
  char name[100];

  for(Int_t i=0;i < 5;i++){
    snprintf(name,99,"hDefault%sPriors",namespecies[i]);
    fDefaultPriorsTPC[i] = (TH2F*) priorsContainer->GetDefaultObject(name);
    if (!fDefaultPriorsTPC[i]) AliFatal(Form("Cannot find priors for %s", namespecies[i]));
  }
}
