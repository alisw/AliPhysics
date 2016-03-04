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
TH2F* AliPIDCombined::fDefaultPriorsTPC[]={0x0};
TH2F* AliPIDCombined::fDefaultPriorsTPCpPb[]={0x0};
Float_t AliPIDCombined::fTOFmismProb = 0;

ClassImp(AliPIDCombined);

AliPIDCombined::AliPIDCombined():
	TNamed("CombinedPID","CombinedPID"),
	fDetectorMask(0),
	fRejectMismatchMask(0x7F),
	fEnablePriors(kFALSE),
	fSelectedSpecies(AliPID::kSPECIES),
	fUseDefaultTPCPriors(kFALSE)
{	
  //
  // default constructor
  //	
  for (Int_t i=0;i<AliPID::kSPECIESC;i++) fPriorsDistributions[i]=NULL;
  AliLog::SetClassDebugLevel("AliPIDCombined",10);
}

//-------------------------------------------------------------------------------------------------	
AliPIDCombined::AliPIDCombined(const TString& name,const TString& title):
	TNamed(name,title),
	fDetectorMask(0),
	fRejectMismatchMask(0x7F),
	fEnablePriors(kFALSE),
	fSelectedSpecies(AliPID::kSPECIES),
	fUseDefaultTPCPriors(kFALSE)
{
  //
  // default constructor with name and title
  //
  for (Int_t i=0;i<AliPID::kSPECIESC;i++) fPriorsDistributions[i]=NULL;
  AliLog::SetClassDebugLevel("AliPIDCombined",10);

}

//-------------------------------------------------------------------------------------------------	
AliPIDCombined::~AliPIDCombined() {

  for(Int_t i=0;i < AliPID::kSPECIESC;i++){
    if(fPriorsDistributions[i])
      delete fPriorsDistributions[i];
  }
}

//-------------------------------------------------------------------------------------------------	
void AliPIDCombined::SetPriorDistribution(AliPID::EParticleType type,TH1F *prior) {
  if ( (type < 0) || ( type >= ((AliPID::EParticleType)AliPID::kSPECIESC) ) ){
    AliError(Form("Invalid EParticleType setting prior (offending type: %d)",type));
    return;
  }
  if(prior) {
    Int_t i = (Int_t)type;
    if (fPriorsDistributions[i]) {
      delete fPriorsDistributions[i]; 
    }
    fPriorsDistributions[i]=new TH1F(*prior);
    SetEnablePriors(kTRUE); 
  }
}

//-------------------------------------------------------------------------------------------------	
UInt_t AliPIDCombined::ComputeProbabilities(const AliVTrack *track, const AliPIDResponse *response, Double_t* bayesProbabilities,Double_t *priorsOwn) const {
  //
  // (1) Get raw probabilities of requested detectors and combine
  // (2) priors passed as argument
  // (3) Compute Bayes probabilities
  //


  // (1) Get raw probabilities of selected detectors and combine
  UInt_t usedMask=0;             // detectors actually used for track
  fTOFmismProb = 0; // reset TOF mismatch weights

  AliPIDResponse::EDetPidStatus status=AliPIDResponse::kDetNoSignal;
  Double_t p[AliPID::kSPECIESC];  // combined probabilities of selected detectors
  Double_t pMismTOF[AliPID::kSPECIESC];  // combined TOF mismatch probabilities using selected detectors
  for (Int_t i=0;i<fSelectedSpecies;i++){ p[i]=1.;pMismTOF[i]=1.;} // no decision
  for (Int_t ibit = 0; ibit < 7 ; ibit++) {
    AliPIDResponse::EDetCode detBit = (AliPIDResponse::EDetCode)(1<<ibit);
    if (fDetectorMask & detBit) {  	    // getting probabilities for requested detectors only
      Double_t detProb[AliPID::kSPECIESC];
      status = response->ComputePIDProbability(detBit,track,fSelectedSpecies,detProb);
      if (status == AliPIDResponse::kDetPidOk) {
	if (fRejectMismatchMask & detBit) {         // mismatch check (currently just for TOF)
	  if (detBit == AliPIDResponse::kDetTOF) {
	    fTOFmismProb = response->GetTOFMismatchProbability(); // mismatch weights computed with TOF probs (no arguments)
	    //Float_t probMis = response->GetTOFMismatchProbability(track); // mismatch compatibility TPC-TOF cut
	    SetCombinedStatus(status,&usedMask,detBit,detProb,fTOFmismProb);
	  } else {
	    SetCombinedStatus(status,&usedMask,detBit);
	  }
	} else {
	  SetCombinedStatus(status,&usedMask,detBit);
	}
	for (Int_t i=0;i<fSelectedSpecies;i++){
	  p[i] *= detProb[i];
	  if(detBit == AliPIDResponse::kDetTOF){
	    Double_t pt = track->Pt();
	    Double_t mismPropagationFactor[] = {1.,1.,1.,1. + TMath::Exp(1. - 1.12*pt),
						1. + 1./(4.71114 - 5.72372*pt + 2.94715*pt*pt),1.,1.,1.,1.,1.}; // correction for kaons and protons which has to be alligned with the one in AliPIDResponse
	    pMismTOF[i] *= fTOFmismProb*mismPropagationFactor[i];
	  }
	  else pMismTOF[i] *= detProb[i];
	}
      }
    }
  }
  // if no detectors available there is no point to go further
  if (usedMask == 0) return usedMask;

  // (2) Get priors and propagate depending on detectors used
  Double_t priors[AliPID::kSPECIESC];
  memset(priors,0,fSelectedSpecies*sizeof(Double_t));

  if(priorsOwn){
    for(Int_t ipr=0;ipr < fSelectedSpecies;ipr++)
      priors[ipr] = priorsOwn[ipr];
  }

  else if (fEnablePriors){
    Bool_t isPPB = (response->GetBeamType() == AliPIDResponse::kPPB);
    GetPriors(track,priors,response->GetCurrentCentrality(),isPPB);

    // We apply the propagation factors of the more external detector
    // 
    // TPC+HMPID    --> apply HMPID propagation factors (TRD and TOF may be present)
    // TPC+EMCAL    --> apply EMCAL propagation factors (TRD and TOF may be present)
    // TPC+TOF      --> apply TOF propagation factors (TRD may be present, HMPID and EMCAL not (if requested))
    // TPC+TRD      --> apply TRD propagation factors (TOF, HMPID and EMCAL not present (if requested) )
    // 
    if(fUseDefaultTPCPriors) {

      Double_t pt=TMath::Abs(track->Pt());
      if ( ( (usedMask & AliPIDResponse::kDetEMCAL)==AliPIDResponse::kDetEMCAL) || ( (usedMask & AliPIDResponse::kDetHMPID)==AliPIDResponse::kDetHMPID) ) {
      // we assume EMCAL and HMPID cannot be simultaneously present
      if ( (usedMask & AliPIDResponse::kDetEMCAL)==AliPIDResponse::kDetEMCAL ) {
	// EMCal case (for the moment only in combination with TPC)
	// propagation factors determined from LHC11d MC (LHC12a15f)
	// v2 clusterizer, dEta < 0.015, dPhi < 0.03, NonLinearityFunction = 6
	       
	Float_t electronEMCALfactor=0.1;
	Float_t kaonEMCALfactor = 0.1;
	Float_t protonEMCALfactor = 0.1;
	       
	if(track->Charge() > 0){
	  // positiv charge (start parametrization at 0.75 GeV/c and stop at 20 GeV/c for the moment)
	  if(pt > 0.75 && pt < 20.0){
	    electronEMCALfactor =  ( 0.214159 * ( 1 - TMath::Exp(-TMath::Power(pt,0.484512)/0.700499)*TMath::Power(pt,-0.669644)) ) /  ( 0.210436 * ( 1 - TMath::Exp(-TMath::Power(pt,-0.219228)/0.947432)*TMath::Power(pt,-0.700792)) );
	    kaonEMCALfactor =  ( 0.208686 * ( 1 - TMath::Exp(-TMath::Power(pt,-3.98149e-05)/1.28447)*TMath::Power(pt,-0.629191)) ) /  ( 0.210436 * ( 1 - TMath::Exp(-TMath::Power(pt,-0.219228)/0.947432)*TMath::Power(pt,-0.700792)) );
	    protonEMCALfactor =  ( 0.27555 * ( 1 - TMath::Exp(-TMath::Power(pt,-1.97226e-05)/1.52719)*TMath::Power(pt,-0.209068)) ) /  ( 0.210436 * ( 1 - TMath::Exp(-TMath::Power(pt,-0.219228)/0.947432)*TMath::Power(pt,-0.700792)) );
		   
	  }
	}
	       
	if(track->Charge() < 0){
	  // negative charge  (start parametrization at 0.75 GeV/c and stop at 20 GeV/c for the moment)
	  if(pt > 0.75 && pt < 20.0){ 		   
	    electronEMCALfactor =  ( 0.216895 * ( 1 - TMath::Exp(-TMath::Power(pt,0.000105924)/0.865938)*TMath::Power(pt,-1.32787)) ) /  ( 0.210385 * ( 1 - TMath::Exp(-TMath::Power(pt,4.41206e-07)/1.08984)*TMath::Power(pt,-0.544375)) );
	    kaonEMCALfactor =  ( 0.204117 * ( 1 - TMath::Exp(-TMath::Power(pt,-1.6853e-05)/1.61765)*TMath::Power(pt,-0.738355)) ) /  ( 0.210385 * ( 1 - TMath::Exp(-TMath::Power(pt,4.41206e-07)/1.08984)*TMath::Power(pt,-0.544375)) );
	    protonEMCALfactor =  ( 0.215679 * ( 1 - TMath::Exp(-TMath::Power(pt,-4.10015e-05)/1.40921)*TMath::Power(pt,-0.533752)) ) /  ( 0.210385 * ( 1 - TMath::Exp(-TMath::Power(pt,4.41206e-07)/1.08984)*TMath::Power(pt,-0.544375)) );
	  }
	}
	priors[0] *= electronEMCALfactor;
	priors[3] *= kaonEMCALfactor;
	priors[4] *= protonEMCALfactor;	      
      } // end of EMCAL case

      else if ( (usedMask & AliPIDResponse::kDetHMPID)==AliPIDResponse::kDetHMPID ) {  // HMPID case

	Float_t kaonHMPIDfactor   = 0.;
	Float_t protonHMPIDfactor = 0.;                     
        if(pt>1. && pt<6.)  kaonHMPIDfactor   = (-0.0729337 + pt*0.0999531 - pt*pt*0.0371803 + pt*pt*pt*0.00706436 - pt*pt*pt*pt*0.000643619 + pt*pt*pt*pt*pt*2.21853e-05)/(-0.00896231+ pt*0.0330702 - pt*pt*0.0109562+ pt*pt*pt*0.00232895 - pt*pt*pt*pt*0.000246143 + pt*pt*pt*pt*pt*9.59812e-06); 
        if(pt>1.4 && pt<6.) protonHMPIDfactor = (-0.0444188 + pt*0.0681506 - pt*pt*0.0231819 + pt*pt*pt*0.00400771 - pt*pt*pt*pt*0.000339315 + pt*pt*pt*pt*pt*1.12616e-05)/(-0.00896231+ pt*0.0330702 - pt*pt*0.0109562+ pt*pt*pt*0.00232895 - pt*pt*pt*pt*0.000246143 + pt*pt*pt*pt*pt*9.59812e-06);        
        if(pt>6. && pt<8.)  kaonHMPIDfactor   = (-0.0729337 + pt*0.0999531 - pt*pt*0.0371803 + pt*pt*pt*0.00706436 - pt*pt*pt*pt*0.000643619 + pt*pt*pt*pt*pt*2.21853e-05)/0.0530456;        
        if(pt>8.)           kaonHMPIDfactor   = 0.0550432/0.0530456; 
        if(pt>6. && pt<8.5) protonHMPIDfactor = (-0.0444188 + pt*0.0681506 - pt*pt*0.0231819 + pt*pt*pt*0.00400771 - pt*pt*pt*pt*0.000339315 + pt*pt*pt*pt*pt*1.12616e-05)/0.0530456;      
        if(pt>8.5)          protonHMPIDfactor = 0.0530071/0.0530456;       
                                 
        if(track->Charge() < 0){ 
         if(pt>0.4 && pt<6.) protonHMPIDfactor = (-0.0351485 + pt*0.0473821 - pt*pt*0.0147947 + pt*pt*pt*0.00254811- pt*pt*pt*pt*0.000224724 + pt*pt*pt*pt*pt*7.9303e-06)/(-0.00896231+ pt*0.0330702 - pt*pt*0.0109562+ pt*pt*pt*0.00232895 - pt*pt*pt*pt*0.000246143 + pt*pt*pt*pt*pt*9.59812e-06);
         if(pt>6. && pt<8.5) protonHMPIDfactor = (-0.0351485 + pt*0.0473821 - pt*pt*0.0147947 + pt*pt*pt*0.00254811- pt*pt*pt*pt*0.000224724 + pt*pt*pt*pt*pt*7.9303e-06)/0.0530456; 
         if(pt>8.5)          protonHMPIDfactor = 0.0457756/0.0530456; 
	} 
      
        priors[3] *= kaonHMPIDfactor;
        priors[4] *= protonHMPIDfactor;
               
      }

    } // end of outer cases: EMCAL/HMPID
      else if ( (usedMask & AliPIDResponse::kDetTOF) == AliPIDResponse::kDetTOF ){
	Float_t kaonTOFfactor = 0.1;
	if(pt > 0.29){
	  kaonTOFfactor = 1 - TMath::Exp(-TMath::Power(pt,4.19618E-07)/4.96502e-01)*TMath::Power(pt,-1.50705);
	}

	Float_t protonTOFfactor = 0.1;
	if(pt > 0.4) protonTOFfactor = 1 - TMath::Exp(-TMath::Power(pt,3.30978)/8.57396E-02)*TMath::Power(pt,-4.42661E-01);

	if(track->Charge() < 0){ // for negative tracks
	  kaonTOFfactor *= 1 - TMath::Exp(-TMath::Power(pt,4.87912E-07)/3.26431E-01)*TMath::Power(pt,-1.22893);
	  protonTOFfactor *= 1 - TMath::Exp(-TMath::Power(pt,2.00575E-07)/4.95605E-01)*TMath::Power(pt,-6.71305E-01);
	}
	// TODO: we may need an electron factor for TOF as well, especially if TRD is present!
	priors[3] *= kaonTOFfactor;
	priors[4] *= protonTOFfactor;
      }  // end of TOF case
      else if ( (usedMask & AliPIDResponse::kDetTRD)==AliPIDResponse::kDetTRD ) {
	Float_t electronTRDfactor=0.1;
	Float_t kaonTRDfactor = 0.1;
	Float_t protonTRDfactor = 0.1;
		 
	if(track->Charge() > 0){
	  // positiv charge
	  if(pt > 0.25) electronTRDfactor =  1 - TMath::Exp(-TMath::Power(pt,5.13315e-03)/2.11145e-01)*TMath::Power(pt,-2.97659e+00);
	  if(pt > 0.35) kaonTRDfactor = 1 - TMath::Exp(-TMath::Power(pt,-4.29549e-02)/4.87989e-01)*TMath::Power(pt,-1.54270e+00);
	  if(pt > 0.35) protonTRDfactor = 1 - TMath::Exp(-TMath::Power(pt,2.81238e+00)/7.57082e-02)*TMath::Power(pt,-8.12595e-01);
	}
		 
	if(track->Charge() < 0){
	  // negative charge
	  if(pt > 0.25) electronTRDfactor =  1 - TMath::Exp(-TMath::Power(pt,2.45537e-02)/1.90397e-01)*TMath::Power(pt,-3.33121e+00);
	  if(pt > 0.35) kaonTRDfactor = 1 - TMath::Exp(-TMath::Power(pt, -3.42831e-03)/5.57013e-01)*TMath::Power(pt,-1.39202e+00);
	  if(pt > 0.35) protonTRDfactor = 1 - TMath::Exp(-TMath::Power(pt,3.36631e+00)/7.18819e-02)*TMath::Power(pt,-2.00577e-01);
	}
	// what about electrons
	priors[0] *= electronTRDfactor;
	priors[3] *= kaonTRDfactor;
	priors[4] *= protonTRDfactor;	      
      } // end of TRD case
    } // end of fUseDefaultTPCPriors
  }   // end of use priors
  else { for (Int_t i=0;i<fSelectedSpecies;i++) priors[i]=1.;}

  //  Compute Bayes probabilities
  ComputeBayesProbabilities(bayesProbabilities,p,priors,pMismTOF);

  // compute TOF probability contribution from mismatch
  fTOFmismProb = 0; 
  for (Int_t i=0;i<fSelectedSpecies;i++) fTOFmismProb += pMismTOF[i];

  return usedMask;
}


//-------------------------------------------------------------------------------------------------
void AliPIDCombined::GetPriors(const AliVTrack *track, Double_t* p,Float_t centrality,Bool_t isPPB) const {
	
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

	  if(! isPPB)
	    for(Int_t i=0;i<AliPID::kSPECIESC;i++) p[i] = fDefaultPriorsTPC[i]->Interpolate(usedCentr,pt);
	  else
	    for(Int_t i=0;i<AliPID::kSPECIESC;i++) p[i] = fDefaultPriorsTPCpPb[i]->Interpolate(usedCentr,pt);

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
	  
	  for(Int_t i=0;i<AliPID::kSPECIESC;i++) p[i] = fDefaultPriorsTPC[i]->Interpolate(usedCentr,pt);

	  // Extra factor if TOF matching was required
	  if(detUsed & AliPIDResponse::kDetTOF){
	    Float_t kaonTOFfactor = 0.1;
	    if(pt > 0.29){
	      kaonTOFfactor = 1 - TMath::Exp(-TMath::Power(pt,4.19618E-07)/4.96502e-01)*TMath::Power(pt,-1.50705);
	    }
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
void AliPIDCombined::ComputeBayesProbabilities(Double_t* probabilities, const Double_t* probDensity, const Double_t* prior, Double_t* probDensityMism) const {


  //
  // calculate Bayesian probabilities
  //
  Double_t sum = 0.;
  for (Int_t i = 0; i < fSelectedSpecies; i++) {
    sum += probDensity[i] * prior[i];
  }
  if (sum <= 0) {

    AliError("Invalid probability densities or priors");
    for (Int_t i = 0; i < fSelectedSpecies; i++) probabilities[i] = 1./fSelectedSpecies;
    return;
  }
  for (Int_t i = 0; i < fSelectedSpecies; i++) {
    probabilities[i] = probDensity[i] * prior[i] / sum;
    if(probDensityMism) probDensityMism[i] *= prior[i] / sum;
  }


}


//----------------------------------------------------------------------------------------
void AliPIDCombined::SetCombinedStatus(AliPIDResponse::EDetPidStatus status, UInt_t *maskDetIn, AliPIDResponse::EDetCode bit) const {
  switch (status) {
  case AliPIDResponse::kDetNoParams:
  case AliPIDResponse::kDetNoSignal:
  case AliPIDResponse::kDetMismatch: // for backward compatibilty, we need then to remove kDetMismatch from AliPIDResponse
    break;
  case AliPIDResponse::kDetPidOk:
    *maskDetIn = *maskDetIn | bit;
    break;
  }

}

//----------------------------------------------------------------------------------------
void AliPIDCombined::SetCombinedStatus(AliPIDResponse::EDetPidStatus status, UInt_t *maskDetIn, AliPIDResponse::EDetCode bit, Double_t* /*p*/, const Float_t /*probMis*/) const {
  switch (status) {
  case AliPIDResponse::kDetNoParams:
  case AliPIDResponse::kDetNoSignal:
  case AliPIDResponse::kDetMismatch: // for backward compatibility, we need then to remove kDeteMismatch from AliPIDResponse
    break;
  case AliPIDResponse::kDetPidOk:
    //      if (probMis > 0.5) for (Int_t j=0;j<fSelectedSpecies;j++) p[j]=1./fSelectedSpecies; // no longer used because mismatch is in the framework now
    //else 
    *maskDetIn = *maskDetIn | bit;
    break;
  }

}



//----------------------------------------------------------------------------------------
void AliPIDCombined::SetDefaultTPCPriors(){
  fEnablePriors=kTRUE;
  fUseDefaultTPCPriors = kTRUE;

  //Check if priors are already initialized
  if (fDefaultPriorsTPC[0]) return;
  
  TString oadbfilename("$ALICE_PHYSICS/OADB/COMMON/PID/data/");
  //  TString oadbfilename("$ALICE_ROOT/OADB/COMMON/PID/data/");
  //  TString oadbfilename(Form("%s/COMMON/PID/data/",AliAnalysisManager::GetOADBPath()));
  oadbfilename += "/PIDdefaultPriors.root";
  TFile * foadb = TFile::Open(oadbfilename.Data());
  if(!foadb || !foadb->IsOpen()) {
    delete foadb;
    AliFatal(Form("Cannot open OADB file %s", oadbfilename.Data()));
    return;
  }
  
  AliOADBContainer * priorsContainer = (AliOADBContainer*) foadb->Get("priorsTPC");
  if (!priorsContainer) AliFatal("Cannot fetch OADB container for Priors");
  
  const char *namespecies[AliPID::kSPECIESC] = {"El","Mu","Pi","Ka","Pr","De","Tr","He3","He3"};
  char name[100];

  for(Int_t i=0;i < AliPID::kSPECIESC;i++){
    snprintf(name,99,"hDefault%sPriors",namespecies[i]);
    fDefaultPriorsTPC[i] = (TH2F*) priorsContainer->GetDefaultObject(name);
    if (!fDefaultPriorsTPC[i]) AliFatal(Form("Cannot find priors for %s", namespecies[i]));
    snprintf(name,99,"hDefault%sPriorsPPb",namespecies[i]);
    fDefaultPriorsTPCpPb[i] = (TH2F*) priorsContainer->GetDefaultObject(name);
    if (!fDefaultPriorsTPCpPb[i]) {
        fDefaultPriorsTPCpPb[i] = fDefaultPriorsTPC[i]; // use PbPb ones
    }
  }

  delete foadb;
}
