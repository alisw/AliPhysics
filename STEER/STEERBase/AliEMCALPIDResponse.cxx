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
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliEMCALPIDResponse                                                  //
//                                                                      //
// EMCAL class to perfom PID                                            //
// This is a prototype and still under development                      //
//                                                                      //
// ---------------------------------------------------------------------//
// GetNumberOfSigmas():                                                 //
//                                                                      //
// Electrons:  Number of Sigmas for E/p value                           //
//             Parametrization of LHC11a (after recalibration)          //   
//                                                                      //
// NON electrons:                                                       //
//             Below or above E/p thresholds ( E/p < 0.5 || E/p > 1.5)  //
//             --> return +/- 999                                       //
//             Otherwise                                                //
//             --> return nsigma (parametrization of LHC10e)            //   
//                                                                      //
// ---------------------------------------------------------------------//
// ComputeEMCALProbability():                                           //
//                                                                      //
// Electrons:  Probability from Gaussian distribution                    //
//                                                                      //
// NON electrons:                                                       //
//             Below or above E/p thresholds ( E/p < 0.5 || E/p > 1.5)  //
//             --> probability to find particles below or above thr.    //
//             Otherwise                                                //
//             -->  Probability from Gaussian distribution              //  
//                  (proper normalization to each other?)               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TF1.h>
#include <TMath.h>

#include "AliEMCALPIDResponse.h"       //class header

#include "AliLog.h"   

ClassImp(AliEMCALPIDResponse)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliEMCALPIDResponse::AliEMCALPIDResponse():
  TObject(),
  fNorm(NULL)
{
  //
  //  The default constructor
  //

  for(Int_t i = 0; i < fNptBins; i++){

    fPtCutMin[i] = 0.0;

    for(Int_t j = 0; j < 2*AliPID::kSPECIES; j++){

      fMeanEoP[j][i]  = 0.0;
      fSigmaEoP[j][i] = 0.0;
      fProbLow[j][i]  = 0.0;
      fProbHigh[j][i] = 0.0;

    }
  } 
  fPtCutMin[fNptBins] = 0.0;

  fNorm = new TF1("fNorm","gaus",-20,20); 

  SetPtBoundary();
  SetParametrizations();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliEMCALPIDResponse::AliEMCALPIDResponse(const AliEMCALPIDResponse &other):
  TObject(other),
  fNorm(other.fNorm)
{
  //
  //  The copy constructor
  //
    for(Int_t i = 0; i < fNptBins; i++)
    {
	fPtCutMin[i] = 0.0;
	for(Int_t j = 0; j < 2*AliPID::kSPECIES; j++)
	{
	    fMeanEoP[j][i]  = 0.0;
	    fSigmaEoP[j][i] = 0.0;
	    fProbLow[j][i]  = 0.0;
	    fProbHigh[j][i] = 0.0;
	}
    }
  
    fPtCutMin[fNptBins] = 0.0;
    SetPtBoundary();
    SetParametrizations();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliEMCALPIDResponse & AliEMCALPIDResponse::operator=( const AliEMCALPIDResponse& other)
{
  //
  //  The assignment operator
  //

  if(this == &other) return *this;
  
  // Make copy
  TObject::operator=(other);
  fNorm = other.fNorm;

  for(Int_t i = 0; i < fNptBins; i++)
    {
	fPtCutMin[i] = 0.0;
	for(Int_t j = 0; j < 2*AliPID::kSPECIES; j++)
	{
	    fMeanEoP[j][i]  = 0.0;
	    fSigmaEoP[j][i] = 0.0;
	    fProbLow[j][i]  = 0.0;
	    fProbHigh[j][i] = 0.0;
	}
    }
  
  fPtCutMin[fNptBins] = 0.0;
  SetPtBoundary();
  SetParametrizations();

  return *this;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliEMCALPIDResponse::~AliEMCALPIDResponse() {

  delete fNorm;

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliEMCALPIDResponse::SetPtBoundary(){
  //
  // Set boundaries for momentum bins
  //
  fPtCutMin[0] = 1.5;
  fPtCutMin[1] = 2.5;
  fPtCutMin[2] = 3.5;
  fPtCutMin[3] = 4.5;
  fPtCutMin[4] = 5.5;
  fPtCutMin[5] = 6.5;

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliEMCALPIDResponse::SetParametrizations(){

// This are the preliminary parametrizations (hard coded)
// For electrons from LHC11a (Deepa Thomas)
// For NON-electrons from LHC10e (TOF/TPC analysis)
  
  // Gaussian mean
  Float_t mean[4][6]    = {
    { 0.932, 0.997, 0.998, 1.001, 1.011, 1.011 },                              // electrons
    { 0.227804, 0.34839, 0.404077, -0.107795, -4.14584, 0.5 },             // NON electrons
    { -2.10377, 0.0582898, 0.0582898, 0.0582898, 0.0582898, 0.0582898  },   // protons
    { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}                                              // anti-protons
  };        
  
  // Gaussian sigma
  Float_t sigma[4][6]= {
    { 0.0866, 0.0693, 0.0664, 0.0583, 0.0488, 0.0515  },                      // electrons
    { 0.310831, 0.267586, 0.404077, 0.381968, 1.46183, 0.314687 },          // NON electrons
    { 0.603209, 0.255332, 0.255332, 0.255332, 0.255332, 0.255332},          // protons
    { 0.516837, 0.351516,0.351516,0.351516,0.351516,0.351516 }              // anti - protons
  };

  // lower probability
  Float_t probL[3][6] = {
    { 0.928689, 0.938455, 0.940448, 0.948496, 0.955981, 0.951923 },         // NON electrons
    { 0.974518, 0.978088, 0.975089, 0.975089, 0.975089,0.975089},           // protons
    { 0.824037, 0.861149, 0.898734, 0.898734, 0.898734, 0.898734},          // anti - protons
  };

  // upper probability
  Float_t probH[3][6] = {
    { 0.00030227, 4.04106e-05, 0.000147406, 0., 0.000956938, 0.00106838 },  // NON electrons
    { 0.000157945, 0., 0., 0., 0., 0. },                                    // protons
    { 0.00343237, 0., 0., 0., 0., 0.}                                       // anti - protons
  };
  
  
  // set parametrizations
  Int_t spec = 0;
  for (Int_t species = 0; species < 2*AliPID::kSPECIES; species++) {          // first negative particles and then positive
    for (Int_t pt = 0; pt < fNptBins; pt++){

      switch(species){
      case 0:      // electrons
	spec = 0;
	break;
      case 4:      // anti - protons
	spec = 3;
	break;
      case 5:      // positrons
	spec = 0;
	break;
      case 9:      // protons
	spec = 2;
	break;
      default:     // NON electrons
	spec = 1;
	break;
      }
    

      fMeanEoP[species][pt]  = mean[spec][pt];    
      fSigmaEoP[species][pt] = sigma[spec][pt];  
      if( spec == 0) {     // electrons have NO lower and upper probability thresholds --> set to 0
	fProbLow[species][pt]  = 0.;
	fProbHigh[species][pt] = 0.;	
      }
      else{
	fProbLow[species][pt]  = probL[spec-1][pt];
	fProbHigh[species][pt] = probH[spec-1][pt];	
      } 
    
    }//loop pt bins
  }//loop species
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliEMCALPIDResponse::GetPtBin(Float_t pt) const {
  //
  // Returns the momentum bin index
  //

  Int_t i = -1;
  while(pt > fPtCutMin[i+1] && i+1 < fNptBins) i++;

  return i;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliEMCALPIDResponse::GetExpectedSignal( Float_t pt, AliPID::EParticleType n, Int_t charge) const  {
  //
  // Calculates the expected PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  //  

  Double_t signal = 0.;

  // Check the charge
  if( charge != -1 && charge != 1){
    
    return signal;
  }
  
  // Get the pt bin
  Int_t ptBin = GetPtBin(pt);
  
  // Get the species (first negative , then positive)
  Int_t species = n + AliPID::kSPECIES * ( charge + 1 ) / 2;

  // Get the signal
  if(species > -1 && species < 2*AliPID::kSPECIES && ptBin > -1 ){
    signal = fMeanEoP[species][ptBin];
  }

  return signal;

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliEMCALPIDResponse::GetExpectedSigma( Float_t pt, AliPID::EParticleType n,  Int_t charge) const  {
  //
  // Calculates the expected sigma of the PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  //  
  //
  
  Double_t sigma = 999.;

  // Check the charge
  if( charge != -1 && charge != 1){
    
    return sigma;
  }  

  // Get the pt bin
  Int_t ptBin = GetPtBin(pt);
  
  // Get the species (first negative , then positive)
  Int_t species = n + AliPID::kSPECIES * ( charge + 1 ) / 2;

  // Get the sigma
  if(species > -1 && species < 2*AliPID::kSPECIES && ptBin > -1 ){
    sigma = fSigmaEoP[species][ptBin];
  }

  return sigma;
 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliEMCALPIDResponse::GetExpectedNorm( Float_t pt, AliPID::EParticleType n,  Int_t charge) const  {
  //
  // Calculates the expected sigma of the PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  //  
  //
  
  Double_t norm = 1.;

  // Check the charge
  if( charge != -1 && charge != 1){
    
    return norm;
  }

  // Get the normalization factor ( Probability in the parametrized area / Integral of parametrized Gauss function in this area )
  fNorm->SetParameters(1./TMath::Sqrt(2*TMath::Pi()*GetExpectedSigma(pt,n,charge)*GetExpectedSigma(pt,n,charge)),GetExpectedSignal(pt,n,charge),GetExpectedSigma(pt,n,charge));
  norm = 1./fNorm->Integral(fLowEoP,fHighEoP)*(1-GetLowProb(pt,n,charge)-GetHighProb(pt,n,charge));

  return norm;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t  AliEMCALPIDResponse::GetNumberOfSigmas( Float_t pt,  Float_t eop, AliPID::EParticleType n,  Int_t charge) const {
      
  Double_t mean  = GetExpectedSignal(pt,n,charge);
  Double_t sigma = GetExpectedSigma(pt,n,charge);

  // if electron
  if(n == AliPID::kElectron){
    return (eop - mean) / sigma;
  }

  // if NON electron
  else{
    if ( eop < fLowEoP )
      return -999.;    // not parametrized 
    else if ( eop > fHighEoP )
      return 999.;     // not parametrized 
    else{
      return (eop - mean) / sigma; 
    }
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliEMCALPIDResponse::GetLowProb( Float_t pt, AliPID::EParticleType n,  Int_t charge) const {
  //
  //
  
  Double_t prob = 0.;

  // Check the charge
  if( charge != -1 && charge != 1){
    
    return prob;
  }
  
  // Get the pt bin
  Int_t ptBin = GetPtBin(pt);
  
  // Get the species (first negative , then positive)
  Int_t species = n + AliPID::kSPECIES * ( charge + 1 ) / 2;

  // Get the probability
  if(species > -1 && species < 2*AliPID::kSPECIES && ptBin > -1 ){
    prob = fProbLow[species][ptBin];
  }

  return prob;
 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliEMCALPIDResponse::GetHighProb( Float_t pt, AliPID::EParticleType n,  Int_t charge) const {
  //
  //
  
  Double_t prob = 0.;

  // Check the charge
  if( charge != -1 && charge != 1){
    
    return prob;
  }
  
  // Get the pt bin
  Int_t ptBin = GetPtBin(pt);
  
  // Get the species (first negative , then positive)
  Int_t species = n + AliPID::kSPECIES * ( charge + 1 ) / 2;

  // Get the probability
  if(species > -1 && species < 2*AliPID::kSPECIES && ptBin > -1 ){
    prob = fProbHigh[species][ptBin];
  }

  return prob;
  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliEMCALPIDResponse::ComputeEMCALProbability(Float_t pt, Float_t eop, Int_t charge, Double_t *pEMCAL) const {
  //
  //
  Double_t fRange  = 5.0;   // hardcoded 
  Double_t nsigma  = 0.0;
  Double_t sigma   = 0.0;
  Double_t proba   = 999.;
  

  // Check the charge
  if( charge != -1 && charge != 1){
    
    return proba;
  }

 
  // default value (will be returned, if pt below threshold)
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    pEMCAL[species] = 999.;
  }

  if( GetPtBin(pt) > -1 ){

    // set E/p range
    if(eop < 0.05) eop = 0.05;
    if(eop > 2.00) eop = 2.00;

    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {

      AliPID::EParticleType type = AliPID::EParticleType(species);

      // get nsigma value for each particle type at this E/p value
      nsigma = GetNumberOfSigmas(pt,eop,type,charge);
      sigma  = GetExpectedSigma(pt,type,charge);

      // electrons (standard Gaussian calculation of probabilities)
      if(type == AliPID::kElectron){
	if (TMath::Abs(nsigma) > fRange) {
	  pEMCAL[species]=TMath::Exp(-0.5*fRange*fRange)/TMath::Sqrt(2*TMath::Pi()*sigma*sigma);
	}
	else{
	  pEMCAL[species]=TMath::Exp(-0.5*(nsigma)*(nsigma))/TMath::Sqrt(2*TMath::Pi()*sigma*sigma);
	}
      }
      //NON electrons
      else{
	// E/p < 0.5  -->  return probability below E/p = 0.5
	if ( nsigma == -999){
	  pEMCAL[species] = GetLowProb(pt,type,charge);
	}
	// E/p > 1.5  -->  return probability above E/p = 1.5
	else if ( nsigma == 999){
	  pEMCAL[species] = GetHighProb(pt,type,charge);
	}
	// in parametrized region --> calculate probability for corresponding Gauss curve
	else{
	  pEMCAL[species]=TMath::Exp(-0.5*(nsigma)*(nsigma))/TMath::Sqrt(2*TMath::Pi()*sigma*sigma);
	
	  // normalize to total probability == 1
	  pEMCAL[species]*=GetExpectedNorm(pt,type,charge);
	}
      }
    }

    // return the electron probability
    proba = pEMCAL[AliPID::kElectron];  

  }

  return proba;

}
