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
//             --> return +/- 99                                        //
//             Otherwise                                                //
//             --> return nsigma (parametrization of LHC10e)            //
//                                                                      //
// NO Parametrization (outside pT range): --> return -999               //
//                                                                      //
// ---------------------------------------------------------------------//
// ComputeEMCALProbability():                                           //
//                                                                      //
// Electrons:  Probability from Gaussian distribution                   //
//                                                                      //
// NON electrons:                                                       //
//             Below or above E/p thresholds ( E/p < 0.5 || E/p > 1.5)  //
//             --> probability to find particles below or above thr.    //
//             Otherwise                                                //
//             -->  Probability from Gaussian distribution              //
//                  (proper normalization to each other?)               //
//                                                                      //
// NO Parametrization (outside pT range): --> return 999                //
//////////////////////////////////////////////////////////////////////////

#include <TF1.h>
#include <TMath.h>

#include "AliEMCALPIDResponse.h"       //class header

#include "AliLog.h"   

ClassImp(AliEMCALPIDResponse)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliEMCALPIDResponse::AliEMCALPIDResponse():
  TObject(),
  fNorm(NULL),
  fkPIDParams(NULL)
{
  //
  //  The default constructor
  //


  fNorm = new TF1("fNorm","gaus",-20,20); 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliEMCALPIDResponse::AliEMCALPIDResponse(const AliEMCALPIDResponse &other):
  TObject(other),
  fNorm(other.fNorm),
  fkPIDParams(other.fkPIDParams)
{
  //
  //  The copy constructor
  //

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
  fkPIDParams = other.fkPIDParams;

 
  return *this;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliEMCALPIDResponse::~AliEMCALPIDResponse() {

  delete fNorm;

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

  // Get the parameters for this particle type and pt
  const TVectorD *params = GetParams(n, pt);

  // IF not in momentum range, NULL is returned --> return default value
  if(!params) return norm;

  Double_t mean     = (*params)[2];   // mean value of Gausiian parametrization
  Double_t sigma    = (*params)[3];   // sigma value of Gausiian parametrization
  Double_t eopMin   = (*params)[4];   // min E/p value for parametrization
  Double_t eopMax   = (*params)[5];   // max E/p value for parametrization
  Double_t probLow  = (*params)[6];   // probability to be below eopMin
  Double_t probHigh = (*params)[7];   // probability to be above eopMax

  // Get the normalization factor ( Probability in the parametrized area / Integral of parametrized Gauss function in this area )
  fNorm->SetParameters(1./TMath::Sqrt(2*TMath::Pi()*sigma*sigma),mean,sigma);
  norm = 1./fNorm->Integral(eopMin,eopMax)*(1-probLow-probHigh);

  return norm;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t  AliEMCALPIDResponse::GetNumberOfSigmas( Float_t pt,  Float_t eop, AliPID::EParticleType n,  Int_t charge) const {
      
  Double_t nsigma = -999.;

  // Check the charge
  if( charge != -1 && charge != 1){
    return nsigma;
  }

  // Get the parameters for this particle type and pt
  const TVectorD *params = GetParams(n, pt);

  // IF not in momentum range, NULL is returned --> return default value
  if(!params) return nsigma;

  Double_t mean     = (*params)[2];   // mean value of Gausiian parametrization
  Double_t sigma    = (*params)[3];   // sigma value of Gausiian parametrization
  Double_t eopMin   = (*params)[4];   // min E/p value for parametrization
  Double_t eopMax   = (*params)[5];   // max E/p value for parametrization

  // if electron
  if(n == AliPID::kElectron){
    if(sigma != 0) nsigma = (eop - mean) / sigma;
  }

  // if NON electron
  else{
    if ( eop < eopMin )
      nsigma = -99;    // not parametrized (smaller than eopMin)
    else if ( eop > eopMax )
      nsigma = 99.;     // not parametrized (bigger than eopMax)
    else{
      if(sigma != 0) nsigma = (eop - mean) / sigma; 
    }
  }

  return nsigma;

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliEMCALPIDResponse::ComputeEMCALProbability(Float_t pt, Float_t eop, Int_t charge, Double_t *pEMCAL) const {
  //
  //
  Double_t fRange  = 5.0;   // hardcoded (???)
  Double_t nsigma  = 0.0;
  Double_t proba   = 999.;
  

  // Check the charge
  if( charge != -1 && charge != 1){
    return proba;
  }

 
  // default value (will be returned, if pt below threshold)
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    pEMCAL[species] = 999.;
  }

  // set E/p range
  if(eop < 0.05) eop = 0.05;
  if(eop > 2.00) eop = 2.00;
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    
    AliPID::EParticleType type = AliPID::EParticleType(species);

    // Get the parameters for this particle type and pt
    const TVectorD *params = GetParams(species, pt);
    
    // IF not in momentum range, NULL is returned --> return default value
    if(!params) return proba;

    Double_t sigma    = (*params)[3];   // sigma value of Gausiian parametrization
    Double_t probLow  = (*params)[6];   // probability to be below eopMin
    Double_t probHigh = (*params)[7];   // probability to be above eopMax

    // get nsigma value for each particle type at this E/p value
    nsigma = GetNumberOfSigmas(pt,eop,type,charge);

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
      // E/p < eopMin  -->  return probability below E/p = eopMin
      if ( nsigma == -99){
	pEMCAL[species] = probLow;
      }
      // E/p > eopMax  -->  return probability above E/p = eopMax
      else if ( nsigma == 99){
	pEMCAL[species] = probHigh;
      }
      // in parametrized region --> calculate probability for corresponding Gauss curve
      else{
	pEMCAL[species]=TMath::Exp(-0.5*(nsigma)*(nsigma))/TMath::Sqrt(2*TMath::Pi()*sigma*sigma);
	
	// normalize to total probability == 1
	pEMCAL[species]*=GetExpectedNorm(pt,type,charge);
      }
    }
    
    // return the electron probability
    proba = pEMCAL[AliPID::kElectron];  

  }

  return proba;

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
const TVectorD* AliEMCALPIDResponse::GetParams(Int_t nParticle, Float_t fPt) const {
  //
  // returns the PID parameters (mean, sigma, probabilities for Hadrons) for a certain particle and pt
  //
  // 0 = momMin
  // 1 = momMax
  // 2 = mean of Gaus
  // 3 = sigma of Gaus
  // 4 = eopLow   
  // 5 = eopHig   
  // 6 = probLow  (not used for electrons)
  // 7 = probHigh (not used for electrons)
  //

  if(nParticle > AliPID::kSPECIES || nParticle <0) return NULL;

  TObjArray * particlePar = dynamic_cast<TObjArray *>(fkPIDParams->At(nParticle));
  if(!particlePar) return NULL;
  
  TIter parIter(particlePar);
  const TVectorD *parameters = NULL;
  Double_t momMin = 0.;
  Double_t momMax = 0.;

  while((parameters = static_cast<const TVectorD *>(parIter()))){

    momMin = (*parameters)[0];
    momMax = (*parameters)[1];

    if( fPt > momMin && fPt < momMax ) return parameters;

  }  
  AliDebug(2, Form("NO params for particle %d and momentum %f \n", nParticle, fPt));

  return parameters;
}
