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

/* $Id$ */

//   Compute PID weights for all the clusters that are in AliESDs.root file
//   the AliESDs.root have to be in the same directory as the class
//
//   and do:    
//   AliEMCALPID *pid = new AliEMCALPID(kFALSE); // this calls the constructor which avoids the call to recparam 
//   pid->SetReconstructor(kFALSE);
//   pid->SetPrintInfo(kTRUE);
//   pid->SetHighFluxParam(); //   pid->SetLowFluxParam(); 
//   
//   then in cluster loop do
//   pid->ComputePID(energy, lambda0);
//  	  
//        Compute PID Weight for all clusters in AliESDs.root file
//	  keep this function for the moment for a simple verification, could be removed
//
//   pid->GetPIDFinal(idx) gives the probabilities
//
//   Double_t PIDFinal[AliPID::kSPECIESN]  is the standard PID for :
//
//
//
//	kElectron :  fPIDFinal[0]
//	kMuon     :  fPIDFinal[1]
//	kPion	  :  fPIDFinal[2]
//	kKaon	  :  fPIDFinal[3]
//	kProton   :  fPIDFinal[4]
//	kPhoton   :  fPIDFinal[5]
//	kPi0	  :  fPIDFinal[6]
//	kNeutron  :  fPIDFinal[7]
//	kKaon0	  :  fPIDFinal[8]
//	kEleCon   :  fPIDFinal[9]
//	kUnknown  :  fPIDFinal[10]
//
//
//    PID[3] is a simple PID for
//      Electron & Photon  PID[0]
//	              Pi0  PID[1]
//		   Hadron  PID[2]
//
// --- standard c ---

// standard C++ includes
//#include <Riostream.h>

// ROOT includes
//#include "TTree.h"
//#include "TVector3.h"
//#include "TBranch.h"
//#include "TClonesArray.h"
//#include "TLorentzVector.h"
#include "TMath.h"
//#include "TRefArray.h"
#include "TArrayD.h"

// STEER includes
#include "AliESDEvent.h"
//#include "AliLog.h"
#include "AliEMCALPID.h"
#include "AliESDCaloCluster.h"
//#include "AliEMCALRecParam.h"
#include "AliEMCALReconstructor.h"

  
ClassImp(AliEMCALPID)
  
//______________________________________________
  AliEMCALPID::AliEMCALPID():
    fPrintInfo(kFALSE), fProbGamma(0.),fProbPiZero(0.),fProbHadron(0.), fWeightHadronEnergy(1.), fWeightGammaEnergy(1.),fWeightPiZeroEnergy(1.),fReconstructor(kTRUE)
{
  //
  // Constructor.
  // Initialize all constant values which have to be used
  // during PID algorithm execution
  //
  
  InitParameters(); 
  
  
}

//______________________________________________
AliEMCALPID::AliEMCALPID(Bool_t reconstructor):
  fPrintInfo(kFALSE), fProbGamma(0.),fProbPiZero(0.),fProbHadron(0.), fWeightHadronEnergy(1.), fWeightGammaEnergy(1.),fWeightPiZeroEnergy(1.),fReconstructor(reconstructor)
{
  //
  // Constructor.
  // Initialize all constant values which have to be used
  // during PID algorithm execution called when used in standalone mode 
  //
  
  InitParameters(); 
  
}

//______________________________________________
void AliEMCALPID::RunPID(AliESDEvent *esd)
{
  //
  // Make the PID for all the EMCAL clusters containedin the ESDs File
  // but just gamma/PiO/Hadron
  //
  // trivial check against NULL object passed
  
  if (esd == 0x0) {
    AliInfo("NULL ESD object passed !!" );
    return ;
  }
  
  Int_t nClusters = esd->GetNumberOfCaloClusters();
  Int_t firstCluster = 0;
  Double_t energy, lambda0;
  for (Int_t iCluster = firstCluster; iCluster < (nClusters + firstCluster); iCluster++) {
    
    AliESDCaloCluster *clust = esd->GetCaloCluster(iCluster);
    if (!clust->IsEMCAL()) continue ; 
    
    energy = clust->E();
    lambda0 = clust->GetM02();
    // verify cluster type
    Int_t clusterType= clust->GetClusterType();
    if (clusterType == AliESDCaloCluster::kEMCALClusterv1 && lambda0 != 0  && energy < 1000) {
      
      //      if (lambda0 != 0  && energy < 1000) {
      
      // reject clusters with lambda0 = 0
      
      
      ComputePID(energy, lambda0);
      
      
      if (fPrintInfo) {
	AliInfo("___________________________________________________");
	AliInfo(Form( "Particle Energy = %f",energy));
	AliInfo(Form( "Particle Lambda0 of the particle = %f", lambda0) );
	AliInfo("PIDWeight of the particle :" );
	AliInfo(Form( " GAMMA  : %f",fPID[0] ));
	AliInfo(Form( " PiZero : %f",fPID[1] ));
	AliInfo(Form( " HADRON : %f", fPID[2] ));
	AliInfo("_________________________________________");
	AliInfo(Form( " kElectron : %f", fPIDFinal[0]) );
	AliInfo(Form( " kMuon     : %f", fPIDFinal[1] ));
	AliInfo(Form( " kPion	    : %f", fPIDFinal[2] ));
	AliInfo(Form( " kKaon	    : %f", fPIDFinal[3] ));
	AliInfo(Form( " kProton   : %f", fPIDFinal[4] ));
	AliInfo(Form( " kPhoton   : %f", fPIDFinal[5] ));
	AliInfo(Form( " kPi0	    : %f", fPIDFinal[6] ));
	AliInfo(Form( " kNeutron  : %f", fPIDFinal[7] ));
	AliInfo(Form( " kKaon0	: %f", fPIDFinal[8] ));
	AliInfo(Form( " kEleCon   : %f", fPIDFinal[9] ));
	AliInfo(Form( " kUnknown  : %f", fPIDFinal[10] ));
	AliInfo("___________________________________________________");
      }
      
      if(fReconstructor){ // In case it is called during reconstruction.
	//	cout << "#############On remplit l esd avec les PIDWeight##########" << endl;
	clust->SetPid(fPIDFinal);}
    } // end if (clusterType...)
  } // end for (iCluster...)
}

//__________________________________________________________
void AliEMCALPID::ComputePID(Double_t energy, Double_t lambda0)
{
//
// This is the main command, which uses the distributions computed and parametrised, 
// and gives the PID by the bayesian method.
//
//   cout << "ENERGY  " <<energy <<" lambda0 "<< lambda0<<  endl;
  
  Double_t weightGammaEnergy  = DistEnergy(energy, 1);
  Double_t weightPiZeroEnergy = DistEnergy(energy, 2);
  Double_t weightHadronEnergy = DistEnergy(energy, 3);
  
  //Double_t weightHadronEnergy = 1.;
  
  Double_t energyhadron=energy;
  if(energyhadron<1.)energyhadron=1.; // no energy dependance of  parametrisation for hadrons below 1 GeV
  if (energy<2){energy =2;} // no energy dependance of parametrisation for gamma and pi0 below 2 GeV
  
  if (energy>55){
    energy =55.;
    energyhadron=55.;
  } // same parametrisation for gamma and hadrons above 55 GeV 
  //   for the pi0 above 55GeV the 2 gammas supperposed no way to distinguish from real gamma  PIDWeight[1]=0
  
  TArrayD paramDistribGamma  = DistLambda0(energy, 1);
  TArrayD paramDistribPiZero = DistLambda0(energy, 2);
  TArrayD paramDistribHadron = DistLambda0(energyhadron, 3);
  
  Bool_t norm = kFALSE;
  
  
  fProbGamma   = TMath::Gaus(lambda0, paramDistribGamma[1], paramDistribGamma[2], norm) * paramDistribGamma[0];
  fProbGamma  += TMath::Landau(((1-paramDistribGamma[4])-lambda0),paramDistribGamma[4],paramDistribGamma[5],norm)* paramDistribGamma[3];
  if(fProbGamma<0.)fProbGamma=0.;
  
  fProbGamma = fProbGamma*weightGammaEnergy;
  
  if(energy>10. || energy < 55.){
    fProbPiZero  = TMath::Gaus(lambda0, paramDistribPiZero[1], paramDistribPiZero[2], norm) * paramDistribPiZero[0];
    fProbPiZero += TMath::Landau(lambda0, paramDistribPiZero[4], paramDistribPiZero[5], norm) * paramDistribPiZero[3];
    if(fProbPiZero<0. || energy<5.)fProbPiZero=0.;
    fProbPiZero = fProbPiZero*weightPiZeroEnergy;
  }
  else {
    fProbPiZero = 0.;
  }
  
  fProbHadron  = TMath::Gaus(lambda0, paramDistribHadron[1], paramDistribHadron[2], norm) * paramDistribHadron[0];
  fProbHadron += TMath::Landau(lambda0, paramDistribHadron[4], paramDistribHadron[5], norm) * paramDistribHadron[3];
  if(fProbHadron<0.)fProbHadron=0.;
  fProbHadron = fProbHadron*weightHadronEnergy; // to take into account the probability for a hadron to have a given reconstructed energy 
  
  // compute PID Weight
  if( (fProbGamma + fProbPiZero + fProbHadron)>0.){
    fPIDWeight[0] = fProbGamma / (fProbGamma + fProbPiZero + fProbHadron);
    fPIDWeight[1] = fProbPiZero / (fProbGamma+fProbPiZero+fProbHadron);
    fPIDWeight[2] = fProbHadron / (fProbGamma+fProbPiZero+fProbHadron);
  }
  else{   
// cases where  energy and lambda0 large,  probably du to 2 clusters folded the clusters PID not assigned to hadron nor Pi0 nor gammas
    fPIDWeight[0] = 0.;
    fPIDWeight[1] = 0.;
    fPIDWeight[2] = 0.;
  }
  
  
  // cout << " PID[0] "<<  fPIDWeight[0] <<  " PID[1] "<<  fPIDWeight[1] <<  " PID[2] "<<  fPIDWeight[2] << endl;
  
  SetPID(fPIDWeight[0], 0);
  SetPID(fPIDWeight[1], 1);
  SetPID(fPIDWeight[2], 2);
  
  // print  pid Weight only for control (= in english ???)
  if (fPrintInfo) {
    AliInfo(Form( "Energy in loop = %f", energy) );
    AliInfo(Form( "Lambda0 in loop = %f", lambda0) );
    AliInfo(Form( "fProbGamma in loop = %f", fProbGamma) );
    AliInfo(Form( "fProbaPiZero = %f", fProbPiZero ));
    AliInfo(Form( "fProbaHadron = %f", fProbHadron) );
    AliInfo(Form( "PIDWeight in loop = %f ||| %f ||| %f",  fPIDWeight[0] , fPIDWeight[1], fPIDWeight[2]) );
    AliInfo("********************************************************" );
  }
  
  fPIDFinal[0]  = fPIDWeight[0]/2; // photon
  fPIDFinal[1]  = fPIDWeight[2]/8;
  fPIDFinal[2]  = fPIDWeight[2]/8;
  fPIDFinal[3]  = fPIDWeight[2]/8;
  fPIDFinal[4]  = fPIDWeight[2]/8;
  fPIDFinal[5]  = fPIDWeight[0]/2; // electron
  fPIDFinal[6]  = fPIDWeight[1]  ; // Pi0
  fPIDFinal[7]  = fPIDWeight[2]/8;
  fPIDFinal[8]  = fPIDWeight[2]/8;
  fPIDFinal[9]  = fPIDWeight[2]/8;
  fPIDFinal[10] = fPIDWeight[2]/8;

}




//________________________________________________________
TArrayD AliEMCALPID::DistLambda0(const Double_t energy, const Int_t type) 
{
  //
  // Compute the values of the parametrised distributions using the data initialised before.
  //
  Double_t constGauss = 0., meanGauss = 0., sigmaGauss = 0.;
  Double_t constLandau=0., mpvLandau=0., sigmaLandau=0.;
  TArrayD  distributionParam(6);
  
  switch (type) {
    
  case 1:
    
    constGauss  = PolynomialMixed2(energy, fGamma[0]);
    meanGauss   = PolynomialMixed2(energy, fGamma[1]);
    sigmaGauss  = PolynomialMixed2(energy, fGamma[2]);
    constLandau = PolynomialMixed2(energy, fGamma[3]);
    mpvLandau   = PolynomialMixed2(energy, fGamma[4]);
    sigmaLandau = PolynomialMixed2(energy, fGamma[5]);
   break;

  case 2:

    constGauss  = PolynomialMixed2(energy, fPiZero[0]);
    meanGauss   = PolynomialMixed2(energy, fPiZero[1]);
    sigmaGauss  = PolynomialMixed2(energy, fPiZero[2]);
    constLandau = PolynomialMixed2(energy, fPiZero[3]);
    mpvLandau   = PolynomialMixed2(energy, fPiZero[4]);
    sigmaLandau = PolynomialMixed2(energy, fPiZero[5]);
    
    break;
  case 3:
    
    constGauss  = PolynomialMixed2(energy, fHadron[0]);
    meanGauss   = PolynomialMixed2(energy, fHadron[1]);
    sigmaGauss  = PolynomialMixed2(energy, fHadron[2]);
    constLandau = PolynomialMixed2(energy, fHadron[3]);
    mpvLandau   = PolynomialMixed2(energy, fHadron[4]);
    sigmaLandau = PolynomialMixed2(energy, fHadron[5]);

    break;
  }
  
  distributionParam[0] = constGauss;
  distributionParam[1] = meanGauss;
  distributionParam[2] = sigmaGauss;
  distributionParam[3] = constLandau;
  distributionParam[4] = mpvLandau;
  distributionParam[5] = sigmaLandau;
  
  return distributionParam;
}

//________________________________________________________
Double_t AliEMCALPID::DistEnergy(const Double_t energy, const Int_t type) 
{
  //
  // Compute the values of the weigh for a given energy the parametrised distribution using the data initialised before.
  //
  Double_t constante = 0.;
  Double_t  energyParam;
  
  switch (type) {
    
  case 1:  
    constante  = 1.;    
    break;
  case 2:
      constante  = 1.;
    break;
  case 3:
    constante  = PowerExp(energy, fHadronEnergyProb);
    break;
  }
  
  energyParam = constante;
  
  // //   cout << "Weight   " << constante << " for energy  "<< energy<< " GeV "<<  endl;
  
  return energyParam;
}


//_______________________________________________________
Double_t AliEMCALPID::Polynomial(const Double_t x, const Double_t *params) const
{
  //
  // Compute a polynomial for a given value of 'x'
  // with the array of parameters passed as the second arg
  //
  
  Double_t y;
  y  = params[0];
  y += params[1] * x;
  y += params[2] * x * x;
  y += params[3] * x * x * x;
  y += params[4] * x * x * x * x;
  y += params[5] * x * x * x * x * x;
  
  return y;
}
//_______________________________________________________
Double_t AliEMCALPID::Polynomial0(const Double_t *params) const 
{
  //
  // Compute a polynomial for a given value of 'x'
  // with the array of parameters passed as the second arg
  //
  
  Double_t y;
  y  = params[0];
  return y;
}

//_______________________________________________________
Double_t AliEMCALPID::Polynomialinv(const Double_t x, const Double_t *params) const
{
  //
  // Compute a polynomial for a given value of 'x'
  // with the array of parameters passed as the second arg
  //
  
  Double_t y;
  if(x>0){
  y  = params[0];
  y += params[1] / x;
  y += params[2] / (x * x);
  y += params[3] / (x * x * x);
  y += params[4] / (x * x * x * x);
  y += params[5] / (x * x * x * x * x);
  }  
  else
    y=0.;
  return y;
  
}
//_______________________________________________________
Double_t AliEMCALPID::PolynomialMixed1(const Double_t x, const Double_t *params) const 
{
  //
  // Compute a polynomial for a given value of 'x'
  // with the array of parameters passed as the second arg
  //
  
  Double_t y;
  if(x>0){
    y  = params[0] / x;
    y += params[1] ;
    y += params[2] * x ;
    //   y += params[3] * 0.;
    //   y += params[4] * 0.;
    //   y += params[5] * 0.;
  }  
  else
    y=0.;
  
  return y;
  
}

//_______________________________________________________
Double_t AliEMCALPID::PolynomialMixed2(const Double_t x, const Double_t *params) const 
{
  //
  // Compute a polynomial for a given value of 'x'
  // with the array of parameters passed as the second arg
  //
  
  Double_t y;
  if(x>0){
    y  = params[0] / ( x * x);
    y += params[1] / x;
    y += params[2] ;
    y += params[3] * x ;
    y += params[4] * x * x ;
    //   y += params[5] * 0.;
  }  
  else
    y=0.;
  //   cout << "y = " << y << endl;
  return y;
  
}

//_______________________________________________________
Double_t AliEMCALPID::PowerExp(const Double_t x, const Double_t *params) const 
{
  //
  // Compute a polynomial for a given value of 'x'
  // with the array of parameters passed as the second arg
  // par[0]*TMath::Power(x[0],par[1])
  // par[0]*TMath::Exp((x[0]-par[1])*par[2]);
  
  Double_t y;
  
  y  = params[0] *TMath::Power( x,params[1]);
  y += params[2] *TMath::Exp((x-params[3])*params[4]);
  
  return y;
  
}


//_______________________________________________________
void AliEMCALPID::InitParameters()
{
  // Initialize PID parameters, depending on the use or not of the reconstructor
  // and the kind of event type if the reconstructor is not used.
  //  fWeightHadronEnergy=0.;
  //  fWeightPiZeroEnergy=0.;
  //  fWeightGammaEnergy=0.;
  
  fPIDWeight[0] = -1;
  fPIDWeight[1] = -1;
  fPIDWeight[2] = -1;
  
  for(Int_t i=0; i<AliPID::kSPECIESN+1; i++)
    fPIDFinal[i]= 0;
  
  const AliEMCALRecParam* recParam = AliEMCALReconstructor::GetRecParam();
  
  if(fReconstructor){
    
    if(!recParam) {
      AliFatal("Reconstruction parameters for EMCAL not set!");
    }
    else {
      
      for(Int_t i=0; i<6; i++){
	for(Int_t j=0; j<6; j++){
	  fGamma[i][j]       = recParam->GetGamma(i,j);
	  fGamma1to10[i][j]  = recParam->GetGamma1to10(i,j);
	  fHadron[i][j]      = recParam->GetHadron(i,j);
	  fHadron1to10[i][j] = recParam->GetHadron1to10(i,j);
	  fPiZero[i][j]      = recParam->GetPiZero(i,j);
	  
	  
	  // 	AliDebug(1,Form("PID parameters (%d, %d): fGamma=%.3f, fPi=%.3f, fHadron=%.3f",
	  // 			i,j, fGamma[i][j],fPiZero[i][j],fHadron[i][j] ));
	  // 	cout << "PID parameters (" << i << " ,"<<j<<") fGamma= "<<  fGamma[i][j]<<" fPi0 ="<<  fPiZero[i][j]<< endl;
	  
	} // end loop j
	fHadronEnergyProb[i] = recParam->GetHadronEnergyProb(i);
	fPiZeroEnergyProb[i] = recParam->GetPiZeroEnergyProb(i);
	fGammaEnergyProb[i]  = recParam->GetGammaEnergyProb(i);
      } //end loop i
      
      
    } // end if !recparam 
    
  } 
  
  else{
    //   init the parameters here instead of from loading from recparam
    //   default parameters are PbPb parameters.
    SetHighFluxParam();
    
  }
  
}


//_______________________________________________________
void AliEMCALPID::SetLowFluxParam()
{
  
  // as a first step, all array elements are initialized to 0.0
  Int_t i, j;
  
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      fGamma[i][j]      = fHadron[i][j] =  fPiZero[i][j] = 0.;
      fGamma1to10[i][j] = fHadron1to10[i][j] = 0.;
    }
       fGammaEnergyProb[i]  =  fGammaEnergyProb[i];
       fPiZeroEnergyProb[i] = fPiZeroEnergyProb[i];
       fHadronEnergyProb[i] = fHadronEnergyProb[i];
  }
  
  // New parametrisation for lambda0^2 (=x): f(x) = normLandau*TMath::Landau(x,mpvLandau,widthLandau)+normgaus*TMath::Gaus(x,meangaus,sigmagaus)
  // See AliEMCALPid (index j) refers to the polynomial parameters of the fit of each parameter vs energy
  // pp

  // paramtype[0][j] = norm gauss
  // paramtype[1][j] = mean gaus
  // paramtype[2][j] = sigma gaus
  // paramtype[3][j] = norm landau
  // paramtype[4][j] = mpv landau
  // paramtype[5][j] = sigma landau

  fGamma[0][0] = -7.656908e-01; 
  fGamma[0][1] =  2.352536e-01; 
  fGamma[0][2] =  1.555996e-02;
  fGamma[0][3] =  2.243525e-04;
  fGamma[0][4] = -2.560087e-06;
  
  fGamma[1][0] =  6.500216e+00;
  fGamma[1][1] = -2.564958e-01;
  fGamma[1][2] =  1.967894e-01;
  fGamma[1][3] = -3.982273e-04;
  fGamma[1][4] =  2.797737e-06;

  fGamma[2][0] =  2.416489e+00;
  fGamma[2][1] = -1.601258e-01;
  fGamma[2][2] =  3.126839e-02;
  fGamma[2][3] =  3.387532e-04;
  fGamma[2][4] = -4.089145e-06;

  fGamma[3][0] =  0.;
  fGamma[3][1] = -2.696008e+00;
  fGamma[3][2] =  6.920305e-01;
  fGamma[3][3] = -2.281122e-03;
  fGamma[3][4] =  0.;

  fGamma[4][0] =  2.281564e-01;
  fGamma[4][1] = -7.575040e-02;
  fGamma[4][2] =  3.813423e-01;
  fGamma[4][3] = -1.243854e-04;
  fGamma[4][4] =  1.232045e-06;

  fGamma[5][0] = -3.290107e-01;
  fGamma[5][1] =  3.707545e-02;
  fGamma[5][2] =  2.917397e-03;
  fGamma[5][3] =  4.695306e-05;
  fGamma[5][4] = -3.572981e-07;

  fHadron[0][0] = 9.482243e-01; 
  fHadron[0][1] =  -2.780896e-01; 
  fHadron[0][2] =  2.223507e-02;
  fHadron[0][3] =  7.294263e-04; 
  fHadron[0][4] =  -5.665872e-06;

  fHadron[1][0] = 0.;
  fHadron[1][1] = 0.;
  fHadron[1][2] = 2.483298e-01;
  fHadron[1][3] = 0.;
  fHadron[1][4] = 0.;

  fHadron[2][0] = -5.601199e+00; 
  fHadron[2][1] =  2.097382e+00; 
  fHadron[2][2] = -2.307965e-01;
  fHadron[2][3] =  9.206871e-03;
  fHadron[2][4] = -8.887548e-05;
 
  fHadron[3][0] =  6.543101e+00;
  fHadron[3][1] =  -2.305203e+00;
  fHadron[3][2] =  2.761673e-01; 
  fHadron[3][3] = -5.465855e-03;
  fHadron[3][4] =  2.784329e-05;
 
  fHadron[4][0] = -2.443530e+01;
  fHadron[4][1] =  8.902578e+00 ;
  fHadron[4][2] = -5.265901e-01;
  fHadron[4][3] = 2.549111e-02;
  fHadron[4][4] =  -2.196801e-04; 

  fHadron[5][0] = 2.102007e-01;
  fHadron[5][1] =  -3.844418e-02;
  fHadron[5][2] =  1.234682e-01;
  fHadron[5][3] = -3.866733e-03;
  fHadron[5][4] = 3.362719e-05 ;

  fPiZero[0][0] =  5.072157e-01;
  fPiZero[0][1] = -5.352747e-01;
  fPiZero[0][2] =  8.499259e-02;
  fPiZero[0][3] = -3.687401e-03;
  fPiZero[0][4] =  5.482280e-05;

  fPiZero[1][0] =  4.590137e+02; 
  fPiZero[1][1] = -7.079341e+01;
  fPiZero[1][2] =  4.990735e+00;
  fPiZero[1][3] = -1.241302e-01;
  fPiZero[1][4] =  1.065772e-03;

  fPiZero[2][0] =  1.376415e+02;
  fPiZero[2][1] = -3.031577e+01;
  fPiZero[2][2] =  2.474338e+00;
  fPiZero[2][3] = -6.903410e-02;
  fPiZero[2][4] =  6.244089e-04;

  fPiZero[3][0] = 0.;
  fPiZero[3][1] =  1.145983e+00;
  fPiZero[3][2] = -2.476052e-01;
  fPiZero[3][3] =  1.367373e-02;
  fPiZero[3][4] = 0.;

  fPiZero[4][0] = -2.097586e+02;
  fPiZero[4][1] =  6.300800e+01;
  fPiZero[4][2] = -4.038906e+00;
  fPiZero[4][3] =  1.088543e-01;
  fPiZero[4][4] = -9.362485e-04;

  fPiZero[5][0] = -1.671477e+01; 
  fPiZero[5][1] =  2.995415e+00;
  fPiZero[5][2] = -6.040360e-02;
  fPiZero[5][3] = -6.137459e-04;
  fPiZero[5][4] =  1.847328e-05;
  
  fHadronEnergyProb[0] = 4.767543e-02;
  fHadronEnergyProb[1] = -1.537523e+00;
  fHadronEnergyProb[2] = 2.956727e-01;
  fHadronEnergyProb[3] = -3.051022e+01;
  fHadronEnergyProb[4] =-6.036931e-02;

  Int_t ii= 0;
  Int_t jj= 3;
  AliDebug(1,Form("PID parameters (%d, %d): fGamma=%.3f, fPi=%.3f, fHadron=%.3f",
 			ii,jj, fGamma[ii][jj],fPiZero[ii][jj],fHadron[ii][jj] ));
  //cout << " LowFlux Parameters fGamma [2][2] = " << fGamma[2][2] << endl;
  //cout << " LowFlux Parameters fHadron [2][2] = " << fHadron[2][2] << endl;
   
  // end for proton-proton  

}

//_______________________________________________________
void AliEMCALPID::SetHighFluxParam()
{
  
  // as a first step, all array elements are initialized to 0.0
  Int_t i, j;
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      fGamma[i][j]      = fHadron[i][j] = fPiZero[i][j] = 0.;
      fGamma1to10[i][j] = fHadron1to10[i][j] = 0.;
    }
    fGammaEnergyProb[i]  = 0.;
    fPiZeroEnergyProb[i] = 0.;
    fHadronEnergyProb[i] = 0.;
  }
  
  // Pb Pb  this goes with inverted landau + gaussian for gammas, landau+gaussian for Pi0 and hadrons
  
  fGamma[0][0] = -7.656908e-01; 
  fGamma[0][1] =  2.352536e-01; 
  fGamma[0][2] =  1.555996e-02;
  fGamma[0][3] =  2.243525e-04;
  fGamma[0][4] = -2.560087e-06;
  
  fGamma[1][0] =  6.500216e+00;
  fGamma[1][1] = -2.564958e-01;
  fGamma[1][2] =  1.967894e-01;
  fGamma[1][3] = -3.982273e-04;
  fGamma[1][4] =  2.797737e-06;

  fGamma[2][0] =  2.416489e+00;
  fGamma[2][1] = -1.601258e-01;
  fGamma[2][2] =  3.126839e-02;
  fGamma[2][3] =  3.387532e-04;
  fGamma[2][4] = -4.089145e-06;
 
  fGamma[3][0] =  0.;
  fGamma[3][1] = -2.696008e+00;
  fGamma[3][2] =  6.920305e-01;
  fGamma[3][3] = -2.281122e-03;
  fGamma[3][4] =  0.;

  fGamma[4][0] =  2.281564e-01;
  fGamma[4][1] = -7.575040e-02;
  fGamma[4][2] =  3.813423e-01;
  fGamma[4][3] = -1.243854e-04;
  fGamma[4][4] =  1.232045e-06;

  fGamma[5][0] = -3.290107e-01;
  fGamma[5][1] =  3.707545e-02;
  fGamma[5][2] =  2.917397e-03;
  fGamma[5][3] =  4.695306e-05;
  fGamma[5][4] = -3.572981e-07;
   
  fHadron[0][0] =   1.519112e-01;
  fHadron[0][1] = -8.267603e-02;
  fHadron[0][2] =  1.914574e-02;
  fHadron[0][3] = -2.677921e-04;
  fHadron[0][4] =  5.447939e-06;

  fHadron[1][0] = 0.;
  fHadron[1][1] = -7.549870e-02; 
  fHadron[1][2] = 3.930087e-01;
  fHadron[1][3] = -2.368500e-03; 
  fHadron[1][4] = 0.;

  fHadron[2][0] = 0.;
  fHadron[2][1] =  -2.463152e-02;
  fHadron[2][2] = 1.349257e-01;
  fHadron[2][3] = -1.089440e-03;
  fHadron[2][4] = 0.;

  fHadron[3][0] = 0.;
  fHadron[3][1] = 5.101560e-01;
  fHadron[3][2] = 1.458679e-01;
  fHadron[3][3] = 4.903068e-04;
  fHadron[3][4] = 0.;

  fHadron[4][0] = 0.;
  fHadron[4][1] = -6.693943e-03; 
  fHadron[4][2] =  2.444753e-01;
  fHadron[4][3] = -5.553749e-05;
  fHadron[4][4] = 0.;

  fHadron[5][0] = -4.414030e-01;
  fHadron[5][1] = 2.292277e-01;
  fHadron[5][2] = -2.433737e-02;
  fHadron[5][3] =  1.758422e-03;
  fHadron[5][4] = -3.001493e-05;
  
  fPiZero[0][0] =  5.072157e-01;
  fPiZero[0][1] = -5.352747e-01;
  fPiZero[0][2] =  8.499259e-02;
  fPiZero[0][3] = -3.687401e-03;
  fPiZero[0][4] =  5.482280e-05;
  
  fPiZero[1][0] =  4.590137e+02; 
  fPiZero[1][1] = -7.079341e+01;
  fPiZero[1][2] =  4.990735e+00;
  fPiZero[1][3] = -1.241302e-01;
  fPiZero[1][4] =  1.065772e-03;
  
  fPiZero[2][0] =  1.376415e+02;
  fPiZero[2][1] = -3.031577e+01;
  fPiZero[2][2] =  2.474338e+00;
  fPiZero[2][3] = -6.903410e-02;
  fPiZero[2][4] =  6.244089e-04;

  fPiZero[3][0] = 0.;
  fPiZero[3][1] =  1.145983e+00;
  fPiZero[3][2] = -2.476052e-01;
  fPiZero[3][3] =  1.367373e-02;
  fPiZero[3][4] = 0.;

  fPiZero[4][0] = -2.097586e+02;
  fPiZero[4][1] =  6.300800e+01;
  fPiZero[4][2] = -4.038906e+00;
  fPiZero[4][3] =  1.088543e-01;
  fPiZero[4][4] = -9.362485e-04;

  fPiZero[5][0] = -1.671477e+01; 
  fPiZero[5][1] =  2.995415e+00;
  fPiZero[5][2] = -6.040360e-02;
  fPiZero[5][3] = -6.137459e-04;
  fPiZero[5][4] =  1.847328e-05;

  // those are the High Flux PbPb ones
  fHadronEnergyProb[0] = 0.;
  fHadronEnergyProb[1] = 0.;
  fHadronEnergyProb[2] =  6.188452e-02;
  fHadronEnergyProb[3] =  2.030230e+00;
  fHadronEnergyProb[4] = -6.402242e-02;

 Int_t ii= 0;
 Int_t jj= 3;
 AliDebug(1,Form("PID parameters (%d, %d): fGamma=%.3f, fPi=%.3f, fHadron=%.3f",
 			ii,jj, fGamma[ii][jj],fPiZero[ii][jj],fHadron[ii][jj] ));
  //cout << " HighFlux Parameters fGamma [2][2] = " << fGamma[2][2] << endl;
  //cout << " HighFlux Parameters fHadron [2][2] = " << fHadron[2][2] << endl;
   
}
