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

//-----------------------------------------------------------------------------
// Container of EMCAL reconstruction parameters
// The purpose of this object is to store it to OCDB
// and retrieve it in the corresponding reconstruction class:
// AliEMCALClusterizer, AliEMCALPID, AliEMCALTracker ...
//
// Author: Yuri Kharlov
//-----------------------------------------------------------------------------

// --- Root header files
//#include "TObjArray.h"

// --- AliRoot header files ---
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliEMCALRecParam.h"

ClassImp(AliEMCALRecParam)
  
TObjArray* AliEMCALRecParam::fgkMaps =0; //ALTRO mappings 

AliEMCALRecParam::AliEMCALRecParam() :
  AliDetectorRecoParam(),
  fClusteringThreshold(0.5),
  fW0(4.5),
  fMinECut(0.045), 
  fUnfold(kFALSE), 
  fLocMaxCut(0.03), 
  fTimeCut(1.),// high value, accept all
  fTimeMin(-1.),// small value, accept all
  fTimeMax(1.),// high value, accept all//clustering
  fClusterizerFlag(AliEMCALRecParam::kClusterizerv1),
  fNRowDiff(1),
  fNColDiff(1),
  fMthCutEta(0.025), 
  fMthCutPhi(0.05),
  fStep(50),
  fTrkCutPt(0.0),
  fTrkCutNITS(0.0),
  fTrkCutNTPC(50.0), //track matching
  fHighLowGainFactor(16.0), 
  fOrderParameter(2), 
  fTau(2.35), 
  fNoiseThreshold(3), 
  fNPedSamples(5), 
  fRemoveBadChannels(kFALSE),
  fFittingAlgorithm(0), 
  fUseFALTRO(kTRUE), 
  fFitLEDEvents(kFALSE)//raw signal
{
  // default reco values
  
  // PID parameters for Pb Pb from Lambda0 distributions fitted by  
  // a landau inverted + Gaussian for Gammas
  // and a Landau +Gaussian for Pi0 and hadrons 
  // New parametrisation for 
  // lambda0^2 (=x): f(x) = normLandau*TMath::Landau(((1-mpvlandau)-x),mpvLandau,widthLandau)+normgaus*TMath::Gaus(x,meangaus,sigmagaus) for gammas
  // lambda0^2 (=x): f(x) = normLandau*TMath::Landau(x,mpvLandau,widthLandau)+normgaus*TMath::Gaus(x,meangaus,sigmagaus) for pi0 & hadrons
  
  // See AliEMCALPid 
  // (index i) refers to each parameters of the f(lambda0^2) 
  // i=0: normGaus
  // i=1: meanGaus
  // i=2: sigmaGaus
  // i=3: normLandau
  // i=4: mpvLandau
  // i=5: sigmaLanda
  // (index j) refers to the polynomial parameters of the fit of each parameter vs energy
  // Pb Pb
  
  // as a first step, all array elements are initialized to 0.0
  Int_t i=0, j=0;
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {      
      fGamma[i][j] =  fPiZero[i][j] = fHadron[i][j] = 0.; 
      fGamma1to10[i][j] =  fHadron1to10[i][j]= 0.;
    }
    fGammaEnergyProb[i] =0.; // not yet implemented
    fHadronEnergyProb[i]=0.; 
    fPiZeroEnergyProb[i]=0.; // not yet implemented
    
    
  }
  // Set here default parameters for Pb+Pb (high flux)
  
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
    
  fHadronEnergyProb[0]= 0.;
  fHadronEnergyProb[1]= 0.;
  fHadronEnergyProb[2]=  6.188452e-02;
  fHadronEnergyProb[3]=  2.030230e+00;
  fHadronEnergyProb[4]= -6.402242e-02;
  
  //unfolding  
  fSSPars[0] = 0.9262;
  fSSPars[1] = 3.365;
  fSSPars[2] = 1.548;
  fSSPars[3] = 0.1625;
  fSSPars[4] = -0.4195;
  fSSPars[5] = 0.;
  fSSPars[6] = 0.;
  fSSPars[7] = 2.332;
  fPar5[0] = 12.31;
  fPar5[1] = -0.007381;
  fPar5[2] = -0.06936;
  fPar6[0] = 0.05452; 
  fPar6[1] = 0.0001228; 
  fPar6[2] = 0.001361; 

}

//-----------------------------------------------------------------------------
AliEMCALRecParam::AliEMCALRecParam(const AliEMCALRecParam& rp) :
  AliDetectorRecoParam(),
  fClusteringThreshold(rp.fClusteringThreshold),
  fW0(rp.fW0),
  fMinECut(rp.fMinECut), 
  fUnfold(rp.fUnfold), 
  fLocMaxCut(rp.fLocMaxCut), 
  fTimeCut(rp.fTimeCut), 
  fTimeMin(rp.fTimeMin),
  fTimeMax(rp.fTimeMax),//clustering
  fClusterizerFlag(rp.fClusterizerFlag),
  fNRowDiff(rp.fNRowDiff),
  fNColDiff(rp.fNColDiff),
  fMthCutEta(rp.fMthCutEta), 
  fMthCutPhi(rp.fMthCutPhi),
  fStep(rp.fStep),
  fTrkCutPt(rp.fTrkCutPt),
  fTrkCutNITS(rp.fTrkCutNITS),
  fTrkCutNTPC(rp.fTrkCutNTPC), // track matching
  fHighLowGainFactor(rp.fHighLowGainFactor), 
  fOrderParameter(rp.fOrderParameter), 
  fTau(rp.fTau), 
  fNoiseThreshold(rp.fNoiseThreshold), 
  fNPedSamples(rp.fNPedSamples), 	
  fRemoveBadChannels(rp.fRemoveBadChannels),
  fFittingAlgorithm(rp.fFittingAlgorithm),  
  fUseFALTRO(rp.fUseFALTRO),
  fFitLEDEvents(rp.fFitLEDEvents) //raw signal
{
  //copy constructor
  
  //PID values
  Int_t i=0, j=0;
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      fGamma[i][j]       = rp.fGamma[i][j];
      fGamma1to10[i][j]  = rp.fGamma1to10[i][j];
      fHadron[i][j]      = rp.fHadron[i][j];
      fHadron1to10[i][j] = rp.fHadron1to10[i][j];
      fPiZero[i][j]      = rp.fPiZero[i][j];
    }
    fGammaEnergyProb[i]  = rp.fGammaEnergyProb[i];
    fPiZeroEnergyProb[i] = rp.fPiZeroEnergyProb[i];
    fHadronEnergyProb[i] = rp.fHadronEnergyProb[i];
    
  }
  
  //unfolding  
  for (i = 0; i < 8; i++) {
    fSSPars[i] = rp.fSSPars[i];
  }
  for (i = 0; i < 3; i++) {
    fPar5[i] = rp.fPar5[i];
    fPar6[i] = rp.fPar6[i];
  }

}

//-----------------------------------------------------------------------------
AliEMCALRecParam& AliEMCALRecParam::operator = (const AliEMCALRecParam& rp)
{
  //assignment operator
  
  if(this != &rp) {
    fClusteringThreshold = rp.fClusteringThreshold;
    fW0        = rp.fW0;
    fMinECut   = rp.fMinECut;
    fUnfold    = rp.fUnfold;
    fLocMaxCut = rp.fLocMaxCut; 
    fTimeCut   = rp.fTimeCut;
    fTimeMax   = rp.fTimeMax;
    fTimeMin   = rp.fTimeMin;//clustering
    fClusterizerFlag   = rp.fClusterizerFlag;
    fNRowDiff  = rp.fNRowDiff;
    fNColDiff  = rp.fNColDiff;
    fMthCutEta         = rp.fMthCutEta;
    fMthCutPhi         = rp.fMthCutPhi;
    fStep              = rp.fStep;
    fTrkCutPt          = rp.fTrkCutPt;
    fTrkCutNITS        = rp.fTrkCutNITS;
    fTrkCutNTPC        = rp.fTrkCutNTPC; //track matching
    fHighLowGainFactor = rp.fHighLowGainFactor; 
    fOrderParameter    = rp.fOrderParameter;
    fTau               = rp.fTau;
    fNoiseThreshold    = rp.fNoiseThreshold;
    fNPedSamples       = rp.fNPedSamples; 
    fRemoveBadChannels = rp.fRemoveBadChannels;
    fFittingAlgorithm  = rp.fFittingAlgorithm;
    fUseFALTRO         = rp.fUseFALTRO;
    fFitLEDEvents      = rp.fFitLEDEvents;//raw signal
	  
    //PID values
    Int_t i=0, j=0;
    for (i = 0; i < 6; i++) {
      for (j = 0; j < 6; j++) {
	fGamma[i][j]       = rp.fGamma[i][j];
	fGamma1to10[i][j]  = rp.fGamma1to10[i][j];
	fHadron[i][j]      = rp.fHadron[i][j];
	fHadron1to10[i][j] = rp.fHadron1to10[i][j];
	fPiZero[i][j]      = rp.fPiZero[i][j];
      }
      fGammaEnergyProb[i]  = rp.fGammaEnergyProb[i];
      fPiZeroEnergyProb[i] = rp.fPiZeroEnergyProb[i];
      fHadronEnergyProb[i] = rp.fHadronEnergyProb[i];
    }
    //unfolding  
    for (i = 0; i < 8; i++) {
      fSSPars[i] = rp.fSSPars[i];
    }
    for (i = 0; i < 3; i++) {
      fPar5[i] = rp.fPar5[i];
      fPar6[i] = rp.fPar6[i];
    }

  }    
  
  return *this;
  
}

//-----------------------------------------------------------------------------
AliEMCALRecParam* AliEMCALRecParam::GetDefaultParameters()
{
  //default parameters for the reconstruction
  AliEMCALRecParam* params = GetLowFluxParam();
  params->SetName("Default - p+p");
  params->SetTitle("Default - p+p");
  return params;
  
}

//-----------------------------------------------------------------------------
AliEMCALRecParam* AliEMCALRecParam::GetCalibParam()
{
  //parameters for the reconstruction of calibration runs
  AliEMCALRecParam* params = GetLowFluxParam();
  //params->SetClusteringThreshold(0.1); // 100 MeV                                             
  //params->SetMinECut(0.01);  //10 MeV       
  params->SetName("Calibration - LED");
  params->SetTitle("Calibration - LED");
  params->SetEventSpecie(AliRecoParam::kCalib);
  
  return params;
  
}

//-----------------------------------------------------------------------------
AliEMCALRecParam* AliEMCALRecParam::GetCosmicParam()
{
  //parameters for the reconstruction of cosmic runs
  AliEMCALRecParam* params = GetLowFluxParam();
  //params->SetClusteringThreshold(0.1); // 100 MeV                                             
  //params->SetMinECut(0.01);  //10 MeV       
  params->SetName("Cosmic");
  params->SetTitle("Cosmic");
  params->SetEventSpecie(AliRecoParam::kCosmic);
  
  return params;
  
}


//-----------------------------------------------------------------------------
AliEMCALRecParam* AliEMCALRecParam::GetLowFluxParam()
{
  //low flux/multiplicity ("p+p") parameters for the reconstruction
  AliEMCALRecParam* params = new AliEMCALRecParam();
  params->SetClusteringThreshold(0.1); // 100 MeV                                             
  params->SetMinECut(0.01);  //10 MeV       
  params->SetName("Low Flux - p+p");
  params->SetTitle("Low Flux - p+p");
  params->SetEventSpecie(AliRecoParam::kLowMult);
  params->SetExtrapolateStep(1);
  
  
  //PID parameters for pp  implemented 
  // as a first step, all array elements are initialized to 0.0
  Int_t i=0, j=0;
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      params->SetGamma(i,j,0.);
      params->SetGamma1to10(i,j,0.);
      params->SetHadron(i,j,0.);
      params->SetHadron1to10(i,j,0.);
      params->SetPiZero(i,j,0.);
      
    }
    params->SetGammaEnergyProb(i,0.); // not yet implemented
    params->SetHadronEnergyProb(i,0.);
    params->SetPiZeroEnergyProb(i,0.); // not yet implemented
  }
  
  
  params->SetGamma(0,0,-7.656908e-01);
  params->SetGamma(0,1,2.352536e-01);
  params->SetGamma(0,2,1.555996e-02);
  params->SetGamma(0,3,2.243525e-04);
  params->SetGamma(0,4,-2.560087e-06);

  params->SetGamma(1,0,6.500216e+00);
  params->SetGamma(1,1,-2.564958e-01);
  params->SetGamma(1,2,1.967894e-01);
  params->SetGamma(1,3,-3.982273e-04);
  params->SetGamma(1,4,2.797737e-06);
  
  params->SetGamma(2,0,2.416489e+00);
  params->SetGamma(2,1,-1.601258e-01);
  params->SetGamma(2,2,3.126839e-02);
  params->SetGamma(2,3,3.387532e-04);
  params->SetGamma(2,4,-4.089145e-06);
 
  params->SetGamma(3,0,0.);
  params->SetGamma(3,1,-2.696008e+00); 
  params->SetGamma(3,2, 6.920305e-01);
  params->SetGamma(3,3,-2.281122e-03);
  params->SetGamma(3,4,0.);

  params->SetGamma(4,0,2.281564e-01); 
  params->SetGamma(4,1,-7.575040e-02);
  params->SetGamma(4,2,3.813423e-01);
  params->SetGamma(4,3,-1.243854e-04);
  params->SetGamma(4,4,1.232045e-06);

  params->SetGamma(5,0,-3.290107e-01);
  params->SetGamma(5,1,3.707545e-02);
  params->SetGamma(5,2,2.917397e-03);
  params->SetGamma(5,3,4.695306e-05);
  params->SetGamma(5,4,-3.572981e-07);

  params->SetHadron(0,0,9.482243e-01); 
  params->SetHadron(0,1,-2.780896e-01);
  params->SetHadron(0,2, 2.223507e-02);
  params->SetHadron(0,3,7.294263e-04);
  params->SetHadron(0,4,-5.665872e-06); 

  params->SetHadron(1,0,0.);
  params->SetHadron(1,1,0.);
  params->SetHadron(1,2,2.483298e-01);
  params->SetHadron(1,3,0.);
  params->SetHadron(1,4,0.);

  params->SetHadron(2,0,-5.601199e+00);
  params->SetHadron(2,1,2.097382e+00);
  params->SetHadron(2,2,-2.307965e-01);
  params->SetHadron(2,3,9.206871e-03);
  params->SetHadron(2,4,-8.887548e-05);

  params->SetHadron(3,0,6.543101e+00);
  params->SetHadron(3,1,-2.305203e+00);
  params->SetHadron(3,2,2.761673e-01);
  params->SetHadron(3,3,-5.465855e-03);
  params->SetHadron(3,4,2.784329e-05);

  params->SetHadron(4,0,-2.443530e+01);
  params->SetHadron(4,1,8.902578e+00);
  params->SetHadron(4,2,-5.265901e-01);
  params->SetHadron(4,3,2.549111e-02);
  params->SetHadron(4,4,-2.196801e-04);

  params->SetHadron(5,0,2.102007e-01);
  params->SetHadron(5,1,-3.844418e-02);
  params->SetHadron(5,2,1.234682e-01);
  params->SetHadron(5,3,-3.866733e-03);
  params->SetHadron(5,4,3.362719e-05);

  params->SetPiZero(0,0,5.07215e-01);
  params->SetPiZero(0,1,-5.35274e-01);
  params->SetPiZero(0,2,8.49925e-02);
  params->SetPiZero(0,3,-3.68740e-03);
  params->SetPiZero(0,4,5.48228e-05);

  params->SetPiZero(1,0,4.590137e+02);
  params->SetPiZero(1,1,-7.079341e+01);
  params->SetPiZero(1,2,4.990735e+00);
  params->SetPiZero(1,3,-1.241302e-01);
  params->SetPiZero(1,4,1.065772e-03);

  params->SetPiZero(2,0,1.376415e+02); 
  params->SetPiZero(2,1,-3.031577e+01);
  params->SetPiZero(2,2,2.474338e+00);
  params->SetPiZero(2,3,-6.903410e-02); 
  params->SetPiZero(2,4,6.244089e-04);

  params->SetPiZero(3,0,0.);
  params->SetPiZero(3,1,1.145983e+00);
  params->SetPiZero(3,2,-2.476052e-01);
  params->SetPiZero(3,3,1.367373e-02);
  params->SetPiZero(3,4,0.);

  params->SetPiZero(4,0,-2.09758e+02);
  params->SetPiZero(4,1,6.30080e+01);
  params->SetPiZero(4,2,-4.03890e+00);
  params->SetPiZero(4,3,1.08854e-01);
  params->SetPiZero(4,4,-9.36248e-04);

  params->SetPiZero(5,0,-1.671477e+01);
  params->SetPiZero(5,1,2.995415e+00);
  params->SetPiZero(5,2,-6.040360e-02);
  params->SetPiZero(5,3,-6.137459e-04);
  params->SetPiZero(5,4,1.847328e-05);

 
//     params->SetHadronEnergyProb(0,0.);
//     params->SetHadronEnergyProb(1,0.);
//     params->SetHadronEnergyProb(2,1.);
//     params->SetHadronEnergyProb(3,0.);
//     params->SetHadronEnergyProb(4,0.);
 
  params->SetHadronEnergyProb(0, 4.767543e-02);
  params->SetHadronEnergyProb(1,-1.537523e+00);
  params->SetHadronEnergyProb(2,2.956727e-01);
  params->SetHadronEnergyProb(3,-3.051022e+01);
  params->SetHadronEnergyProb(4,-6.036931e-02);

//   Int_t ii= 0;
//   Int_t jj= 3;
// 	AliDebug(1,Form("PID parameters (%d, %d): fGamma=%.3f, fPi=%.3f, fHadron=%.3f",
//  			ii,jj, params->GetGamma(ii,jj), params->GetPiZero(ii,jj), params->GetHadron(ii,jj)));
// 	cout << " Low Flux Parameters fGamma [2][2] = " << params->GetGamma(2,2) << endl;
// 	cout << " Low Flux Parameters fHadron [2][2] = " << params->GetHadron(2,2) << endl;
   
  return params;

}


//-----------------------------------------------------------------------------
AliEMCALRecParam* AliEMCALRecParam::GetHighFluxParam()
{
  //high flux/multiplicity ("Pb+Pb") parameters for the reconstruction
  AliEMCALRecParam* params = new AliEMCALRecParam();
  //For now, same as default
  //if later these need to be modified, here's where it is done
  params->SetName("High Flux - Pb+Pb");
  params->SetTitle("High Flux - Pb+Pb");
  params->SetEventSpecie(AliRecoParam::kHighMult);
  params->SetTrkCutPt(0.15);//This value can be higher if necessary
  
  return params;
  
}

//-----------------------------------------------------------------------------
void AliEMCALRecParam::Print(Option_t * opt) const
{
  // Print reconstruction parameters to stdout
  // if nothing is specified print all, if "reco" just clusterization/track matching
  // if "pid", just PID parameters, if "raw", just raw utils parameters.
  if(!strcmp("",opt) || !strcmp("reco",opt)){
    AliInfo(Form("Clusterizer selected: %d", fClusterizerFlag));
    AliInfo(Form("Clusterization parameters :\n fClusteringThreshold=%.3f,\n fW0=%.3f,\n fMinECut=%.3f,\n fUnfold=%d,\n fLocMaxCut=%.3f,\n fTimeCut=%2.1f ns\n fTimeMin=%2.1f ns\n fTimeMax=%2.1f ns\n",
		 fClusteringThreshold,fW0,fMinECut,fUnfold,fLocMaxCut,fTimeCut*1.e9,fTimeMin*1e9,fTimeMax*1e9));
    
    AliInfo(Form("Track-matching cuts :\n dEta<%f, dPhi<%f, step=%f[cm], pT>%f, NITS>%f, NTPC>%f\n", 
		 fMthCutEta, fMthCutPhi, fStep, fTrkCutPt, fTrkCutNITS,fTrkCutNTPC));

    AliInfo(Form("Unfolding parameters, Shower shape function :\n")); 
    for(Int_t i = 0; i < 8; i++){
	printf(" %f, ", fSSPars[i]);
    }
    printf("\n Parameter 5 : ");
    for(Int_t i = 0; i < 3; i++){
      printf(" %f, ", fPar5[i]);
    }
    printf("\n Parameter 6 : ");
    for(Int_t i = 0; i < 3; i++){
      printf(" %f, ", fPar6[i]);
    }
    printf("\n");
  }
  
  if(!strcmp("",opt) || !strcmp("pid",opt)){
    AliInfo(Form("PID parameters, Gamma :\n"));
    for(Int_t i = 0; i < 6; i++){
      for(Int_t j = 0; j < 6; j++){
	printf(" %f, ", fGamma[i][j]);
      }
      printf("\n");
    }
    
    
    AliInfo(Form("PID parameters, Hadron :\n"));
    for(Int_t i = 0; i < 6; i++){
      for(Int_t j = 0; j < 6; j++){
	printf(" %f, ", fHadron[i][j]);
      }
      printf("\n");
    }
    
    printf("\n");
    
    AliInfo(Form("PID parameters, Pi0zero :\n"));
    for(Int_t i = 0; i < 6; i++){
      for(Int_t j = 0; j < 6; j++){
	printf(" %f, ", fPiZero[i][j]);
      }
      printf("\n");
    }
    
    printf("\n");
    
  }

  if(!strcmp("",opt) || !strcmp("raw",opt)){
    AliInfo(Form("Raw signal parameters: \n gain factor=%f, order=%d, tau=%f, noise threshold=%d, nped samples=%d \n",
		 fHighLowGainFactor,fOrderParameter,fTau,fNoiseThreshold,fNPedSamples));
    AliInfo(Form("Raw signal: remove bad channels? %d, \n \t with fitting algorithm %d, \n \t Use FALTRO %d, Fit LED events %d \n",
		 fRemoveBadChannels, fFittingAlgorithm, fUseFALTRO, fFitLEDEvents));
  }
}

//-----------------------------------------------------------------------------
const TObjArray* AliEMCALRecParam::GetMappings()
{
  //Returns array of AliAltroMappings for RCU0..RCUX.
  //If not found, read it from OCDB.                                           
  
  //Quick check as follows:                                                   
  //  root [0] AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT"
  //  root [1] AliCDBManager::Instance()->SetRun(1);                             
  //  root [2] TObjArray* maps = AliEMCALRecParam::GetMappings();                
  //  root [3] maps->Print();                                                    
  
  if(fgkMaps) return fgkMaps;
  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("EMCAL/Calib/Mapping");
  if(entry)
    fgkMaps = (TObjArray*)entry->GetObject();
  
  return fgkMaps;
  
}

