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
// --- AliRoot header files ---
#include "TObjArray.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliEMCALRecParam.h"
#include "AliLog.h"

ClassImp(AliEMCALRecParam)

TObjArray* AliEMCALRecParam::fgkMaps =0; //ALTRO mappings 

//-----------------------------------------------------------------------------
// Container of EMCAL reconstruction parameters
// The purpose of this object is to store it to OCDB
// and retrieve it in AliEMCALClusterizerv1
// Author: Yuri Kharlov
//-----------------------------------------------------------------------------

AliEMCALRecParam::AliEMCALRecParam() :
  AliDetectorRecoParam(),
  fClusteringThreshold(0.5),
  fW0(4.5),
  fMinECut(0.45), 
  fUnfold(kFALSE), 
  fLocMaxCut(0.03), //clustering
  fTrkCutX(6.0), 
  fTrkCutY(6.0), 
  fTrkCutZ(6.0),  
  fTrkCutR(10.0),
  fTrkCutAlphaMin(-50.0), 
  fTrkCutAlphaMax(50.0), 
  fTrkCutAngle(10000.0),
  fTrkCutNITS(3.0),
  fTrkCutNTPC(20.0), //track matching
  fHighLowGainFactor(16.0), 
  fOrderParameter(2), 
  fTau(2.35), 
  fNoiseThreshold(3), 
  fNPedSamples(5) //raw signal
{
  // default reco values

  //PID parameters (Guenole)

  // as a first step, all array elements are initialized to 0.0
  Int_t i, j;
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      fGamma[i][j] = fHadron[i][j] = fPiZero5to10[i][j] = fPiZero10to60[i][j] = 0.;
    }
  }

  // then, only the ones which must be not zero are initialized
  // while the others will remain to the value 0.0

  fGamma[0][0] =  0.038022;
  fGamma[0][1] = -0.0001883;
  fGamma[0][2] =  5.449e-06;

  fGamma[1][0] =  0.207313;
  fGamma[1][1] = -0.000978;
  fGamma[1][2] =  0.00001634;

  fGamma[2][0] =  0.043364;
  fGamma[2][1] = -0.0002048;
  fGamma[2][2] =  8.661e-06;
  fGamma[2][3] = -1.353e-07;

  fGamma[3][0] =  0.265004;
  fGamma[3][1] =  0.061298;
  fGamma[3][2] = -0.003203;
  fGamma[3][3] =  4.73e-05;

  fGamma[4][0] =  0.243579;
  fGamma[4][1] = -1.614e-05;

  fGamma[5][0] =  0.002942;
  fGamma[5][1] = -3.976e-05;

  fHadron[0][0] =  0.011945 / 3.;
  fHadron[0][1] =  0.000386 / 3.;
  fHadron[0][2] = -0.000014 / 3.;
  fHadron[0][3] =  1.336e-07 / 3.;

  fHadron[1][0] =  0.496544;
  fHadron[1][1] = -0.003226;
  fHadron[1][2] =  0.00001678;

  fHadron[2][0] =  0.144838;
  fHadron[2][1] = -0.002954;
  fHadron[2][2] =  0.00008754;
  fHadron[2][3] = -7.587e-07;

  fHadron[3][0] =  1.264461 / 7.;
  fHadron[3][1] =  0.002097 / 7.;

  fHadron[4][0] =  0.261950;
  fHadron[4][1] = -0.001078;
  fHadron[4][2] =  0.00003237;
  fHadron[4][3] = -3.241e-07;
  fHadron[4][4] =  0.;
  fHadron[4][5] =  0.;
  fHadron[5][0] =  0.010317;
  fHadron[5][1] =  0.;
  fHadron[5][2] =  0.;
  fHadron[5][3] =  0.;
  fHadron[5][4] =  0.;
  fHadron[5][5] =  0.;

  fPiZero5to10[0][0] = 0.009138;
  fPiZero5to10[0][1] = 0.0006377;

  fPiZero5to10[1][0] = 0.08;

  fPiZero5to10[2][0] = -0.061119;
  fPiZero5to10[2][1] =  0.019013;

  fPiZero5to10[3][0] =  0.2;

  fPiZero5to10[4][0] =  0.252044;
  fPiZero5to10[4][1] = -0.002315;

  fPiZero5to10[5][0] =  0.002942;
  fPiZero5to10[5][1] = -3.976e-05;
  
  fPiZero10to60[0][0] =  0.009138;
  fPiZero10to60[0][1] =  0.0006377;

  fPiZero10to60[1][0] =  1.272837;
  fPiZero10to60[1][1] = -0.069708;
  fPiZero10to60[1][2] =  0.001568;
  fPiZero10to60[1][3] = -1.162e-05;

  fPiZero10to60[2][0] =  0.139703;
  fPiZero10to60[2][1] =  0.003687;
  fPiZero10to60[2][2] = -0.000568;
  fPiZero10to60[2][3] =  1.498e-05;
  fPiZero10to60[2][4] = -1.174e-07;

  fPiZero10to60[3][0] = -0.826367;
  fPiZero10to60[3][1] =  0.096951;
  fPiZero10to60[3][2] = -0.002215;
  fPiZero10to60[3][3] =  2.523e-05;

  fPiZero10to60[4][0] =  0.249890;
  fPiZero10to60[4][1] = -0.000063;

  fPiZero10to60[5][0] =  0.002942;
  fPiZero10to60[5][1] = -3.976e-05;

}

//-----------------------------------------------------------------------------
AliEMCALRecParam::AliEMCALRecParam(const AliEMCALRecParam& rp) :
  AliDetectorRecoParam(),
  fClusteringThreshold(rp.fClusteringThreshold),
  fW0(rp.fW0),
  fMinECut(rp.fMinECut), 
  fUnfold(rp.fUnfold), 
  fLocMaxCut(rp.fLocMaxCut), //clustering
  fTrkCutX(rp.fTrkCutX), 
  fTrkCutY(rp.fTrkCutY), 
  fTrkCutZ(rp.fTrkCutZ),  
  fTrkCutR(rp.fTrkCutR),
  fTrkCutAlphaMin(rp.fTrkCutAlphaMin), 
  fTrkCutAlphaMax(rp.fTrkCutAlphaMax), 
  fTrkCutAngle(rp.fTrkCutAngle),
  fTrkCutNITS(rp.fTrkCutNITS),
  fTrkCutNTPC(rp.fTrkCutNTPC), // track matching
  fHighLowGainFactor(rp.fHighLowGainFactor), 
  fOrderParameter(rp.fOrderParameter), 
  fTau(rp.fTau), 
  fNoiseThreshold(rp.fNoiseThreshold), 
  fNPedSamples(rp.fNPedSamples) //raw signal
{
  //copy constructor

  //PID values
  Int_t i, j;
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      fGamma[i][j] = rp.fGamma[i][j];
      fHadron[i][j] = rp.fHadron[i][j];
      fPiZero5to10[i][j] = rp.fPiZero5to10[i][j];
      fPiZero10to60[i][j] = rp.fPiZero10to60[i][j];
    }
  }

}

//-----------------------------------------------------------------------------
AliEMCALRecParam& AliEMCALRecParam::operator = (const AliEMCALRecParam& rp)
{
  //assignment operator

  if(this != &rp) {
    fClusteringThreshold = rp.fClusteringThreshold;
    fW0 = rp.fW0;
    fMinECut = rp.fMinECut;
    fUnfold = rp.fUnfold;
    fLocMaxCut = rp.fLocMaxCut; //clustering
    fTrkCutX = rp.fTrkCutX;
    fTrkCutY = rp.fTrkCutY;
    fTrkCutZ = rp.fTrkCutZ;
    fTrkCutR = rp.fTrkCutR;
    fTrkCutAlphaMin = rp.fTrkCutAlphaMin;
    fTrkCutAlphaMax = rp.fTrkCutAlphaMax;
    fTrkCutAngle = rp.fTrkCutAngle;
    fTrkCutNITS = rp.fTrkCutNITS;
    fTrkCutNTPC = rp.fTrkCutNTPC; //track matching
    fHighLowGainFactor = rp.fHighLowGainFactor; 
    fOrderParameter = rp.fOrderParameter;
    fTau = rp.fTau;
    fNoiseThreshold = rp.fNoiseThreshold;
    fNPedSamples = rp.fNPedSamples; //raw signal

    //PID values
    Int_t i, j;
    for (i = 0; i < 6; i++) {
      for (j = 0; j < 6; j++) {
	fGamma[i][j] = rp.fGamma[i][j];
	fHadron[i][j] = rp.fHadron[i][j];
	fPiZero5to10[i][j] = rp.fPiZero5to10[i][j];
	fPiZero10to60[i][j] = rp.fPiZero10to60[i][j];
      }
    }

  }    
  
  return *this;

}

//-----------------------------------------------------------------------------
AliEMCALRecParam* AliEMCALRecParam::GetDefaultParameters()
{
  //default parameters for the reconstruction
  AliEMCALRecParam* params = new AliEMCALRecParam();
  params->SetName("Default - Pb+Pb");
  params->SetTitle("Default - Pb+Pb");
  return params;

}


//-----------------------------------------------------------------------------
AliEMCALRecParam* AliEMCALRecParam::GetLowFluxParam()
{
  //low flux/multiplicity ("p+p") parameters for the reconstruction
  AliEMCALRecParam* params = new AliEMCALRecParam();
  params->SetClusteringThreshold(0.2); // 200 MeV                                             
  params->SetMinECut(0.01);  //10 MeV       
  params->SetName("Low Flux - p+p");
  params->SetTitle("Low Flux - p+p");
  params->SetEventSpecie(AliRecoParam::kLowMult);

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

  return params;

}

//-----------------------------------------------------------------------------
void AliEMCALRecParam::Print(Option_t *) const
{
  // Print reconstruction parameters to stdout
  AliInfo(Form("Clusterization parameters :\n fClusteringThreshold=%.3f,\n fW0=%.3f,\n fMinECut=%.3f,\n fUnfold=%d,\n fLocMaxCut=%.3f \n",
	       fClusteringThreshold,fW0,fMinECut,fUnfold,fLocMaxCut));

  AliInfo(Form("Track-matching cuts :\n x %f, y %f, z %f, R %f \n alphaMin %f, alphaMax %f, Angle %f, NITS %f, NTPC %f\n", fTrkCutX, fTrkCutY, fTrkCutZ, fTrkCutR,fTrkCutAlphaMin,fTrkCutAlphaMax, fTrkCutAngle,fTrkCutNITS,fTrkCutNTPC));

  AliInfo(Form("PID parameters, Gamma :\n"));
  for(Int_t i = 0; i < 6; i++){
    for(Int_t j = 0; j < 6; j++){
      printf(" %f, ", fGamma[i][j]);
    }
    printf("\n");
  }

  printf("\n");

  AliInfo(Form("PID parameters, Hadron :\n"));
  for(Int_t i = 0; i < 6; i++){
    for(Int_t j = 0; j < 6; j++){
      printf(" %f, ", fHadron[i][j]);
    }
    printf("\n");
  }

  printf("\n");

  AliInfo(Form("PID parameters, Pi0zero5to10 :\n"));
  for(Int_t i = 0; i < 6; i++){
    for(Int_t j = 0; j < 6; j++){
      printf(" %f, ", fPiZero5to10[i][j]);
    }
    printf("\n");
  }

  printf("\n");

  AliInfo(Form("PID parameters, Pi0zero10to60 :\n"));
  for(Int_t i = 0; i < 6; i++){
    for(Int_t j = 0; j < 6; j++){
      printf(" %f, ", fPiZero10to60[i][j]);
    }
    printf("\n");
  }

  printf("\n");

  AliInfo(Form("Raw signal parameters: \n gain factor=%f, order=%d, tau=%f, noise threshold=%d, nped samples=%d \n",
	       fHighLowGainFactor,fOrderParameter,fTau,fNoiseThreshold,fNPedSamples));

}

//-----------------------------------------------------------------------------
const TObjArray* AliEMCALRecParam::GetMappings()
{
  //Returns array of AliAltroMappings for RCU0..RCUX.
  //If not found, read it from OCDB.                                           

  //Quick check as follows:                                                   
  //  root [0] AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB"
  //  root [1] AliCDBManager::Instance()->SetRun(1);                             
  //  root [2] TObjArray* maps = AliEMCALRecParam::GetMappings();                
  //  root [3] maps->Print();                                                    

  if(fgkMaps) return fgkMaps;

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("EMCAL/Calib/Mapping");
  if(entry)
    fgkMaps = (TObjArray*)entry->GetObject();

  return fgkMaps;

}

