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
/* History of cvs commits:
 *
 * $Log$
 *
 */
//    to compute PID for all the clusters in ESDs.root file
//     the ESDs.root have to be in the same directory as the class
//
//
//
//
//
//	AliEMCALPID::CalculPID(Energy,Lambda0)
//	  Calcul PID for all clusters in AliESDs.root file
//	  keep this function for the moment for a simple verification, could be removed
//
//
//
//   AliEMCALPID::CalculPID(Energy,Lambda0)
//    calcul PID Weght for a cluster with Energy, Lambda0 .
//    Double_t PIDFinal[AliPID::kSPECIESN]  is the standard PID for :
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
//  
//
//
//
// --- ROOT system ---

// standard C++ includes
#include <Riostream.h>
// #include <cstdlib>
// #include <cmath>

// ROOT includes
#include "TTree.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TParticle.h"

// // STEER includes
// #include "AliRun.h"
// #include "AliRunLoader.h"
// #include "AliHeader.h"
// #include "AliLoader.h"
// #include "AliStack.h"
// #include "AliESDtrack.h"
// #include "AliESD.h"
#include "AliLog.h"
#include "AliEMCALPID.h"
  
ClassImp(AliEMCALPID)

AliEMCALPID::AliEMCALPID()
{
//
// Constructor.
// Initialize all constant values which have to be used
// during PID algorithm execution
//

	// set flag for printing to FALSE by default
	fPrintInfo = kFALSE;
	
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
	
	fPIDWeight[0] = -1;
	fPIDWeight[1] = -1;
	fPIDWeight[2] = -1;
	fReconstructor = kFALSE;
}
//
//
void AliEMCALPID::RunPID(AliESD *esd)
{
//
// Make the PID for all the EMCAL clusters containedin the ESDs File
// but just gamma/PiO/Hadron
//
	// trivial check against NULL object passed

  if (esd == 0x0) {
    AliInfo("NULL ESD object passed!!" );
    return ;
  }
  Int_t nClusters = esd->GetNumberOfEMCALClusters();
  Int_t firstCluster = esd->GetFirstEMCALCluster();
  Double_t energy, lambda0;
  for (Int_t iCluster = firstCluster; iCluster < (nClusters + firstCluster); iCluster++) {
    
    AliESDCaloCluster *clust = esd->GetCaloCluster(iCluster);
    energy = clust->GetClusterEnergy();
    lambda0 = clust->GetM02();
    // verify cluster type
    Int_t clusterType= clust->GetClusterType();
    if (clusterType == AliESDCaloCluster::kClusterv1 && lambda0 != 0 && energy > 5 && energy < 1000) {
      // reject clusters with lambda0 = 0
      // reject clusters with energy < 5 GeV
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
      if(fReconstructor) // In case it is called during reconstruction.
	clust->SetPid(fPIDFinal);
    } // end if (clusterType...)
  } // end for (iCluster...)
}
//
//
void AliEMCALPID::ComputePID(Double_t energy, Double_t lambda0)
{
//
// This is the main command, which uses the distributions computed and parametrised, 
// and gives the PID by the bayesian method.
//
	TArrayD paramDistribGamma  = DistLambda0(energy, 1);
	TArrayD paramDistribPiZero = DistLambda0(energy, 2);
	TArrayD paramDistribHadron = DistLambda0(energy, 3);
		
	Bool_t norm = kFALSE;
	
	fProbGamma   = TMath::Gaus(lambda0, paramDistribGamma[1], paramDistribGamma[2], norm) * paramDistribGamma[0];
	fProbGamma  += TMath::Landau(lambda0, paramDistribGamma[4], paramDistribGamma[5], norm) * paramDistribGamma[3];
	fProbPiZero  = TMath::Gaus(lambda0, paramDistribPiZero[1], paramDistribPiZero[2], norm) * paramDistribPiZero[0];
	fProbPiZero += TMath::Landau(lambda0, paramDistribPiZero[4], paramDistribPiZero[5], norm) * paramDistribPiZero[3];
	fProbHadron  = TMath::Gaus(lambda0, paramDistribHadron[1], paramDistribHadron[2], norm) * paramDistribHadron[0];
	fProbHadron += TMath::Landau(lambda0, paramDistribHadron[4], paramDistribHadron[5], norm) * paramDistribHadron[3];
	
	// compute PID Weight
	fPIDWeight[0] = fProbGamma / (fProbGamma + fProbPiZero + fProbHadron);
	fPIDWeight[1] = fProbPiZero / (fProbGamma+fProbPiZero+fProbHadron);
	fPIDWeight[2] = fProbHadron / (fProbGamma+fProbPiZero+fProbHadron);
	
	SetPID(fPIDWeight[0], 0);
	SetPID(fPIDWeight[1], 1);
	SetPID(fPIDWeight[2], 2);
	
	// sortie ecran pid Weight only for control (= in english ???)
	if (fPrintInfo) {
	  AliInfo(Form( "Energy in loop = %f", energy) );
	  AliInfo(Form( "Lambda0 in loop = %f", lambda0) );
	  AliInfo(Form( "fProbGamma in loop = %f", fProbGamma) );
		//	AliInfo(Form( "fParametresDistribGamma[2] = %f", fParamDistribGamma[2]) );
	  AliInfo(Form( "fProbaPiZero = %f", fProbPiZero ));
	  AliInfo(Form( "fProbaHadron = %f", fProbHadron) );
	  AliInfo(Form( "PIDWeight in loop = %f ||| %f ||| %f",  fPIDWeight[0] , fPIDWeight[1], fPIDWeight[2]) );
	  AliInfo(Form( "fGamma[2][2] = %f", fGamma[2][2] ));
	  AliInfo("********************************************************" );
	}
	
	fPIDFinal[0]  = fPIDWeight[0]/2;
	fPIDFinal[1]  = fPIDWeight[2]/8;
	fPIDFinal[2]  = fPIDWeight[2]/8;
	fPIDFinal[3]  = fPIDWeight[2]/8;
	fPIDFinal[4]  = fPIDWeight[2]/8;
	fPIDFinal[5]  = fPIDWeight[0]/2;
	fPIDFinal[6]  = fPIDWeight[1]/2;
	fPIDFinal[7]  = fPIDWeight[2]/8;
	fPIDFinal[8]  = fPIDWeight[2]/8;
	fPIDFinal[9]  = fPIDWeight[2]/8;
	fPIDFinal[10] = fPIDWeight[2]/8;
	fPIDFinal[11] = 0;
}
//
//
TArrayD AliEMCALPID::DistLambda0(Double_t energy, Int_t type)
{
//
// Compute the values of the parametrised distributions using the data initialised before.
//
	Double_t constGauss = 0., meanGauss = 0., sigmaGauss = 0.;
	Double_t constLandau=0., mpvLandau=0., sigmaLandau=0.;
	TArrayD  distributionParam(6);
	
	switch (type) {
		case 1:
			constGauss  = Polynomial(energy, fGamma[0]);
			meanGauss   = Polynomial(energy, fGamma[1]);
			sigmaGauss  = Polynomial(energy, fGamma[2]);
			constLandau = Polynomial(energy, fGamma[3]);
			mpvLandau   = Polynomial(energy, fGamma[4]);
			sigmaLandau = Polynomial(energy, fGamma[5]);
			break;
		case 2:
			if (energy < 10) {
				constGauss  = Polynomial(energy, fPiZero5to10[0]);
				meanGauss   = Polynomial(energy, fPiZero5to10[1]);
				sigmaGauss  = Polynomial(energy, fPiZero5to10[2]);
				constLandau = Polynomial(energy, fPiZero5to10[3]);
				mpvLandau   = Polynomial(energy, fPiZero5to10[4]);
				sigmaLandau = Polynomial(energy, fPiZero5to10[5]);
			}
			else {
				constGauss  = Polynomial(energy, fPiZero10to60[0]);
				meanGauss   = Polynomial(energy, fPiZero10to60[1]);
				sigmaGauss  = Polynomial(energy, fPiZero10to60[2]);
				constLandau = Polynomial(energy, fPiZero10to60[3]);
				mpvLandau   = Polynomial(energy, fPiZero10to60[4]);
				sigmaLandau = Polynomial(energy, fPiZero10to60[5]);
			}
			break;
		case 3:
			constGauss  = Polynomial(energy, fHadron[0]);
			meanGauss   = Polynomial(energy, fHadron[1]);
			sigmaGauss  = Polynomial(energy, fHadron[2]);
			constLandau = Polynomial(energy, fHadron[3]);
			mpvLandau   = Polynomial(energy, fHadron[4]);
			sigmaLandau = Polynomial(energy, fHadron[5]);
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
//
//
Double_t AliEMCALPID::Polynomial(Double_t x, Double_t *params)
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
