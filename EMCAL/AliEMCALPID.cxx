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
 * Revision 1.16  2007/11/23 13:39:05  gustavo
 * Track matching and PID parameters added to AliEMCALRecParam
 *
 * Revision 1.15  2007/10/09 08:46:10  hristov
 * The data members fEMCALClusterCluster and fPHOSCluster are removed from AliESDCaloCluster, the fClusterType is used to select PHOS or EMCAL clusters. Changes, needed to use correctly the new AliESDCaloCluster. (Christian)
 *
 * Revision 1.14  2007/07/26 16:54:53  morsch
 * Changes in AliESDEvent fwd declarartions.
 *
 * Revision 1.13  2007/07/11 13:43:29  hristov
 * New class AliESDEvent, backward compatibility with the old AliESD (Christian)
 *
 * Revision 1.12  2007/06/11 20:43:06  hristov
 * Changes required by the updated AliESDCaloCluster (Gustavo)
 *
 * Revision 1.11  2007/03/30 13:50:34  gustavo
 * PID for particles with E < 5 GeV was not done, temporal solution found (Guenole)
 *
 * Revision 1.10  2007/03/09 14:34:11  gustavo
 * Correct probability calculation, added missing initialization of data members
 *
 * Revision 1.9  2007/02/20 20:17:43  hristov
 * Corrected array size, removed warnings (icc)
 *
 * Revision 1.8  2006/12/19 08:49:35  gustavo
 * New PID class for EMCAL, bayesian analysis done with ESD data, PID information filled when calling AliEMCALPID in AliEMCALReconstructor::FillESD()
 *
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

// STEER includes
#include "AliLog.h"
#include "AliEMCALPID.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALRecParam.h"
#include "AliEMCALReconstructor.h"

  
ClassImp(AliEMCALPID)

//______________________________________________
  AliEMCALPID::AliEMCALPID():
    fPrintInfo(kFALSE), fProbGamma(0.),fProbPiZero(0.),fProbHadron(0.),fReconstructor(kFALSE)
{
  //
  // Constructor.
  // Initialize all constant values which have to be used
  // during PID algorithm execution
  //
 
  fPIDWeight[0] = -1;
  fPIDWeight[1] = -1;
  fPIDWeight[2] = -1;

  for(Int_t i=0; i<AliPID::kSPECIESN+1; i++)
    fPIDFinal[i]= 0;

  const AliEMCALRecParam* recParam = AliEMCALReconstructor::GetRecParam();
  if(!recParam) {
    AliFatal("Reconstruction parameters for EMCAL not set!");
  }
  else {
    for(Int_t i=0; i<6; i++){
      for(Int_t j=0; j<6; j++){
	fGamma[i][j] = recParam->GetGamma(i,j);
	fHadron[i][j] = recParam->GetHadron(i,j);
	fPiZero5to10[i][j] = recParam->GetPiZero5to10(i,j);
	fPiZero10to60[i][j] = recParam->GetPiZero10to60(i,j);
	AliDebug(1,Form("PID parameters (%d, %d): fGamma=%.3f, fPi=%.3f, fHadron=%.3f",
			i,j, fGamma[i][j],fPiZero5to10[i][j],fHadron[i][j] ));
      }
    }
    
  }
  
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

      if(fReconstructor) // In case it is called during reconstruction.
	clust->SetPid(fPIDFinal);
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

if (energy<5){energy =6;}


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
  fPIDFinal[6]  = fPIDWeight[1]  ;
  fPIDFinal[7]  = fPIDWeight[2]/8;
  fPIDFinal[8]  = fPIDWeight[2]/8;
  fPIDFinal[9]  = fPIDWeight[2]/8;
  fPIDFinal[10] = fPIDWeight[2]/8;
}

//________________________________________________________
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

//_______________________________________________________
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
