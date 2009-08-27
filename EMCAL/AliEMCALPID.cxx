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
// Author: Genole Bourdaud 2007 (SUBATECH)
//         Marie Germain 07/2009 (SUBATECH), new parametrization for low and high flux environment
//         Gustavo Conesa 08/2009 (LNF), divide class in AliEMCALPID and AliEMCALPIDUtils, PIDUtils belong to library EMCALUtils 
// --- standard c ---

// standard C++ includes
//#include <Riostream.h>

// ROOT includes

// STEER includes
#include "AliESDEvent.h"
#include "AliEMCALPID.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALReconstructor.h"

  
ClassImp(AliEMCALPID)
  
//______________________________________________
  AliEMCALPID::AliEMCALPID()
	: AliEMCALPIDUtils(), fReconstructor(kTRUE)
{
  //
  // Constructor.
  // Initialize all constant values which have to be used
  // during PID algorithm execution
  //
  
  InitParameters(); 
  
  
}

//______________________________________________
AliEMCALPID::AliEMCALPID(Bool_t reconstructor)
: AliEMCALPIDUtils(), fReconstructor(reconstructor)
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
	//	cout << "############# Fill ESDs with PIDWeight ##########" << endl;
	clust->SetPid(fPIDFinal);}
    } // end if (clusterType...)
  } // end for (iCluster...)
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

