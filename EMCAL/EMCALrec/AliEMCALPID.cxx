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

// STEER includes
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"

// EMCAL includes
#include "AliEMCALReconstructor.h"
#include "AliEMCALPID.h"

/// \cond CLASSIMP
ClassImp(AliEMCALPID) ;
/// \endcond
 
///
/// Constructor.
/// Initialize all constant values which have to be used
/// during PID algorithm execution.
///
//______________________________________________
AliEMCALPID::AliEMCALPID()
: AliEMCALPIDUtils(), fReconstructor(kTRUE)
{  
  InitParameters(); 
}

///
/// Constructor.
/// Initialize all constant values which have to be used
/// during PID algorithm execution called when used in standalone mode. 
///
/// \param reconstructor: bool with reconstruction/analysis mode
///
//______________________________________________
AliEMCALPID::AliEMCALPID(Bool_t reconstructor)
: AliEMCALPIDUtils(), fReconstructor(reconstructor)
{  
  InitParameters(); 
}

///
/// Make the PID for all the EMCAL clusters containedin the ESDs File
/// but just gamma/PiO/Hadron.
///
/// Generalize to VEvents? Move to AliEMCALPIDUtils?
///
/// \param esd: pointer to ESD event 
///
//______________________________________________
void AliEMCALPID::RunPID(AliESDEvent *esd)
{
  // Trivial check against NULL object passed.
  if (esd == 0x0) 
  {
    AliInfo("NULL ESD object passed !!" );
    return ;
  }
  
  Int_t nClusters = esd->GetNumberOfCaloClusters();
  Int_t firstCluster = 0;
  Double_t energy=0., lambda0=0.;
  
  for (Int_t iCluster = firstCluster; iCluster < (nClusters + firstCluster); iCluster++) 
  {
    AliESDCaloCluster *clust = esd->GetCaloCluster(iCluster);
    if (!clust->IsEMCAL()) continue ; 
    
    energy = clust->E();
    lambda0 = clust->GetM02();
    
    if (lambda0 != 0  && energy < 1000) 
    {
      // reject clusters with lambda0 = 0
      
      ComputePID(energy, lambda0);
      
      if ( fPrintInfo ) 
      {
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
      
      if ( fReconstructor )
      { // In case it is called during reconstruction.
        //	cout << "############# Fill ESDs with PIDWeight ##########" << endl;
        clust->SetPID(fPIDFinal);
      }
      
    } // end if (lambda0 != 0  && energy < 1000)
  } // end for (iCluster...)
}

///
/// Initialize PID parameters, depending on the use or not of the reconstructor
/// and the kind of event type if the reconstructor is not used.
///
//______________________________________________
void AliEMCALPID::InitParameters()
{
  //  fWeightHadronEnergy=0.;
  //  fWeightPiZeroEnergy=0.;
  //  fWeightGammaEnergy=0.;
  
  fPIDWeight[0] = -1;
  fPIDWeight[1] = -1;
  fPIDWeight[2] = -1;
  
  for(Int_t i=0; i<AliPID::kSPECIESCN+1; i++)
    fPIDFinal[i]= 0;
  
  const AliEMCALRecParam* recParam = AliEMCALReconstructor::GetRecParam();
  
  if(fReconstructor)
  {
    if(!recParam) 
    {
      AliFatal("Reconstruction parameters for EMCAL not set!");
    }
    else 
    {
      for(Int_t i=0; i<6; i++)
      {
        for(Int_t j=0; j<6; j++)
        {
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
  
  else
  {
    //   init the parameters here instead of from loading from recparam
    //   default parameters are PbPb parameters.
    SetHighFluxParam();
  }
  
}

