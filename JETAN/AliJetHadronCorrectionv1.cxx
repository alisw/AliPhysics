/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
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

//===============================================================
// To be modified for hadron correction using particle ID and track merging
// Author : magali.estienne@subatech.in2p3.fr
//===============================================================
// Author : Mark Horner (LBL/UCT)

// --- Standard library ---
#include "Riostream.h"
#include "TMath.h"

// --- AliRoot header files ---
#include "AliAODTrack.h"
#include "AliJetDummyGeo.h"
#include "AliJetHadronCorrectionv1.h"

static Double_t etaGrid[HCPARAMETERS]={ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.67};

ClassImp(AliJetHadronCorrectionv1)

Double_t AliJetHadronCorrectionv1::fgParLookup[HCPARAMETERS][HCPARAMETERSETS] = 
{  

    {7.52624e-2  , 9.80119e-2   , 1.29086e-2},
    {-2.07449e-2 , 2.22951e-2   , 8.31971e-2},
    {-3.64861e-4 , -3.65680e-03 , -3.07148e-02},
    {5.35252e-3  , 5.12234e-03  , 2.25559e-02},
    {1.27106e-1  , 1.25248e-01  , 1.06352e-01},
    {-6.72909e-2 , 9.78677e-01  , -1.35909e-01},
    {-3.80660e-5 , -2.77522e-05 , -2.05218e-05},
    {1.54855e-7  , 1.00604e-07  , 1.02481e-07}
};

AliJetHadronCorrectionv1* AliJetHadronCorrectionv1::fgHadrCorr = 0;

AliJetHadronCorrectionv1::AliJetHadronCorrectionv1(const char *name,const char *title) 
    :AliJetHadronCorrection(name, title), 
     fSamplingFraction(0)
{
  fgHadrCorr = this;
  for (Int_t i = 0; i < 8; i++) fPar[i] = 0.;
}

AliJetHadronCorrectionv1*
AliJetHadronCorrectionv1::Instance()
{
  // return pointer to global instance. Instantiate if needed
  if (! fgHadrCorr) fgHadrCorr = new AliJetHadronCorrectionv1();
  return fgHadrCorr;
}

void AliJetHadronCorrectionv1::SetGeometry2(const AliJetDummyGeo *geometry)
{
  // Initialise EMCAL geometry
    if (!geometry)
    {
	SetParameters();
//      commented out since geometry is 0x0 after null check
//	fSamplingFraction = geometry->GetSampling();
    }else
    {
        SetParameters(geometry->GetName());
	fSamplingFraction = geometry->GetSampling();
    }
    return;	
}	

/*
###########################################################################
###########################################################################

This will have to be modified with the inclusion of the new EMCAL geometry
That means I guess a new study of the hadron correction...

Geometry : SHISH...

###########################################################################
###########################################################################
*/

void AliJetHadronCorrectionv1::SetGeometry(TString /*name*/,Double_t fs)
{
  // Initialise EMCAL geometry
  fSamplingFraction = fs;	

  /*
  if ( name == ""              ||
       name == "EMCAL_5655_21" ||
       name == "EMCALArch1a"   ||
       name == "EMCALArch1aN"  ||
       name == "G56_2_55_19"   ||
       name == "G56_2_55_19_104_14" )
    { // set parameters to this hadron correction
      cout<<"HC parameters!"<<endl;
      for (Int_t i=0;i<6;i++)
	{
	  fPar[i] = fgParLookup[i][0];  
	  cout <<fPar[i]<<endl;
	}
    }else if( name == "EMCAL_6564_21" ||  
	      name == "G65_2_64_19" )
    {	  
      cout<<"HC parameters!"<<endl;
      for (Int_t i=0;i<6;i++)
	{
	  fPar[i] = fgParLookup[i][1];  
	  cout <<fPar[i]<<endl;
	}
    }else
    {
      printf("Geometry not defined in hadron correction\n"); 
    }	  
  */

}	

// Double_t AliJetHadronCorrectionv1::GetEnergy(Double_t pmom, Double_t eta, Int_t /*gid*/)
Double_t AliJetHadronCorrectionv1::GetEnergy(Double_t pmom, Double_t eta, Int_t /*gid*/)
{
  // Return parametrised energy response
  Double_t etai = TMath::Abs(eta); 

  if(etai < etaGrid[1]) {
    for (Int_t i=0;i<8;i++)
      {
	fPar[i] = fgParLookup[i][1];  
	//	cout << "fPar[" << i << "]: " << fPar[i] << endl;
      }
  } else if(etai >= etaGrid[1] && etai <= etaGrid[HCPARAMETERS-2]) {
    for (Int_t i=0;i<8;i++)
      {
	fPar[i] = fgParLookup[i][0];  
	//	cout << "fPar[" << i << "]: " << fPar[i] << endl;
      }
  } else {
    for (Int_t i=0;i<8;i++)
      {
	fPar[i] = fgParLookup[i][2];  
	//	cout << "fPar[" << i << "]: " << fPar[i] << endl;
      }

  }

  Double_t value = fPar[5]*pow(etai,3) + 
                   fPar[0]*pow(etai,2) + 
                   fPar[1]*etai        +
                   fPar[2]*etai*pmom   +
                   fPar[3]*pmom        + 
                   fPar[4]             +
                   fPar[6]*pow(pmom,2) +
                   fPar[7]*pow(pmom,3);

  return fSamplingFraction*value;
   
}

void AliJetHadronCorrectionv1::TrackPositionEMCal(const AliAODTrack* track,Double_t &eta, Double_t &phi)
{
// Return track position on EMCal
  AliAODPid*    pid   = (AliAODPid*) track->GetDetPid();
    
  if(pid) {
    Double_t emcpos[3];
    pid->GetEMCALPosition(emcpos);      
    TVector3 tpos(emcpos[0],emcpos[1],emcpos[2]);
    
    eta = tpos.Eta();
    phi = tpos.Phi();

  }

}
