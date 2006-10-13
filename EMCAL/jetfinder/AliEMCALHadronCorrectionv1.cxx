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


// Author : Mark Horner (LBL/UCT)

#include "AliEMCALHadronCorrectionv1.h"
#include "AliEMCALGeometry.h"
#include "Riostream.h"
#include "TMath.h"

ClassImp(AliEMCALHadronCorrectionv1)

Double_t AliEMCALHadronCorrectionv1::fgParLookup[HCPARAMETERS][HCPARAMETERSETS] = 
{  
    {-2.82271e-4 , -2.39954e-4},   
    {2.50796e-2  , 2.07172e-2},   
    {1.02861e-3  , 1.48576e-3},   
    {2.11539e-2  , -1.38473-2},   
    {2.27003e-2  , 2.78252e-2},   
    {1.65078e-6  , 1.51821e-6}
};

AliEMCALHadronCorrectionv1* AliEMCALHadronCorrectionv1::fgHadrCorr = 0;

//______________________________________________________________________
void AliEMCALHadronCorrectionv1::SetGeometry(AliEMCALGeometry *geometry)
{
  // Initialise EMCAL geometry
    if (!geometry)
    {
	SetParameters();
	fSamplingFraction = geometry->GetSampling();
	cout<<"Setting the sampling fraction to :"<<fSamplingFraction<<endl; 
    }else
    {
        SetParameters(geometry->GetName());
	fSamplingFraction = geometry->GetSampling();
	cout<<"Setting the sampling fraction to :"<<fSamplingFraction<<endl; 
    }
    return;	
}	
	
//______________________________________________________________________
void AliEMCALHadronCorrectionv1::SetGeometry(TString name,Double_t fs)
{
  // Initialise EMCAL geometry
  cout << "Setting sampling fraction to "<<fSamplingFraction<<endl;	
   fSamplingFraction = fs;	
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
	
}	

//______________________________________________________________________	
AliEMCALHadronCorrectionv1::AliEMCALHadronCorrectionv1(const char *name,const char *title) 
  :AliEMCALHadronCorrection(name, title),
   fSamplingFraction(0.)
{
  fgHadrCorr = this;
  fSamplingFraction = 1.0;
  for(Int_t i = 0; i < 6; i++) fPar[i] = 0.;

}

/*
AliEMCALHadronCorrectionv1::AliEMCALHadronCorrectionv1(const char *name,const char *title,AliEMCALGeometry *geometry)
{

  fHadrCorr = this;
  SetGeometry(geometry);  
	
}	
*/

//______________________________________________________________________
AliEMCALHadronCorrectionv1* AliEMCALHadronCorrectionv1::Instance()
{
  // return pointer to global instance. Instantiate if needed
  if (! fgHadrCorr) new AliEMCALHadronCorrectionv1();
  return fgHadrCorr;
}

//______________________________________________________________________
Double_t AliEMCALHadronCorrectionv1::GetEnergy(Double_t pmom, Double_t eta, Int_t /*gid*/)
{
  // Return parametrised energy response
  Double_t etai = TMath::Abs(eta); 
  Double_t value =  fPar[5]*pmom*pmom*pmom+ fPar[0]*pmom*pmom+fPar[1]*pmom +fPar[2]*pmom*etai +fPar[3]*etai + fPar[4];
  return fSamplingFraction*value;
   
}
