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
static Double_t par_look_up[HCPARAMETERS][HCPARAMETERSETS] = 
{  
    {-1.68945e-04 , -1.68945e-04},   
    { 2.09684e-02 ,  2.09684e-02},   
    {-1.45683e-04 , -1.45683e-04},   
    { 6.64803e-03 ,  6.64803e-03},   
    { 2.00834e-02 ,  2.00834e-02},   
    { 9.22074e-07 ,  9.22074e-07}
};



ClassImp(AliEMCALHadronCorrectionv1)

AliEMCALHadronCorrectionv1* AliEMCALHadronCorrectionv1::fHadrCorr = 0;

void AliEMCALHadronCorrectionv1::SetGeometry(AliEMCALGeometry *geometry)
{
    if (!geometry)
    {
	SetParameters();
    }else
    {
        SetParameters(geometry->GetName());
    }
    return;	
}	
	
void AliEMCALHadronCorrectionv1::SetGeometry(TString name)
{
  if ( name == ""              ||
       name == "EMCAL_5655_21" ||
       name == "EMCALArch1a"   ||
       name == "EMCALArch1aN"  ||
       name == "G56_2_55_19"   ||
       name == "G56_2_55_19_104_14" )
  { // set parameters to this hadron correction
     for (Int_t i=0;i<6;i++)
     {
	   fPar[i] = par_look_up[i][0];  
     }
  }else if( name == "EMCAL_6564_21" ||  
	    name == "G65_2_64_19" )
  {	  
     for (Int_t i=0;i<6;i++)
     {
	   fPar[i] = par_look_up[i][1];  
     }
  }else
  {
    printf("Geometry not defined in hadron correction\n"); 
  }	  
	
}	

	
AliEMCALHadronCorrectionv1::AliEMCALHadronCorrectionv1(const char *name,const char *title) 
                           :AliEMCALHadronCorrection(name, title)
{
  fHadrCorr = this;
}

/*
AliEMCALHadronCorrectionv1::AliEMCALHadronCorrectionv1(const char *name,const char *title,AliEMCALGeometry *geometry)
{

  fHadrCorr = this;
  SetGeometry(geometry);  
	
}	
*/

AliEMCALHadronCorrectionv1*
AliEMCALHadronCorrectionv1::Instance()
{
  if (! fHadrCorr) new AliEMCALHadronCorrectionv1();
  return fHadrCorr;
}

Double_t 
AliEMCALHadronCorrectionv1::GetEnergy(const Double_t pmom,const Double_t eta,const Int_t gid)
{

  Double_t value =  fPar[5]*pmom*pmom*pmom+ fPar[0]*pmom*pmom+fPar[1]*pmom +fPar[2]*pmom*eta +fPar[3]*eta + fPar[4];
  return value;
   
}
