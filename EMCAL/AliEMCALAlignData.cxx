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

////////////////////////////////////////////////
//  class for EMCAL alignment parameters       //
////////////////////////////////////////////////

#include "AliEMCALAlignData.h"

ClassImp(AliEMCALAlignData)

//________________________________________________________________
AliEMCALAlignData::AliEMCALAlignData()
{
  // Default constructor
  Reset();
}

//________________________________________________________________
AliEMCALAlignData::AliEMCALAlignData(const char* name)
{
  // Constructor
  TString namst = "Align_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
}

//________________________________________________________________
AliEMCALAlignData::AliEMCALAlignData(const AliEMCALAlignData& alignda) :
  TNamed(alignda)
{
  // copy constructor
  SetName(alignda.GetName());
  SetTitle(alignda.GetName());
  Reset();
  fNSuperModules = alignda.GetNSuperModules();
  for(Int_t module=0; module<fNSuperModules; module++) {
    for (Int_t axis=0; axis<3; axis++) {
      fSuperModuleCenter[module][axis] = 
	alignda.GetSuperModuleCenter(module,axis);
      for (Int_t angle=0; angle<2; angle++) {
	fSuperModuleAngle[module][axis][angle] = 
	  alignda.GetSuperModuleAngle(module,axis,angle);
      }
    }
  }
}

//________________________________________________________________
AliEMCALAlignData &AliEMCALAlignData::operator =(const AliEMCALAlignData& alignda)
{
  // assignment operator
  SetName(alignda.GetName());
  SetTitle(alignda.GetName());
  Reset();
  fNSuperModules = alignda.GetNSuperModules();
  for(Int_t module=0; module<fNSuperModules; module++) {
    for (Int_t axis=0; axis<3; axis++) {
      fSuperModuleCenter[module][axis] = 
	alignda.GetSuperModuleCenter(module,axis);
      for (Int_t angle=0; angle<2; angle++) {
	fSuperModuleAngle[module][axis][angle] = 
	  alignda.GetSuperModuleAngle(module,axis,angle);
      }
    }
  }
  return *this;
}

//________________________________________________________________
AliEMCALAlignData::~AliEMCALAlignData()
{
  // Destructor
}

//________________________________________________________________
void AliEMCALAlignData::Reset()
{
  // Set all to default values
  fNSuperModules = 12;
  memset(fSuperModuleCenter,0,12*3*sizeof(Float_t));
  memset(fSuperModuleAngle ,0,12*3*2*sizeof(Float_t));
}

//________________________________________________________________
void  AliEMCALAlignData::Print(Option_t */*option =""*/) const
{
  // Print alignment data

  printf("EMCAL alignment object\n");
  printf("     Number of modules: %d\n",fNSuperModules);
}

//________________________________________________________________
