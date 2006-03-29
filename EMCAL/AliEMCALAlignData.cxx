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
#include "AliAlignObjMatrix.h"

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
    fSuperModuleMatrix[module] = alignda.fSuperModuleMatrix[module];
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
    fSuperModuleMatrix[module] = new AliAlignObjMatrix(*alignda.fSuperModuleMatrix[module]);
  }
  return *this;
}

//________________________________________________________________
AliEMCALAlignData::~AliEMCALAlignData()
{
  // Destructor
  for(Int_t module=0; module<fNSuperModules; module++) {
    if(fSuperModuleMatrix[module]) delete fSuperModuleMatrix[module];
  }
}

//________________________________________________________________
void AliEMCALAlignData::Reset()
{
  // Set all to default values
  fNSuperModules = 12;
  memset(fSuperModuleMatrix,0,12*sizeof(AliAlignObjMatrix*));
  for(Int_t module=0; module<fNSuperModules; module++) fSuperModuleMatrix[module] = 0;
}

//________________________________________________________________
void  AliEMCALAlignData::Print(Option_t */*option =""*/) const
{
  // Print alignment data

  printf("EMCAL alignment object\n");
  printf("     Number of modules: %d\n",fNSuperModules);
}
//________________________________________________________________
