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
//  class for PHOS alignment parameters       //
////////////////////////////////////////////////

#include "AliPHOSAlignData.h"

ClassImp(AliPHOSAlignData)

//________________________________________________________________
AliPHOSAlignData::AliPHOSAlignData()
{
  // Default constructor
  Reset();
}

//________________________________________________________________
AliPHOSAlignData::AliPHOSAlignData(const char* name)
{
  // Constructor
  TString namst = "Align_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
}

//________________________________________________________________
AliPHOSAlignData::AliPHOSAlignData(const AliPHOSAlignData& alignda) :
  TNamed(alignda)
{
  // copy constructor
  SetName(alignda.GetName());
  SetTitle(alignda.GetName());
  Reset();
  fNModules = alignda.GetNModules();
  for(Int_t module=0; module<fNModules; module++) {
    for (Int_t axis=0; axis<3; axis++) {
      fModuleCenter[module][axis] = 
	alignda.GetModuleCenter(module,axis);
      for (Int_t angle=0; angle<2; angle++) {
	fModuleAngle[module][axis][angle] = 
	  alignda.GetModuleAngle(module,axis,angle);
      }
    }
  }
}

//________________________________________________________________
AliPHOSAlignData &AliPHOSAlignData::operator =(const AliPHOSAlignData& alignda)
{
  // assignment operator
  SetName(alignda.GetName());
  SetTitle(alignda.GetName());
  Reset();
  fNModules = alignda.GetNModules();
  for(Int_t module=0; module<fNModules; module++) {
    for (Int_t axis=0; axis<3; axis++) {
      fModuleCenter[module][axis] = 
	alignda.GetModuleCenter(module,axis);
      for (Int_t angle=0; angle<2; angle++) {
	fModuleAngle[module][axis][angle] = 
	  alignda.GetModuleAngle(module,axis,angle);
      }
    }
  }
  return *this;
}

//________________________________________________________________
AliPHOSAlignData::~AliPHOSAlignData()
{
  // Destructor
}

//________________________________________________________________
void AliPHOSAlignData::Reset()
{
  // Set all to default values
  fNModules = 5;
  memset(fModuleCenter,0,5*3*sizeof(Float_t));
  memset(fModuleAngle ,0,5*3*2*sizeof(Float_t));
}

//________________________________________________________________
void  AliPHOSAlignData::Print(Option_t */*option =""*/) const
{
  // Print alignment data

  printf("PHOS alignment object\n");
  printf("     Number of modules: %d\n",fNModules);
  printf("     Module centers in MARS:\n");
  Int_t iModule;
  for (iModule=0; iModule<fNModules; iModule++) {
    printf("     (x=%.3f, y=%.3f, z=%.3f cm\n",
	   fModuleCenter[iModule][0],
	   fModuleCenter[iModule][1],
	   fModuleCenter[iModule][2]);
  }
  printf("     Module orientation angles:\n");
  for (iModule=0; iModule<fNModules; iModule++) {
    printf("     (theta1=%.3f, phi1=%.3f degrees\n",
	   fModuleAngle[iModule][0][0],
	   fModuleAngle[iModule][0][1]);
    printf("     (theta2=%.3f, phi2=%.3f degrees\n",
	   fModuleAngle[iModule][1][0],
	   fModuleAngle[iModule][1][1]);
    printf("     (theta3=%.3f, phi3=%.3f degrees\n",
	   fModuleAngle[iModule][2][0],
	   fModuleAngle[iModule][2][1]);
  }
}

//________________________________________________________________
