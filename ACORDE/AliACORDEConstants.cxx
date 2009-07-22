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

////////////////////////////////////////////////////////////////////////
//
// AliACORDEConstants class
//
// This class serves to group constants needed by ACORDE detector in 1
// easily accessible place. All constants are public const static data 
// members. The class is never instatiated.
// Authors: Arturo Fernandez, Enrique Gamez, Mario Rodriguez Cahuantzi, Eleazar Cuautle(ICN-UNAM) 
//         FCFM-UAP, Mexico.
// Last update: Nov. 24th 08
////////////////////////////////////////////////////////////////////////

#include "AliACORDEConstants.h"

AliACORDEConstants* AliACORDEConstants::fgInstance = 0;

const Float_t AliACORDEConstants::fgkModuleLength          = 300.0;
const Float_t AliACORDEConstants::fgkModuleWidth           = 26.0;
const Float_t AliACORDEConstants::fgkModuleHeight          =  10.0;
const Float_t AliACORDEConstants::fgkPlasticLength = 190.0;
const Float_t AliACORDEConstants::fgkPlasticWidth  =  20.0;
const Float_t AliACORDEConstants::fgkPlasticHeight =   1.0;
const Float_t AliACORDEConstants::fgkProfileWidth =    3.8;
const Float_t AliACORDEConstants::fgkProfileThickness = 0.3;
const Float_t AliACORDEConstants::fgkDepth               =4420; 

const Float_t AliACORDEConstants::fgkHitEnergyThreshold = 1.52; // MeV
const Float_t AliACORDEConstants::fgkMaxHitTimeDifference = 40.0; // ns
const Int_t AliACORDEConstants::fgkMultiMuonThreshold = 2;
const Float_t AliACORDEConstants::fgkMultiMuonWindow = 25;
const Float_t AliACORDEConstants::fgkModulePositionX[60] = {
  641, 641, 641, 641, 641, 641, 641, 641, 641, 641,
  426, 426, 426, 426, 426, 426, 426, 426, 426, 426,
  153, 153, 153, 153, 153, 153, 153, 153, 153, 153,
  -153, -153, -153, -153, -153, -153, -153, -153, -153,
  -153, -426, -426, -426, -426, -426, -426, -426, -426,
  -426, -426, -644, -644, -644, -644, -644, -619, -623,
  -641, -641, -641};
const Float_t AliACORDEConstants::fgkModulePositionY[60] = {
  582, 574, 574, 574, 574, 574, 574, 574, 574, 582,
  789, 789, 789, 789, 789, 789, 789, 789, 789, 789,
  850, 850, 850, 850, 850, 850, 850, 850, 850, 850,
  850, 850, 850, 850, 850, 850, 850, 850, 850, 850,
  789, 789, 789, 789, 789, 789, 789, 789, 789, 789,
  582, 574, 574, 574, 574, 601, 597, 574, 574, 582};
const Float_t AliACORDEConstants::fgkModulePositionZ[60] = {
  450, 350, 250, 150, 50, -50, -120, -280, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 104, 50, -85, -184, -258, -350, -450};


const Float_t AliACORDEConstants::fgkSupportModulePositionX[60] = {
  641, 641, 641, 641, 641, 641, 641, 641, 641, 641,
  426, 426, 426, 426, 426, 426, 426, 426, 426, 426,
  153, 153, 153, 153, 153, 153, 153, 153, 153, 153,
  -153, -153, -153, -153, -153, -153, -153, -153, -153,
  -153, -426, -426, -426, -426, -426, -426, -426, -426,
  -426, -426, -644, -644, -644, -644, -644, -619, -623,
  -641, -641, -641};
const Float_t AliACORDEConstants::fgkSupportModulePositionY[60] = {
  582, 582, 582, 582, 582, 582, 582, 582, 582, 582,
  797, 797, 797, 797, 797, 797, 797, 797, 797, 797,
  850, 850, 850, 850, 850, 850, 850, 850, 850, 850,
  850, 850, 850, 850, 850, 850, 850, 850, 850, 850,
  797, 797, 797, 797, 797, 797, 797, 797, 797, 797,
  582, 582, 582, 582, 582, 609, 605, 582, 582, 582};
const Float_t AliACORDEConstants::fgkSupportModulePositionZ[60] = {
  450, 350, 250, 150, 50, -50, -120, -280, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 104, 50, -85, -176, -250, -350, -450};
  

const Float_t AliACORDEConstants::fgkExtraModulePositionZ[4] = {93.0, 18., -18, -93};
const Float_t AliACORDEConstants::fgkExtraModulePositionX = 0.0;
const Float_t AliACORDEConstants::fgkExtraModulePositionY = 850.0;
const Int_t AliACORDEConstants::fgkModuleElectronicChannel[60] = {
/* DCS 0_0 ITS-1*/ 10,
/* DCS 0_1 */ 4,
/* DCS 0_2 */ 8,
/* DCS 0_3 */ 7,
/* DCS 0_4 */ 6,
/* DCS 0_5 */ 5,
/* DCS 0_6 */ 9,
/* DCS 0_7 */ 3,
/* DCS 0_8 */ 2,
/* DCS 0_9 ITS-2*/ 42,
/* DCS 1_0 */ 20,
/* DCS 1_1 */ 19,
/* DCS 1_2 */ 18,
/* DCS 1_3 */ 17,
/* DCS 1_4 */ 16,
/* DCS 1_5 */ 15,
/* DCS 1_6 */ 14,
/* DCS 1_7 */ 13,
/* DCS 1_8 */ 12,
/* DCS 1_9 */ 11,
/* DCS 2_0 */ 60,
/* DCS 2_1 */ 59,
/* DCS 2_2 */ 58,
/* DCS 2_3 */ 57,
/* DCS 2_4 */ 56,
/* DCS 2_5 */ 55,
/* DCS 2_6 */ 54,
/* DCS 2_7 */ 53,
/* DCS 2_8 */ 52,
/* DCS 2_9 */ 51,
/* DCS 3_0 */ 40,
/* DCS 3_1 */ 39,
/* DCS 3_2 */ 38,
/* DCS 3_3 */ 37,
/* DCS 3_4 */ 36,
/* DCS 3_5 */ 35,
/* DCS 3_6 */ 34,
/* DCS 3_7 */ 33,
/* DCS 3_8 */ 32,
/* DCS 3_9 */ 31,
/* DCS 4_0 */ 30,
/* DCS 4_1 */ 29,
/* DCS 4_2 */ 28,
/* DCS 4_3 */ 27,
/* DCS 4_4 */ 26,
/* DCS 4_5 */ 25,
/* DCS 4_6 */ 24,
/* DCS 4_7 */ 23,
/* DCS 4_8 */ 22,
/* DCS 4_9 */ 21,
/* DCS 5_0 ITS-3*/ 1,
/* DCS 5_1 */ 49,
/* DCS 5_2 */ 48,
/* DCS 5_3 */ 47,
/* DCS 5_4 */ 46,
/* DCS 5_5 */ 45,
/* DCS 5_6 */ 44,
/* DCS 5_7 */ 43,
/* DCS 5_8 */ 50,
/* DCS 5_9 ITS-4*/ 41
};



ClassImp(AliACORDEConstants)

//_____________________________________________________________________________
AliACORDEConstants::AliACORDEConstants()
  : TObject()
{
  // Default constructor
}


//_____________________________________________________________________________
AliACORDEConstants* AliACORDEConstants::Instance()
{
  if ( !fgInstance ) {
    fgInstance = new AliACORDEConstants;
  }
  return fgInstance;
}

//_____________________________________________________________________________
AliACORDEConstants::~AliACORDEConstants()
{
  fgInstance = 0;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::ModulePositionX(Int_t i) const
{
  // Module lenght
  return fgkModulePositionX[i];
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::ModulePositionY(Int_t i) const
{
  // Module lenght
  return fgkModulePositionY[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::ModulePositionZ(Int_t i) const
{
  // Module lenght
  return fgkModulePositionZ[i];
}


//_____________________________________________________________________________
Float_t AliACORDEConstants::SupportModulePositionX(Int_t i) const
{
  // Module lenght
  return fgkSupportModulePositionX[i];
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::SupportModulePositionY(Int_t i) const
{
  // Module lenght
  return fgkSupportModulePositionY[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::SupportModulePositionZ(Int_t i) const
{
  // Module lenght
  return fgkSupportModulePositionZ[i];
}



Float_t AliACORDEConstants::ExtraModulePositionX() const
{
  // Module lenght
  return fgkExtraModulePositionX;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::ExtraModulePositionY() const
{
  // Module lenght
  return fgkExtraModulePositionY;
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::ExtraModulePositionZ(Int_t i) const
{
  // Module lenght
  return fgkExtraModulePositionZ[i];
}
//_____________________________________________________________________________
Int_t AliACORDEConstants::ModuleElectronicChannel(Int_t i) const
{
	// return de ID (electronic channel in ACORDE) of each module
	// acording to the current match between DCS and Electronic nomenclature
	return fgkModuleElectronicChannel[i];
}
//_____________________________________________________________________________
Float_t AliACORDEConstants::ModuleLength() const
{
  // Module lenght
  return fgkModuleLength;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::ModuleWidth() const
{
  // Module width
  return fgkModuleWidth;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::ModuleHeight() const
{
  // Module height
  return fgkModuleHeight;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::PlasticLength() const
{
  // Length of the scintillator active zone for a single counter
  return fgkPlasticLength;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::PlasticWidth() const
{
  // Width of the scintillator active zone for a single counter
  return fgkPlasticWidth;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::PlasticHeight() const
{
  // Height of the scintillator active zone for a single counter
  return fgkPlasticHeight;
}

Float_t AliACORDEConstants::ProfileWidth() const
{
  // Width of the profile of the Al box
  return fgkProfileWidth;
}

Float_t AliACORDEConstants::ProfileThickness() const
{
  // Thickness of the profile of the Al box
  return fgkProfileThickness;
}


//_____________________________________________________________________________
Float_t AliACORDEConstants::Depth() const
{
  // Alice IP depth
  return fgkDepth;
}
