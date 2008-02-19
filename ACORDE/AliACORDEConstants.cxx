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
// Authors: Arturo Fernandez, Enrique Gamez, Mario Rodríguez Cahuantzi, Eleazar Cuautle(ICN-UNAM) 
//         FCFM-UAP, Mexico.
//
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
  582, 582, 582, 582, 582, 582, 582, 582, 582, 582,
  797, 797, 797, 797, 797, 797, 797, 797, 797, 797,
  850, 850, 850, 850, 850, 850, 850, 850, 850, 850,
  850, 850, 850, 850, 850, 850, 850, 850, 850, 850,
  797, 797, 797, 797, 797, 797, 797, 797, 797, 797,
  582, 582, 582, 582, 582, 609, 605, 582, 582, 582};
const Float_t AliACORDEConstants::fgkModulePositionZ[60] = {
  450, 350, 250, 150, 50, -50, -120, -280, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 150, 50, -50, -150, -250, -350, -450,
  450, 350, 250, 104, 50, -76, -176, -250, -350, -450};

const Float_t AliACORDEConstants::fgkExtraModulePositionZ[4] = {93.0, 18., -18, -93};
const Float_t AliACORDEConstants::fgkExtraModulePositionX = 0.0;
const Float_t AliACORDEConstants::fgkExtraModulePositionY = 850.0;

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
